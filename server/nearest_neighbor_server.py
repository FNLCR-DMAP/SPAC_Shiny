"""
Nearest neighbor visualization server module for SPAC Shiny application.
"""

from shiny import ui, render, reactive
import io
import tempfile
import multiprocessing
import matplotlib.pyplot as plt
from utils.plot_manager import PlotManager


def run_nn_worker(queue, adata, source_label, target_labels, image_id,
                annotation, plot_method, plot_type, log_scale, facet_plot,
                x_axis_rotation, shared_x_title, x_title_fontsize,
                color_mapping, figure_width, figure_height, figure_dpi,
                font_size):
    try:
        plt.clf()
        plt.close('all')

        from utils.template_wrapper import (
            register_memory_object,
            unregister_memory_object
        )
        from spac.templates.visualize_nearest_neighbor_template import (
            run_from_json
        )

        virtual_path = register_memory_object(adata)

        params = {
            "Upstream_Analysis": virtual_path,
            "Annotation": annotation,
            "Source_Anchor_Cell_Label": source_label,
            "Target_Cell_Label": (
                ",".join(target_labels) if target_labels else "All"
            ),
            "ImageID": image_id or "None",
            "Plot_Method": plot_method,
            "Plot_Type": plot_type,
            "Nearest_Neighbor_Associated_Table": "spatial_distance",
            "Log_Scale": log_scale,
            "Facet_Plot": facet_plot,
            "X_Axis_Label_Rotation": x_axis_rotation,
            "Shared_X_Axis_Title_": shared_x_title,
            "X_Axis_Title_Font_Size": x_title_fontsize if x_title_fontsize else "None",
            "Defined_Color_Mapping": color_mapping or "None",
            "Figure_Width": figure_width,
            "Figure_Height": figure_height,
            "Figure_DPI": figure_dpi,
            "Font_Size": font_size
        }

        try:
            figs, df_data = run_from_json(
                json_path=params,
                save_results=False,
                show_plot=False
            )
        finally:
            unregister_memory_object(virtual_path)

        fig = figs[0] if isinstance(figs, list) and len(figs) > 0 else figs
        if fig is None:
            queue.put("Error: No figure generated")
            return

        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=figure_dpi,
                    bbox_inches='tight', pad_inches=0.2)
        buf.seek(0)
        img_bytes = buf.read()
        plt.close(fig)

        queue.put((img_bytes, df_data))

    except Exception as e:
        import traceback
        traceback.print_exc()
        queue.put(f"Error: {str(e)}")


def nearest_neighbor_server(input, output, session, shared):
    pm = PlotManager('nn', shared, plot_type='process', data_key='df_nn')

    @reactive.calc
    def get_adata():
        return shared['adata_main'].get()

    @reactive.calc
    def process_target_labels():
        target_labels = input.nn_target_label()
        return list(target_labels) if target_labels and len(target_labels) > 0 else None

    @reactive.calc
    def get_plot_type():
        method = input.nn_plot_method()
        return input.nn_plot_type_numeric() if method == "numeric" else input.nn_plot_type_distribution()

    @reactive.calc
    def get_image_id():
        image_id = input.nn_image_id()
        return None if image_id == "None" else image_id

    @reactive.calc
    def get_color_mapping():
        try:
            color_mapping = input.nn_color_mapping()
            return None if color_mapping == "None" else color_mapping
        except Exception:
            return None

    @reactive.calc
    def detect_annotation():
        """Auto-detect annotation column matching spatial_distance phenotypes."""
        adata = get_adata()
        if adata is None:
            return None

        spatial_distance_key = "spatial_distance"
        distance_df = None
        if spatial_distance_key in adata.obsm:
            distance_df = adata.obsm[spatial_distance_key]
        elif spatial_distance_key in adata.uns:
            distance_df = adata.uns[spatial_distance_key]

        if distance_df is not None and hasattr(distance_df, 'columns'):
            spatial_phenotypes = set(distance_df.columns)
            for col in adata.obs.columns:
                is_categorical = (adata.obs[col].dtype == 'object' or
                                adata.obs[col].dtype.name == 'category')
                if is_categorical:
                    obs_phenotypes = set(adata.obs[col].unique())
                    overlap = spatial_phenotypes.intersection(obs_phenotypes)
                    if len(overlap) >= len(spatial_phenotypes) * 0.8:
                        return col

            # Fallback: first categorical column
            for col in adata.obs.columns:
                if (adata.obs[col].dtype == 'object' or
                        adata.obs[col].dtype.name == 'category'):
                    return col

        return None

    @output
    @render.ui
    def nn_color_mapping_ui():
        adata = get_adata()
        choices = {"None": "None (Auto)"}
        if adata is not None:
            if hasattr(adata, 'uns') and adata.uns is not None:
                for key in adata.uns.keys():
                    if key.endswith('_color_map') or 'color' in key.lower():
                        choices[key] = key
        return ui.input_select(
            "nn_color_mapping",
            ui.tags.span(
                "Defined Color Mapping",
                ui.tags.span(
                    "\u24D8",
                    title=(
                        "Color map from loaded data. "
                        "Use 'None (Auto)' for automatic coloring."
                    ),
                    tabindex="0",
                    class_="accessible-tooltip",
                    style="margin-left:5px; cursor:help; color:#007bff;"
                )
            ),
            choices=choices,
            selected="None"
        )

    @reactive.Effect
    @reactive.event(input.go_nn_viz, ignore_none=True)
    def start_nn_task():
        if pm.is_calculating.get():
            return

        adata = get_adata()
        if adata is None:
            return

        source_label = input.nn_source_label()
        if not source_label:
            return

        annotation = detect_annotation()
        if not annotation:
            return

        target_labels = process_target_labels()
        image_id = get_image_id()
        plot_method = input.nn_plot_method()
        plot_type = get_plot_type()
        log_scale = input.nn_log_scale()
        facet_plot = input.nn_facet_plot()
        x_axis_rotation = input.nn_x_axis_rotation()
        shared_x_title = input.nn_shared_x_title()
        x_title_fontsize = input.nn_x_title_fontsize()
        color_mapping = get_color_mapping()
        figure_width = input.nn_figure_width()
        figure_height = input.nn_figure_height()
        figure_dpi = input.nn_figure_dpi()
        font_size = input.nn_font_size()

        pm.start_process(
            run_nn_worker,
            args=(adata, source_label, target_labels, image_id,
                annotation, plot_method, plot_type, log_scale, facet_plot,
                x_axis_rotation, shared_x_title, x_title_fontsize,
                color_mapping, figure_width, figure_height, figure_dpi,
                font_size)
        )

    @reactive.Effect
    def check_status():
        def on_result(res):
            img_bytes, df_data = res
            pm.result.set(img_bytes)
            shared['df_nn'].set(df_data)
        pm.check_process(on_result)

    @output
    @render.image
    def nn_visualization_plot():
        img_bytes = pm.result.get()
        if img_bytes is None:
            return None
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            tmp.write(img_bytes)
            tmp_path = tmp.name
        return {
            "src": tmp_path,
            "contentType": "image/png",
            "style": "max-width: 100%; height: auto;"
        }

    @render.ui
    def nn_stop_button_ui():
        return pm.stop_button_ui('stop_nn')

    @render.download(filename="nearest_neighbor_data.csv")
    def download_df_nn():
        df = shared['df_nn'].get()
        if df is not None:
            return df.to_csv(index=False).encode("utf-8"), "text/csv"
        return None

    @render.ui
    def download_button_ui_nn():
        return pm.download_button_ui('download_df_nn')

    @render.ui
    def download_nn_plot_button_ui():
        return pm.plot_download_button_ui('download_nn_plot')

    @render.download(filename="nearest_neighbor_plot.png")
    def download_nn_plot():
        return pm.create_plot_download_handler()()