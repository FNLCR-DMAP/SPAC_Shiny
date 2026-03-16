from shiny import ui, render, reactive
import anndata as ad
import pandas as pd
import io
import tempfile
import multiprocessing
import matplotlib.pyplot as plt
import spac.visualization
from utils.plot_manager import PlotManager


def run_scatterplot_worker(queue, x, y, color_values, x_label, y_label, title):
    try:
        plt.clf()
        plt.close('all')

        if color_values is not None:
            fig, ax = spac.visualization.visualize_2D_scatter(
                x, y, labels=color_values
            )
            for a in fig.axes:
                if hasattr(a, "get_ylabel") and a != ax:
                    a.set_ylabel(f"Colored by: {y_label}")
        else:
            fig, ax = spac.visualization.visualize_2D_scatter(x, y)

        ax.set_title(title, fontsize=14)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)

        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=140,
                    bbox_inches='tight', pad_inches=0.2)
        buf.seek(0)
        img_bytes = buf.read()
        plt.close(fig)

        queue.put(img_bytes)

    except Exception as e:
        queue.put(f"Error: {str(e)}")


def scatterplot_server(input, output, session, shared):
    pm = PlotManager('scatterplot', shared, plot_type='process', data_key=None)

    @reactive.calc
    def get_scatterplot_names():
        obsm_list = shared['obsm_names'].get()
        var_list = shared['var_names'].get()
        if obsm_list is None or var_list is None:
            return []
        return {
            "Annotated Tables": {item: item for item in obsm_list},
            "Features": {item: item for item in var_list},
        }

    @reactive.Effect
    def update_select_input_layer_scatter():
        choices = get_scatterplot_names()
        ui.update_select("scatter_x", choices=choices)
        ui.update_select("scatter_y", choices=choices)

    @reactive.calc
    def get_scatterplot_coordinates_x():
        adata = ad.AnnData(
            X=shared['X_data'].get(),
            var=pd.DataFrame(shared['var_data'].get()),
            obsm=shared['obsm_data'].get(),
            layers=shared['layers_data'].get()
        )
        obsm_names = shared['obsm_names'].get()
        features = shared['var_names'].get()
        layer_selection = input.scatter_layer()
        selection = input.scatter_x()

        if selection in obsm_names:
            return adata.obsm[selection][:, 0]
        if selection in features:
            col_idx = adata.var_names.get_loc(selection)
            if layer_selection == "Original":
                return adata.X[:, col_idx]
            else:
                return adata.layers[layer_selection][:, col_idx]
        return None

    @reactive.calc
    def get_scatterplot_coordinates_y():
        adata = shared['adata_main'].get()
        obsm_names = shared['obsm_names'].get()
        features = shared['var_names'].get()
        layer_selection = input.scatter_layer()
        selection = input.scatter_y()

        if selection in obsm_names:
            return adata.obsm[selection][:, 1]
        if selection in features:
            col_idx = adata.var_names.get_loc(selection)
            if layer_selection == "Original":
                return adata.X[:, col_idx]
            else:
                return adata.layers[layer_selection][:, col_idx]
        return None

    scatter_ui_initialized = reactive.Value(False)

    @reactive.effect
    def scatter_reactivity():
        btn = input.scatter_color_check()
        if btn and not scatter_ui_initialized.get():
            ui.insert_ui(
                ui.div(
                    {"id": "inserted-scatter_dropdown"},
                    ui.input_select(
                        "scatter_color",
                        "Select Feature",
                        choices=shared['var_names'].get()
                    )
                ),
                selector="#main-scatter_dropdown",
                where="beforeEnd",
            )
            scatter_ui_initialized.set(True)
        elif not btn and scatter_ui_initialized.get():
            ui.remove_ui("#inserted-scatter_dropdown")
            scatter_ui_initialized.set(False)

    @reactive.calc
    def get_color_values():
        selected_feature = input.scatter_color()
        if selected_feature is None:
            return None
        adata = ad.AnnData(
            X=shared['X_data'].get(),
            var=pd.DataFrame(shared['var_data'].get())
        )
        if selected_feature in adata.var_names:
            col_idx = adata.var_names.get_loc(selected_feature)
            return adata.X[:, col_idx]
        return None

    @reactive.Effect
    @reactive.event(input.go_scatter, ignore_none=True)
    def start_scatterplot_task():
        if pm.is_calculating.get():
            return

        # Capture all reactive values before starting process
        x = get_scatterplot_coordinates_x()
        y = get_scatterplot_coordinates_y()
        if x is None or y is None:
            return

        color_enabled = input.scatter_color_check()
        color_values = get_color_values() if color_enabled else None
        x_label = input.scatter_x()
        y_label = input.scatter_y()
        title = f"Scatterplot: {x_label} vs {y_label}"

        pm.start_process(
            run_scatterplot_worker,
            args=(x, y, color_values, x_label, y_label, title)
        )

    @reactive.Effect
    def check_status():
        def on_result(res):
            pm.result.set(res)
        pm.check_process(on_result)

    @output
    @render.image
    def spac_Scatter():
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
    def scatterplot_stop_button_ui():
        return pm.stop_button_ui('stop_scatterplot')

    @render.ui
    def download_scatter_plot_button_ui():
        return pm.plot_download_button_ui('download_scatter_plot')

    @render.download(filename="scatterplot.png")
    def download_scatter_plot():
        return pm.create_plot_download_handler()()