"""
Nearest neighbor visualization server module for SPAC Shiny application.

This module handles the server-side logic for visualizing precomputed nearest
neighbor distances using the visualize_nearest_neighbor_template functionality.
"""

from shiny import ui, render, reactive, req
from utils.plot_utils import fig_to_png_bytes, png_bytes_to_figure


def nearest_neighbor_server(input, output, session, shared):
    """
    Server logic for nearest neighbor visualization feature.

    Parameters
    ----------
    input : shiny.session.Inputs
        Shiny input object
    output : shiny.session.Outputs
        Shiny output object
    session : shiny.session.Session
        Shiny session object
    shared : dict
        Shared reactive values across server modules
    """

    @reactive.calc
    def get_adata():
        """Get the main AnnData object from shared state."""
        return shared['adata_main'].get()

    @reactive.calc
    def process_target_labels():
        """
        Process target label selection.

        Returns
        -------
        list or None
            List of target phenotypes or None for 'All'
        """
        target_labels = input.nn_target_label()
        if target_labels and len(target_labels) > 0:
            return target_labels
        return None

    @reactive.calc
    def get_plot_type():
        """Get the appropriate plot type based on method selection."""
        method = input.nn_plot_method()
        if method == "numeric":
            return input.nn_plot_type_numeric()
        else:
            return input.nn_plot_type_distribution()

    @reactive.calc
    def get_image_id():
        """Process ImageID selection."""
        image_id = input.nn_image_id()
        return None if image_id == "None" else image_id

    @reactive.calc
    def get_color_mapping():
        """Process color mapping selection."""
        try:
            color_mapping = input.nn_color_mapping()
            return None if color_mapping == "None" else color_mapping
        except Exception:
            # Return None if input not available yet (dynamic UI not rendered)
            return None

    @reactive.calc
    def get_font_size():
        """Process font size, returning None if using default."""
        font_size = input.nn_x_title_fontsize()
        return font_size if font_size != 12 else None

    @output
    @render.ui
    def nn_color_mapping_ui():
        """
        Create dynamic color mapping select input from available data.

        Returns
        -------
        shiny.ui element
            Select input with available color mappings or None option
        """
        adata = get_adata()
        choices = {"None": "None (Auto)"}

        if adata is not None:
            # Extract available color mappings from uns
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

    @output
    @render.plot
    @reactive.event(input.go_nn_viz, ignore_none=True)
    def nn_visualization_plot():
        """
        Generate the nearest neighbor visualization plot.

        Returns
        -------
        matplotlib.figure.Figure
            The generated plot figure
        """
        adata = get_adata()
        if adata is None:
            return None

        source_label = input.nn_source_label()
        if not source_label:
            return None

        # Auto-detect annotation column matching spatial_distance phenotypes
        annotation = None
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
                        annotation = col
                        break

            if annotation is None:
                # Fallback: use the first categorical column
                for col in adata.obs.columns:
                    is_obj = adata.obs[col].dtype == 'object'
                    is_cat = adata.obs[col].dtype.name == 'category'
                    if is_obj or is_cat:
                        annotation = col
                        break

        if not annotation:
            return None

        cache = shared['cache']
        version = shared['dataset_version'].get()

        target_labels = process_target_labels()
        image_id = get_image_id()
        color_mapping = get_color_mapping()
        font_size_val = input.nn_x_title_fontsize()

        params = {
            "Annotation": annotation,
            "Source_Anchor_Cell_Label": source_label,
            "Target_Cell_Label": (
                tuple(target_labels) if target_labels else None
            ),
            "ImageID": image_id or "None",
            "Plot_Method": input.nn_plot_method(),
            "Plot_Type": get_plot_type(),
            "Nearest_Neighbor_Associated_Table": "spatial_distance",
            "Log_Scale": input.nn_log_scale(),
            "Facet_Plot": input.nn_facet_plot(),
            "X_Axis_Label_Rotation": input.nn_x_axis_rotation(),
            "Shared_X_Axis_Title_": input.nn_shared_x_title(),
            "X_Axis_Title_Font_Size": (
                font_size_val if font_size_val else "None"
            ),
            "Defined_Color_Mapping": color_mapping or "None",
            "Figure_Width": input.nn_figure_width(),
            "Figure_Height": input.nn_figure_height(),
            "Figure_DPI": input.nn_figure_dpi(),
            "Font_Size": input.nn_font_size(),
        }

        def compute():
            try:
                from utils.template_wrapper import (
                    register_memory_object,
                    unregister_memory_object
                )
                from spac.templates.visualize_nearest_neighbor_template import (
                    run_from_json
                )

                virtual_path = register_memory_object(adata)

                nn_params = {
                    "Upstream_Analysis": virtual_path,
                    **params,
                    "Target_Cell_Label": (
                        ",".join(params["Target_Cell_Label"])
                        if params["Target_Cell_Label"] else "All"
                    ),
                }

                try:
                    figs, df_data = run_from_json(
                        json_path=nn_params,
                        save_results=False,
                        show_plot=False
                    )
                finally:
                    unregister_memory_object(virtual_path)

                if isinstance(figs, list):
                    fig = figs[0] if len(figs) > 0 else None
                else:
                    fig = figs

                if fig is None:
                    return None, None

                return fig_to_png_bytes(fig), df_data

            except Exception:
                import traceback
                traceback.print_exc()
                return None, None

        img_bytes, df = cache.get_or_compute(
            'nearest_neighbor', version, params, compute
        )

        if img_bytes is None:
            return None

        shared['df_nn'].set(df)
        return png_bytes_to_figure(img_bytes)

    @render.download(filename="nearest_neighbor_data.csv")
    def download_df_nn():
        """
        Download the nearest neighbor data as CSV.

        Returns
        -------
        tuple
            CSV bytes and content type
        """
        df = shared['df_nn'].get()
        if df is not None:
            csv_string = df.to_csv(index=False)
            csv_bytes = csv_string.encode("utf-8")
            return csv_bytes, "text/csv"
        return None

    @render.ui
    @reactive.event(input.go_nn_viz, ignore_none=True)
    def download_button_ui_nn():
        """
        Show download button when data is available.

        Returns
        -------
        shiny.ui element or None
            Download button UI or None if no data
        """
        if shared['df_nn'].get() is not None:
            return ui.download_button(
                "download_df_nn",
                "Download Data",
                class_="btn-warning"
            )
        return None
