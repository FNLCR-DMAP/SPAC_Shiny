"""
Feature vs Annotation heatmap visualization module for SPAC Shiny application.

This module handles the server-side logic for generating heatmaps that
visualize features (genes/proteins) against cell annotations using the
hierarchical_heatmap_template functionality.
"""

from shiny import ui, render, reactive, req
import numpy as np
import logging
from typing import Tuple, Any

from utils.template_wrapper import (
    register_memory_object,
    unregister_memory_object,
)
from utils.plot_utils import abbreviate_labels, apply_axis_style
from spac.templates.hierarchical_heatmap_template import run_from_json


# Set up logger
logger = logging.getLogger(__name__)


def feat_vs_anno_server(input, output, session, shared):
    """
    Server logic for feature vs annotation heatmap visualization.

    This implementation registers the in-memory AnnData object with the
    memory registry and delegates plotting to
    `spac.templates.hierarchical_heatmap_template.run_from_json`.

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

    @render.ui
    def hm1_features_ui():
        """Create dynamic feature multi-select from loaded data."""
        adata = get_adata()
        if adata is None:
            return None

        feature_choices = ["All"] + sorted(adata.var_names.tolist())
        return ui.input_selectize(
            "hm1_features",
            "Feature(s)",
            choices=feature_choices,
            selected=["All"],
            multiple=True
        )

    @reactive.calc
    def get_layer():
        """
        Get the selected layer name.

        Returns
        -------
        str
            Layer name, or 'Original' to use adata.X.
        """
        return input.hm1_layer() if input.hm1_layer() != "Original" else "Original"

    @reactive.calc
    def get_dendrogram_settings():
        """
        Check if dendrogram is enabled and return the appropriate values.

        Returns
        -------
        tuple of (bool, bool)
            (cluster_annotations, cluster_features)
        """
        if input.hm1_dendogram():
            return (input.hm1_anno_dendro(), input.hm1_feat_dendro())
        return (False, False)

    @reactive.calc
    def get_vmin_vmax():
        """
        Process min/max value inputs.

        Returns
        -------
        tuple of (str or float, str or float)
            (vmin, vmax) values, or "None" if not available
        """
        try:
            vmin = input.hm1_min_select()
            vmax = input.hm1_max_select()
        except Exception:
            vmin = "None"
            vmax = "None"
        return (vmin if vmin is not None else "None",
                vmax if vmax is not None else "None")

    @output
    @render.plot(alt="Heatmap Plot")
    @reactive.event(input.go_hm1, ignore_none=True)
    def spac_Heatmap():
        """
        Render heatmap of features vs annotations.

        This function generates a clustered heatmap showing the relationship
        between selected features (columns) and cell annotations (rows).

        Returns
        -------
        matplotlib.figure.Figure or None
            Heatmap figure with optional dendrograms, or None if
            generation fails
        """
        # Validation: Ensure required inputs are present
        req(input.hm1_anno())
        req(input.hm1_layer())

        adata = get_adata()
        if adata is None:
            return None

        annotation = input.hm1_anno()
        layer = get_layer()
        cmap = input.hm1_cmap()
        cluster_annotations, cluster_features = get_dendrogram_settings()
        vmin, vmax = get_vmin_vmax()

        # Register adata in memory registry and call run_from_json
        try:
            virtual_path = register_memory_object(adata)

            params = {
                # Core parameters
                "Upstream_Analysis": virtual_path,
                "Annotation": annotation,
                "Table_to_Visualize": layer,
                "Feature_s_": list(input.hm1_features()) if input.hm1_features() else ["All"],
                "Feature_Dendrogram": cluster_features,
                "Annotation_Dendrogram": cluster_annotations,

                # Plot configuration — read from UI inputs
                "Z_Score": input.hm1_z_score(),
                "Standard_Scale_": input.hm1_standard_scale(),
                "Swap_Axes": input.hm1_swap_axes(),
                "Color_Map": cmap,
                "Value_Min": vmin,
                "Value_Max": vmax,

                # Figure configuration — read from UI inputs
                "Figure_Title": input.hm1_figure_title(),
                "Figure_Width": input.hm1_figure_width(),
                "Figure_Height": input.hm1_figure_height(),
                "Figure_DPI": input.hm1_figure_dpi(),
                "Font_Size": input.hm1_font_size(),
                "Matrix_Plot_Ratio": input.hm1_matrix_ratio(),

                # Dendrogram ratios — keep defaults
                "Rotate_Label_": False,
                "Horizontal_Dendrogram_Display_Ratio": 0.2,
                "Vertical_Dendrogram_Display_Ratio": 0.2,
            }

            # Call template to get ClusterGrid and dataframe in-memory.
            # Requires the fixed hierarchical_heatmap_template that
            # returns (clustergrid, mean_intensity) when
            # save_results_flag=False. See SCSAWorkflow PR #425.
            result: Tuple[Any, Any] = run_from_json(
                json_path=params,
                save_results_flag=False,
                show_plot=False
            )
            if result is None:
                return None

            clustergrid, df = result

        except ValueError as e:
            logger.error(
                "Heatmap generation failed with invalid parameters: %s", e
            )
            return None
        except Exception as e:
            logger.error(
                "Unexpected error during heatmap generation: %s", e
            )
            return None
        finally:
            try:
                unregister_memory_object(virtual_path)
            except Exception:
                pass

        if clustergrid is None or not hasattr(clustergrid, "ax_heatmap"):
            logger.error("Invalid figure structure.")
            return None

        # Apply colormap
        if cmap != "viridis":
            clustergrid.ax_heatmap.collections[0].set_cmap(cmap)

        shared['df_heatmap'].set(df)

        # Rotate X and Y axis labels
        clustergrid.ax_heatmap.set_xticklabels(
            clustergrid.ax_heatmap.get_xticklabels(),
            rotation=input.hm1_x_label_rotation(),
            horizontalalignment='right'
        )
        clustergrid.ax_heatmap.set_yticklabels(
            clustergrid.ax_heatmap.get_yticklabels(),
            rotation=input.hm1_y_label_rotation(),
            verticalalignment='center'
        )

        # Abbreviate labels if enabled
        if input.hm1_enable_abbreviation():
            limit = input.hm1_label_char_limit()
            abbreviated_xticks = abbreviate_labels(
                clustergrid.ax_heatmap.get_xticklabels(), limit)
            clustergrid.ax_heatmap.set_xticklabels(
                abbreviated_xticks,
                rotation=input.hm1_x_label_rotation())
            abbreviated_yticks = abbreviate_labels(
                clustergrid.ax_heatmap.get_yticklabels(), limit)
            clustergrid.ax_heatmap.set_yticklabels(
                abbreviated_yticks,
                rotation=input.hm1_y_label_rotation())

        # Set font size for axis labels
        axis_fontsize = input.hm1_axis_label_fontsize()
        apply_axis_style(
            clustergrid.ax_heatmap.get_xticklabels(), axis_fontsize)
        apply_axis_style(
            clustergrid.ax_heatmap.get_yticklabels(), axis_fontsize)

        # Adjust figure layout with small margins to prevent label clipping
        LAYOUT_RECT = (0.02, 0.02, 0.98, 0.98)
        clustergrid.fig.tight_layout(rect=LAYOUT_RECT)
        clustergrid.fig.subplots_adjust(bottom=0.15, left=0)
        return clustergrid

    @render.download(filename="heatmap_data.csv")
    def download_df_hm1():
        df = shared['df_heatmap'].get()
        if df is not None:
            csv_string = df.to_csv(index=False)
            csv_bytes = csv_string.encode("utf-8")
            return csv_bytes, "text/csv"
        return None

    @render.ui
    @reactive.event(input.go_hm1, ignore_none=True)
    def download_button_ui_hm1():
        if shared['df_heatmap'].get() is not None:
            return ui.download_button(
                "download_df_hm1", "Download Data", class_="btn-warning")
        return None

    @reactive.effect
    @reactive.event(input.hm1_layer)
    def update_min_max():
        req(input.hm1_anno())
        req(input.hm1_layer())

        adata = get_adata()
        if adata is None:
            return None

        try:
            # Determine layer data source
            if input.hm1_layer() == "Original":
                layer_data = adata.X
            else:
                if input.hm1_layer() not in adata.layers:
                    return None
                layer_data = adata.layers[input.hm1_layer()]

            if input.hm1_anno() not in adata.obs:
                return None

            mask = adata.obs[input.hm1_anno()].notna()
            layer_data = layer_data[mask]

            if layer_data.size == 0:
                return None

            min_val = round(float(np.min(layer_data)), 2)
            max_val = round(float(np.max(layer_data)), 2)

            ui.remove_ui("#inserted-hm1_min_num")
            ui.remove_ui("#inserted-hm1_max_num")

            min_num = ui.input_numeric(
                "hm1_min_select",
                "Minimum",
                min_val,
                min=min_val,
                max=max_val
            )
            ui.insert_ui(
                ui.div({"id": "inserted-hm1_min_num"}, min_num),
                selector="#main-hm1_min_num",
                where="beforeEnd",
            )

            max_num = ui.input_numeric(
                "hm1_max_select",
                "Maximum",
                max_val,
                min=min_val,
                max=max_val
            )
            ui.insert_ui(
                ui.div({"id": "inserted-hm1_max_num"}, max_num),
                selector="#main-hm1_max_num",
                where="beforeEnd",
            )
        except Exception as e:
            logger.error("Error updating min/max values: %s", e)
            return None
