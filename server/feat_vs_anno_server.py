"""
Feature vs Annotation heatmap visualization module for SPAC Shiny application.

This module handles the server-side logic for generating heatmaps that
visualize features (genes/proteins) against cell annotations using the
hierarchical_heatmap function.
"""

from shiny import ui, render, reactive, req
import anndata as ad
import numpy as np
import pandas as pd
import logging
import spac.visualization
from utils.plot_utils import abbreviate_labels, apply_axis_style
from utils.download_naming import build_download_filename


# Set up logger
logger = logging.getLogger(__name__)


def feat_vs_anno_server(input, output, session, shared):
    """
    Server logic for feature vs annotation heatmap visualization.

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

    def on_layer_check():
        """
        Get the selected layer name or None for original data.

        Returns
        -------
        str or None
            Layer name if not "Original", otherwise None
        """
        return input.hm1_layer() if input.hm1_layer() != "Original" else None

    def on_dendro_check():
        """
        Check if dendrogram is enabled and return the appropriate values.

        Returns
        -------
        tuple of (bool, bool) or (None, None)
            Annotation dendrogram and feature dendrogram flags.
            Returns (None, None) if dendrogram is disabled.
        """
        return (
            (input.hm1_anno_dendro(), input.hm1_feat_dendro())
            if input.hm1_dendogram()
            else (None, None)
        )

    @reactive.calc
    def get_adata():
        """
        Get the main AnnData object from shared state.

        Returns
        -------
        anndata.AnnData or None
            AnnData object reconstructed from shared data components
            Returns None if data is not loaded
        """
        x_data = shared['X_data'].get()

        # STOP THE CRASH: If data isn't loaded, don't access .dtype
        if x_data is None:
            return None

        return ad.AnnData(
            X=x_data,
            obs=pd.DataFrame(shared['obs_data'].get()),
            var=pd.DataFrame(shared['var_data'].get()),
            layers=shared['layers_data'].get(),
            dtype=x_data.dtype
        )

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

        vmin = input.hm1_min_select()
        vmax = input.hm1_max_select()
        kwargs = {"vmin": vmin, "vmax": vmax}
        cluster_annotations, cluster_features = on_dendro_check()

        # Error Handling: Catch and log specific errors
        try:
            df, fig, ax = spac.visualization.hierarchical_heatmap(
                adata,
                annotation=input.hm1_anno(),
                layer=on_layer_check(),
                z_score=None,
                cluster_annotations=cluster_annotations,
                cluster_feature=cluster_features,
                **kwargs
            )
        except ValueError as e:
            error_msg = ("Heatmap generation failed with invalid "
                        f"parameters: {e}")
            logger.error(error_msg)
            return None
        except Exception as e:
            error_msg = ("Unexpected error during heatmap "
                        f"generation: {e}")
            logger.error(error_msg)
            return None

        if fig is None or not hasattr(fig, "ax_heatmap"):
            logger.error("Invalid figure structure.")
            return None

        # Apply colormap
        cmap = input.hm1_cmap()
        if cmap != "viridis":
            fig.ax_heatmap.collections[0].set_cmap(cmap)

        shared['df_heatmap'].set(df)

        # Rotate X and Y axis labels
        fig.ax_heatmap.set_xticklabels(
            fig.ax_heatmap.get_xticklabels(),
            rotation=input.hm1_x_label_rotation(),
            horizontalalignment='right'
        )
        fig.ax_heatmap.set_yticklabels(
            fig.ax_heatmap.get_yticklabels(),
            rotation=input.hm1_y_label_rotation(),
            verticalalignment='center'
        )

        # Abbreviate labels if enabled
        if input.hm1_enable_abbreviation():
            limit = input.hm1_label_char_limit()
            abbreviated_xticks = abbreviate_labels(
                fig.ax_heatmap.get_xticklabels(), limit)
            fig.ax_heatmap.set_xticklabels(
                abbreviated_xticks, rotation=input.hm1_x_label_rotation())
            abbreviated_yticks = abbreviate_labels(
                fig.ax_heatmap.get_yticklabels(), limit)
            fig.ax_heatmap.set_yticklabels(
                abbreviated_yticks, rotation=input.hm1_y_label_rotation())

        # Set font size for axis labels
        axis_fontsize = input.hm1_axis_label_fontsize()
        apply_axis_style(fig.ax_heatmap.get_xticklabels(), axis_fontsize)
        apply_axis_style(fig.ax_heatmap.get_yticklabels(), axis_fontsize)

        # Adjust figure layout with small margins to prevent label clipping
        # rect format: [left, bottom, right, top] as fraction of figure size
        LAYOUT_RECT = (0.02, 0.02, 0.98, 0.98)
        fig.fig.tight_layout(rect=LAYOUT_RECT)
        fig.fig.subplots_adjust(bottom=0.15, left=0)
        return fig

    def get_heatmap_filename():
        return build_download_filename(shared, "heatmap", mime_type="text/csv")

    @render.download(filename=get_heatmap_filename)
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
                # Check if layer exists in AnnData
                if input.hm1_layer() not in adata.layers:
                    return None
                layer_data = adata.layers[input.hm1_layer()]

            # Check if annotation exists in obs
            if input.hm1_anno() not in adata.obs:
                return None

            # Filter layer data based on valid annotations
            mask = adata.obs[input.hm1_anno()].notna()
            layer_data = layer_data[mask]

            # Avoid calculation on empty data
            if layer_data.size == 0:
                return None

            min_val = round(float(np.min(layer_data)), 2)
            max_val = round(float(np.max(layer_data)), 2)

            # UI Update
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
            logger.error(f"Error updating min/max values: {e}")
            return None
