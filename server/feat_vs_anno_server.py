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
import io
import tempfile
import multiprocessing
import matplotlib.pyplot as plt
import logging
import spac.visualization
from utils.plot_manager import PlotManager
from utils.plot_utils import abbreviate_labels, apply_axis_style


logger = logging.getLogger(__name__)


def run_heatmap_worker(queue, adata, annotation, layer, cluster_annotations,
                       cluster_features, vmin, vmax, cmap,
                       x_rotation, y_rotation,
                       enable_abbrev, char_limit, axis_fontsize):
    try:
        plt.clf()
        plt.close('all')

        df, fig, ax = spac.visualization.hierarchical_heatmap(
            adata,
            annotation=annotation,
            layer=layer,
            z_score=None,
            cluster_annotations=cluster_annotations,
            cluster_feature=cluster_features,
            vmin=vmin,
            vmax=vmax,
        )

        if fig is None or not hasattr(fig, "ax_heatmap"):
            queue.put("Error: Invalid figure structure")
            return

        if cmap != "viridis":
            fig.ax_heatmap.collections[0].set_cmap(cmap)

        fig.ax_heatmap.set_xticklabels(
            fig.ax_heatmap.get_xticklabels(),
            rotation=x_rotation,
            horizontalalignment='right'
        )
        fig.ax_heatmap.set_yticklabels(
            fig.ax_heatmap.get_yticklabels(),
            rotation=y_rotation,
            verticalalignment='center'
        )

        if enable_abbrev and char_limit:
            abbreviated_xticks = abbreviate_labels(
                fig.ax_heatmap.get_xticklabels(), char_limit)
            fig.ax_heatmap.set_xticklabels(
                abbreviated_xticks, rotation=x_rotation)
            abbreviated_yticks = abbreviate_labels(
                fig.ax_heatmap.get_yticklabels(), char_limit)
            fig.ax_heatmap.set_yticklabels(
                abbreviated_yticks, rotation=y_rotation)

        apply_axis_style(fig.ax_heatmap.get_xticklabels(), axis_fontsize)
        apply_axis_style(fig.ax_heatmap.get_yticklabels(), axis_fontsize)

        LAYOUT_RECT = (0.02, 0.02, 0.98, 0.98)
        fig.fig.tight_layout(rect=LAYOUT_RECT)
        fig.fig.subplots_adjust(bottom=0.15, left=0)

        buf = io.BytesIO()
        fig.fig.savefig(buf, format='png', dpi=140,
                        bbox_inches='tight', pad_inches=0.2)
        buf.seek(0)
        img_bytes = buf.read()
        plt.close(fig.fig)

        queue.put((img_bytes, df))
    except Exception as e:
        queue.put(f"Error: {str(e)}")


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
    pm = PlotManager('hm1', shared, plot_type='process', data_key='df_heatmap')

    def on_layer_check():
        return input.hm1_layer() if input.hm1_layer() != "Original" else None

    def on_dendro_check():
        return (
            (input.hm1_anno_dendro(), input.hm1_feat_dendro())
            if input.hm1_dendogram()
            else (None, None)
        )

    @reactive.calc
    def get_adata():
        x_data = shared['X_data'].get()
        if x_data is None:
            return None
        return ad.AnnData(
            X=x_data,
            obs=pd.DataFrame(shared['obs_data'].get()),
            var=pd.DataFrame(shared['var_data'].get()),
            layers=shared['layers_data'].get(),
            dtype=x_data.dtype
        )

    @reactive.Effect
    @reactive.event(input.go_hm1, ignore_none=True)
    def start_heatmap_task():
        if pm.is_calculating.get():
            return

        req(input.hm1_anno())
        req(input.hm1_layer())

        adata = get_adata()
        if adata is None:
            return

        cluster_annotations, cluster_features = on_dendro_check()
        enable_abbrev = input.hm1_enable_abbreviation()

        pm.start_process(
            run_heatmap_worker,
            args=(
                adata,
                input.hm1_anno(),
                on_layer_check(),
                cluster_annotations,
                cluster_features,
                input.hm1_min_select(),
                input.hm1_max_select(),
                input.hm1_cmap(),
                input.hm1_x_label_rotation(),
                input.hm1_y_label_rotation(),
                enable_abbrev,
                input.hm1_label_char_limit() if enable_abbrev else None,
                input.hm1_axis_label_fontsize(),
            )
        )

    @reactive.Effect
    def check_status():
        def on_result(res):
            img_bytes, df = res
            pm.result.set(img_bytes)
            shared['df_heatmap'].set(df)
        pm.check_process(on_result)

    @output
    @render.image
    def spac_Heatmap():
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
    def heatmap_stop_button_ui():
        return pm.stop_button_ui('stop_hm1')

    @render.ui
    def download_button_ui_hm1():
        return pm.download_button_ui('download_df_hm1')

    @render.download(filename="heatmap_data.csv")
    def download_df_hm1():
        df = shared['df_heatmap'].get()
        if df is not None:
            return df.to_csv(index=False).encode("utf-8"), "text/csv"
        return None

    @render.ui
    def download_heatmap_plot_button_ui():
        return pm.plot_download_button_ui('download_heatmap_plot')

    @render.download(filename="heatmap_plot.png")
    def download_heatmap_plot():
        return pm.create_plot_download_handler()()

    @reactive.effect
    @reactive.event(input.hm1_layer)
    def update_min_max():
        req(input.hm1_anno())
        req(input.hm1_layer())

        adata = get_adata()
        if adata is None:
            return None

        try:
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

            ui.insert_ui(
                ui.div(
                    {"id": "inserted-hm1_min_num"},
                    ui.input_numeric(
                        "hm1_min_select", "Minimum",
                        min_val, min=min_val, max=max_val
                    )
                ),
                selector="#main-hm1_min_num",
                where="beforeEnd",
            )
            ui.insert_ui(
                ui.div(
                    {"id": "inserted-hm1_max_num"},
                    ui.input_numeric(
                        "hm1_max_select", "Maximum",
                        max_val, min=min_val, max=max_val
                    )
                ),
                selector="#main-hm1_max_num",
                where="beforeEnd",
            )
        except Exception as e:
            logger.error(f"Error updating min/max values: {e}")
            return None
