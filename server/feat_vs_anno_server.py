from shiny import ui, render, reactive
import anndata as ad
import numpy as np
import pandas as pd
import io
import tempfile
import multiprocessing
import matplotlib.pyplot as plt
import spac.visualization
from utils.plot_manager import PlotManager


def run_heatmap_worker(queue, adata, annotation, layer, cluster_annotations,
                       cluster_features, vmin, vmax, cmap, x_rotation):
    try:
        plt.clf()
        plt.close('all')

        kwargs = {"vmin": vmin, "vmax": vmax}

        df, fig, ax = spac.visualization.hierarchical_heatmap(
            adata,
            annotation=annotation,
            layer=layer,
            z_score=None,
            cluster_annotations=cluster_annotations,
            cluster_feature=cluster_features,
            **kwargs
        )

        if cmap != "viridis":
            fig.ax_heatmap.collections[0].set_cmap(cmap)

        fig.ax_heatmap.set_xticklabels(
            fig.ax_heatmap.get_xticklabels(),
            rotation=x_rotation,
            horizontalalignment='right'
        )
        fig.fig.subplots_adjust(bottom=0.4)
        fig.fig.subplots_adjust(left=0.1)

        buf = io.BytesIO()
        fig.fig.savefig(buf, format='png', dpi=140, bbox_inches='tight', pad_inches=0.2)
        buf.seek(0)
        img_bytes = buf.read()
        plt.close(fig.fig)

        queue.put((img_bytes, df))
    except Exception as e:
        queue.put(f"Error: {str(e)}")


def feat_vs_anno_server(input, output, session, shared):
    pm = PlotManager('hm1', shared, plot_type='process', data_key='df_heatmap')

    def on_layer_check():
        return input.hm1_layer() if input.hm1_layer() != "Original" else None

    def on_dendro_check():
        return (
            (input.h2_anno_dendro(), input.h2_feat_dendro())
            if input.dendogram()
            else (None, None)
        )

    @reactive.Effect
    @reactive.event(input.go_hm1, ignore_none=True)
    def start_heatmap_task():
        if pm.is_calculating.get():
            return

        adata = ad.AnnData(
            X=shared['X_data'].get(),
            obs=pd.DataFrame(shared['obs_data'].get()),
            var=pd.DataFrame(shared['var_data'].get()),
            layers=shared['layers_data'].get(),
            dtype=shared['X_data'].get().dtype
        )
        if adata is None:
            return

        vmin = input.min_select()
        vmax = input.max_select()
        cmap = input.hm1_cmap()
        cluster_annotations, cluster_features = on_dendro_check()
        annotation = input.hm1_anno()
        layer = on_layer_check()
        x_rotation = input.hm_x_label_rotation()

        pm.start_process(
            run_heatmap_worker,
            args=(adata, annotation, layer, cluster_annotations,
                  cluster_features, vmin, vmax, cmap, x_rotation)
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
    def download_button_ui():
        return pm.download_button_ui('download_df')

    @render.download(filename="heatmap_data.csv")
    def download_df():
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

    heatmap_ui_initialized = reactive.Value(False)

    @reactive.effect
    def heatmap_reactivity():
        btn = input.dendogram()
        ui_initialized = heatmap_ui_initialized.get()

        if btn and not ui_initialized:
            ui.insert_ui(
                ui.div(
                    {"id": "inserted-check"},
                    ui.input_checkbox("h2_anno_dendro", "Annotation Cluster", value=False)
                ),
                selector="#main-hm1_check",
                where="beforeEnd",
            )
            ui.insert_ui(
                ui.div(
                    {"id": "inserted-check1"},
                    ui.input_checkbox("h2_feat_dendro", "Feature Cluster", value=False)
                ),
                selector="#main-hm2_check",
                where="beforeEnd",
            )
            heatmap_ui_initialized.set(True)

        elif not btn and ui_initialized:
            ui.remove_ui("#inserted-check")
            ui.remove_ui("#inserted-check1")
            heatmap_ui_initialized.set(False)

    @reactive.effect
    @reactive.event(input.hm1_layer)
    def update_min_max():
        adata = ad.AnnData(
            X=shared['X_data'].get(),
            obs=pd.DataFrame(shared['obs_data'].get()),
            var=pd.DataFrame(shared['var_data'].get()),
            layers=shared['layers_data'].get()
        )
        if input.hm1_layer() == "Original":
            layer_data = adata.X
        else:
            layer_data = adata.layers[input.hm1_layer()]
        mask = adata.obs[input.hm1_anno()].notna()
        layer_data = layer_data[mask]
        min_val = round(float(np.min(layer_data)), 2)
        max_val = round(float(np.max(layer_data)), 2)

        ui.remove_ui("#inserted-min_num")
        ui.remove_ui("#inserted-max_num")

        ui.insert_ui(
            ui.div(
                {"id": "inserted-min_num"},
                ui.input_numeric("min_select", "Minimum", min_val, min=min_val, max=max_val)
            ),
            selector="#main-min_num",
            where="beforeEnd",
        )
        ui.insert_ui(
            ui.div(
                {"id": "inserted-max_num"},
                ui.input_numeric("max_select", "Maximum", max_val, min=min_val, max=max_val)
            ),
            selector="#main-max_num",
            where="beforeEnd",
        )