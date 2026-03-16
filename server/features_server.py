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


def run_features_histogram_worker(queue, adata, feature, layer, x_log_scale,
                                   y_log_scale, group_by, together, multiple, rotation):
    try:
        plt.clf()
        plt.close('all')

        kwargs = {
            "adata": adata,
            "feature": feature,
            "layer": layer,
            "x_log_scale": x_log_scale,
            "y_log_scale": y_log_scale,
        }
        if group_by:
            kwargs["group_by"] = group_by
            kwargs["together"] = together
            if multiple:
                kwargs["multiple"] = multiple

        fig, ax, df = spac.visualization.histogram(**kwargs).values()

        axes = ax if isinstance(ax, (list, np.ndarray)) else [ax]
        for a in axes:
            a.tick_params(axis='x', rotation=rotation, labelsize=10)

        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=140, bbox_inches=None, pad_inches=0.2)
        buf.seek(0)
        img_bytes = buf.read()
        plt.close(fig)

        queue.put((img_bytes, df))
    except Exception as e:
        queue.put(f"Error: {str(e)}")


def features_server(input, output, session, shared):
    pm = PlotManager('features', shared, plot_type='process', data_key='df_histogram1')

    def on_layer_check():
        return input.h1_layer() if input.h1_layer() != "Original" else None

    @reactive.Effect
    @reactive.event(input.go_h1, ignore_none=True)
    def start_features_task():
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

        feature = input.h1_feat()
        layer = on_layer_check()
        x_log_scale = input.h1_log_x()
        y_log_scale = input.h1_log_y()
        rotation = input.feat_slider()
        is_grouped = input.h1_group_by_check()
        group_by = input.h1_anno() if is_grouped else None
        together = input.h1_together_check() if is_grouped else False
        multiple = input.h1_together_drop() if (is_grouped and input.h1_together_check()) else None

        pm.start_process(
            run_features_histogram_worker,
            args=(adata, feature, layer, x_log_scale, y_log_scale,
                  group_by, together, multiple, rotation)
        )

    @reactive.Effect
    def check_status():
        def on_result(res):
            img_bytes, df = res
            pm.result.set(img_bytes)
            shared['df_histogram1'].set(df)
        pm.check_process(on_result)

    @output
    @render.image
    def spac_Histogram_1():
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
    def features_stop_button_ui():
        return pm.stop_button_ui('stop_features')

    @render.ui
    def download_histogram1_button_ui():
        return pm.download_button_ui('download_histogram1_df')

    @render.download(filename="features_histogram_data.csv")
    def download_histogram1_df():
        df = shared['df_histogram1'].get()
        if df is not None:
            return df.to_csv(index=False).encode("utf-8"), "text/csv"
        return None
    @render.ui
    def download_features_plot_button_ui():
        return pm.plot_download_button_ui('download_histogram1_plot')

    @render.download(filename="features_plot.png")
    def download_histogram1_plot():
        return pm.create_plot_download_handler()()

    histogram_ui_initialized = reactive.Value(False)

    @reactive.effect
    def histogram_reactivity():
        btn = input.h1_group_by_check()
        ui_initialized = histogram_ui_initialized.get()

        if btn and not ui_initialized:
            ui.insert_ui(
                ui.div(
                    {"id": "inserted-dropdown"},
                    ui.input_select(
                        "h1_anno",
                        "Select an Annotation",
                        choices=shared['obs_names'].get()
                    )
                ),
                selector="#main-h1_dropdown",
                where="beforeEnd",
            )
            ui.insert_ui(
                ui.div(
                    {"id": "inserted-check"},
                    ui.input_checkbox("h1_together_check", "Plot Together", value=True)
                ),
                selector="#main-h1_check",
                where="beforeEnd",
            )
            histogram_ui_initialized.set(True)

        elif not btn and ui_initialized:
            ui.remove_ui("#inserted-dropdown")
            ui.remove_ui("#inserted-check")
            ui.remove_ui("#inserted-dropdown_together")
            histogram_ui_initialized.set(False)

    @reactive.effect
    @reactive.event(input.h1_together_check)
    def update_stack_type_dropdown():
        if input.h1_together_check():
            ui.insert_ui(
                ui.div(
                    {"id": "inserted-dropdown_together"},
                    ui.input_select(
                        "h1_together_drop",
                        "Select Stack Type",
                        choices=['stack', 'layer', 'dodge', 'fill'],
                        selected='stack'
                    )
                ),
                selector="#main-h1_together_drop",
                where="beforeEnd",
            )
        else:
            ui.remove_ui("#inserted-dropdown_together")