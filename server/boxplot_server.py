from shiny import render, reactive
from shinywidgets import render_widget
import anndata as ad
import pandas as pd
import spac.visualization
from utils.plot_manager import PlotManager


def boxplot_server(input, output, session, shared):
    pm = PlotManager('boxplot', shared, plot_type='thread', data_key='df_boxplot')

    def on_outlier_check():
        return None if input.bp_outlier_check() == "none" else input.bp_outlier_check()

    def on_orient_check():
        return "h" if input.bp_orient() else "v"

    def on_layer_check():
        return input.bp_layer() if input.bp_layer() != "Original" else None

    def on_anno_check():
        return input.bp_anno() if input.bp_anno() != "No Annotation" else None

    @reactive.Effect
    @reactive.event(input.go_bp, ignore_none=True)
    def start_boxplot_task():
        if pm.is_calculating.get():
            return

        adata = ad.AnnData(
            X=shared['X_data'].get(),
            obs=pd.DataFrame(shared['obs_data'].get()),
            var=pd.DataFrame(shared['var_data'].get()),
            layers=shared['layers_data'].get(),
            dtype=shared['X_data'].get().dtype
        )
        output_type = input.bp_output_type()
        annotation = on_anno_check()
        layer = on_layer_check()
        features = list(input.bp_features())
        showfliers = on_outlier_check()
        log_scale = input.bp_log_scale()
        orient = on_orient_check()

        def worker():
            try:
                if adata is None or adata.var is None:
                    pm.result_queue.put("Error: Invalid data")
                    return
                fig, df = spac.visualization.boxplot_interactive(
                    adata,
                    annotation=annotation,
                    layer=layer,
                    features=features,
                    showfliers=showfliers,
                    log_scale=log_scale,
                    orient=orient,
                    figure_height=3,
                    figure_width=4.8,
                    figure_type="interactive" if output_type else "static"
                ).values()
                pm.result_queue.put((fig, df, output_type))
            except Exception as e:
                pm.result_queue.put(f"Error: {str(e)}")

        pm.start_thread(worker)

    @reactive.Effect
    def check_status():
        def on_result(res):
            fig, df, output_type = res
            pm.result.set((fig, output_type))
            shared['df_boxplot'].set(df)
        pm.check_thread(on_result)

    @render_widget
    def spac_Boxplot():
        result = pm.result.get()
        if result is None:
            return None
        fig, output_type = result
        return fig if output_type else None

    @render_widget
    def boxplot_static():
        result = pm.result.get()
        if result is None:
            return None
        fig, output_type = result
        return fig if not output_type else None

    @render.ui
    def boxplot_stop_button_ui():
        return pm.stop_button_ui('stop_boxplot')

    @render.download(filename="boxplot_data.csv")
    def download_boxplot():
        df = shared['df_boxplot'].get()
        if df is not None:
            return df.to_csv(index=False).encode("utf-8"), "text/csv"
        return None

    @render.ui
    def download_button_ui1():
        return pm.download_button_ui('download_boxplot')