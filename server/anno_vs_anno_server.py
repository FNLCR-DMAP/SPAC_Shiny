from shiny import render, reactive
from shinywidgets import render_widget
import anndata as ad
import pandas as pd
import spac.visualization
from utils.plot_manager import PlotManager


def anno_vs_anno_server(input, output, session, shared):
    sankey_pm = PlotManager('sankey', shared, plot_type='thread', data_key=None)
    relational_pm = PlotManager('relational', shared, plot_type='thread', data_key='df_relational')

    @reactive.Effect
    @reactive.event(input.go_sk1, ignore_none=True)
    def start_sankey_task():
        if sankey_pm.is_calculating.get():
            return

        adata = ad.AnnData(
            X=shared['X_data'].get(),
            obs=pd.DataFrame(shared['obs_data'].get()),
            layers=shared['layers_data'].get(),
            dtype=shared['X_data'].get().dtype
        )
        source_anno = input.sk1_anno1()
        target_anno = input.sk1_anno2()

        def worker():
            try:
                fig = spac.visualization.sankey_plot(
                    adata,
                    source_annotation=source_anno,
                    target_annotation=target_anno
                )
                sankey_pm.result_queue.put(fig)
            except Exception as e:
                sankey_pm.result_queue.put(f"Error: {str(e)}")

        sankey_pm.start_thread(worker)

    @reactive.Effect
    def check_sankey():
        def on_result(res):
            sankey_pm.result.set(res)
        sankey_pm.check_thread(on_result)

    @render_widget
    def spac_Sankey():
        return sankey_pm.result.get()

    @render.ui
    def sankey_stop_button_ui():
        return sankey_pm.stop_button_ui('stop_sankey')

    @reactive.Effect
    @reactive.event(input.go_rhm1, ignore_none=True)
    def start_relational_task():
        if relational_pm.is_calculating.get():
            return

        adata = ad.AnnData(
            X=shared['X_data'].get(),
            obs=pd.DataFrame(shared['obs_data'].get())
        )
        source_anno = input.rhm_anno1()
        target_anno = input.rhm_anno2()

        def worker():
            try:
                result = spac.visualization.relational_heatmap(
                    adata,
                    source_annotation=source_anno,
                    target_annotation=target_anno
                )
                relational_pm.result_queue.put((result['figure'], result['data']))
            except Exception as e:
                relational_pm.result_queue.put(f"Error: {str(e)}")

        relational_pm.start_thread(worker)

    @reactive.Effect
    def check_relational():
        def on_result(res):
            fig, data = res
            relational_pm.result.set(fig)
            shared['df_relational'].set(data)
        relational_pm.check_thread(on_result)

    @render_widget
    def spac_Relational():
        return relational_pm.result.get()

    @render.ui
    def relational_stop_button_ui():
        return relational_pm.stop_button_ui('stop_relational')

    @render.download(filename="relational_data.csv")
    def download_df_1():
        df = shared['df_relational'].get()
        if df is not None:
            return df.to_csv(index=False).encode("utf-8"), "text/csv"
        return None

    @render.ui
    def download_button_ui_1():
        return relational_pm.download_button_ui('download_df_1')