from shiny import ui, render, reactive
from shinywidgets import render_widget
import anndata as ad
import pandas as pd
import spac.visualization


def anno_vs_anno_server(input, output, session, shared):

    @output
    @render_widget
    @reactive.event(input.go_sk1, ignore_none=True)
    def spac_Sankey():
        x_data = shared['X_data'].get()
        if x_data is None:
            return None

        cache = shared['cache']
        version = shared['dataset_version'].get()

        params = {
            'source_annotation': input.sk1_anno1(),
            'target_annotation': input.sk1_anno2(),
        }

        def compute():
            adata = ad.AnnData(
                X=x_data,
                obs=pd.DataFrame(shared['obs_data'].get()),
                layers=shared['layers_data'].get(),
                dtype=x_data.dtype
            )
            fig = spac.visualization.sankey_plot(
                adata,
                source_annotation=params['source_annotation'],
                target_annotation=params['target_annotation']
            )
            return fig, None

        fig, _ = cache.get_or_compute('sankey', version, params, compute)
        return fig

    @output
    @render_widget
    @reactive.event(input.go_rhm1, ignore_none=True)
    def spac_Relational():
        x_data = shared['X_data'].get()
        if x_data is None:
            return None

        cache = shared['cache']
        version = shared['dataset_version'].get()

        params = {
            'source_annotation': input.rhm_anno1(),
            'target_annotation': input.rhm_anno2(),
        }

        def compute():
            adata = ad.AnnData(
                X=x_data,
                obs=pd.DataFrame(shared['obs_data'].get())
            )
            result = spac.visualization.relational_heatmap(
                adata,
                source_annotation=params['source_annotation'],
                target_annotation=params['target_annotation']
            )
            return result['figure'], result['data']

        fig, df = cache.get_or_compute(
            'relational_heatmap', version, params, compute
        )

        if fig is None:
            return None

        shared['df_relational'].set(df)
        return fig

    @render.download(filename="relational_data.csv")
    def download_df_1():
        df = shared['df_relational'].get()
        if df is not None:
            csv_string = df.to_csv(index=False)
            csv_bytes = csv_string.encode("utf-8")
            return csv_bytes, "text/csv"
        return None

    @render.ui
    @reactive.event(input.go_rhm1, ignore_none=True)
    def download_button_ui_1():
        if shared['df_relational'].get() is not None:
            return ui.download_button(
                "download_df_1", "Download Data", class_="btn-warning")
        return None
