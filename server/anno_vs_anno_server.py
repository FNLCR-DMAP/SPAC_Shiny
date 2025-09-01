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
        adata = ad.AnnData(
            X=shared['X_data'].get(),
            obs=pd.DataFrame(shared['obs_data'].get()),
            layers=shared['layers_data'].get(),
            dtype=shared['X_data'].get().dtype
        )
        if adata is not None:
            fig = spac.visualization.sankey_plot(
                adata,
                source_annotation=input.sk1_anno1(),
                target_annotation=input.sk1_anno2()
            )
            
            font_size = input.sankey_font_size()
            
            # Modified... 
            # Applying font size directly to the Sankey trace for node and 
            # label text, as the global layout font can sometimes be ignored.
            fig.update_layout(font=dict(size=font_size))
            fig.update_traces(textfont=dict(size=font_size), selector=dict(type='sankey'))

            return fig
        return None

    @output
    @render_widget
    @reactive.event(input.go_rhm1, ignore_none=True)
    def spac_Relational():
        adata = ad.AnnData(
            X=shared['X_data'].get(),
            obs=pd.DataFrame(shared['obs_data'].get())
        )
        if adata is not None:
            result = spac.visualization.relational_heatmap(
                adata,
                source_annotation=input.rhm_anno1(),
                target_annotation=input.rhm_anno2()
            )
            shared['df_relational'].set(result['data'])

            font_size = input.heatmap_font_size()

            # Modified...
            # Applying font size directly to the heatmap axes and color bar,
            # as these elements often have their own font settings.
            result['figure'].update_layout(
                font=dict(size=font_size),
                xaxis=dict(tickfont=dict(size=font_size)),
                yaxis=dict(tickfont=dict(size=font_size))
            )
            result['figure'].update_coloraxes(colorbar_tickfont_size=font_size)

            return result['figure']
        return None

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
            return ui.download_button("download_df_1", "Download Data", class_="btn-warning")
        return None