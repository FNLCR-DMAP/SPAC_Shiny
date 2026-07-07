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
        if adata is None:
            return None
        fig = spac.visualization.sankey_plot(
            adata,
            source_annotation=input.sk1_anno1(),
            target_annotation=input.sk1_anno2()
        )
        shared['sankey_fig'].set(fig)  # Store figure for HTML download
        return fig

    @output
    @render_widget
    @reactive.event(input.go_rhm1, ignore_none=True)
    def spac_Relational():
        adata = ad.AnnData(
            X=shared['X_data'].get(),
            obs=pd.DataFrame(shared['obs_data'].get())
        )
        if adata is None:
            return None
        result = spac.visualization.relational_heatmap(
            adata,
            source_annotation=input.rhm_anno1(),
            target_annotation=input.rhm_anno2()
        )
        shared['df_relational'].set(result['data'])
        return result['figure']


    @render.download(filename="relational_data.csv")
    def download_df_1():
        df = shared['df_relational'].get()
        if df is None:
            return None
        csv_string = df.to_csv(index=False)
        csv_bytes = csv_string.encode("utf-8")
        return csv_bytes, "text/csv"


    @render.ui
    @reactive.event(input.go_rhm1, ignore_none=True)
    def download_button_ui_1():
        if shared['df_relational'].get() is None:
            return None
        return ui.download_button("download_df_1", "Download Data", class_="btn-warning")
    
    def get_sankey_html_filename():
        """Generate HTML download filename for sankey."""
        input_filename = shared['input_filename'].get()
        if input_filename:
            return f"{input_filename}_sankey.html"
        return "sankey.html"
    
    @render.download(filename=get_sankey_html_filename)
    def download_sankey_html():
        fig = shared['sankey_fig'].get()
        if fig is None:
            return None
        html_string = fig.to_html(include_plotlyjs='cdn')
        html_bytes = html_string.encode("utf-8")
        return html_bytes, "text/html"
    
    @render.ui
    @reactive.event(input.go_sk1, ignore_none=True)
    def download_button_ui_sankey():
        if shared['sankey_fig'].get() is None:
            return None
        return ui.input_action_button(
            "show_download_modal_sankey",
            "Download Data",
            class_="btn-warning"
        )
    
    @reactive.Effect
    @reactive.event(input.show_download_modal_sankey)
    def show_download_modal_sankey():
        m = ui.modal(
            ui.div(
                ui.download_button("download_sankey_html", "HTML", class_="btn-primary"),
                style="display: flex; gap: 10px; justify-content: center;"
            ),
            title="Select a Format:",
            easy_close=True,
            footer=None
        )
        ui.modal_show(m)
