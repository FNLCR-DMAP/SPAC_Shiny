from shiny import ui, render, reactive
from shinywidgets import render_widget
import anndata as ad
import pandas as pd
import spac.visualization


def boxplot_server(input, output, session, shared):
   # Helper functions for reusability
    def on_outlier_check():
        selected_choice = input.bp_outlier_check()
        return None if selected_choice == "none" else selected_choice

    def on_orient_check():
        return "h" if input.bp_orient() else "v"

    def on_layer_check():
        return input.bp_layer() if input.bp_layer() != "Original" else None

    def on_anno_check():
        return input.bp_anno() if input.bp_anno() != "No Annotation" else None


    @output
    @render_widget
    @reactive.event(input.go_bp, ignore_none=True)
    def spac_Boxplot():
        """
        This function produces an interactive (Plotly) boxplot figure.
        """
        # Only run this function if both conditions are met

        if not input.bp_output_type():
            return None
        else: 

            adata = ad.AnnData(
                X=shared['X_data'].get(), 
                obs=pd.DataFrame(shared['obs_data'].get()), 
                var=pd.DataFrame(shared['var_data'].get()), 
                layers=shared['layers_data'].get(), 
                dtype=shared['X_data'].get().dtype
            )

            # Proceed only if adata is valid
            if adata is not None and adata.var is not None:

                fig, df = spac.visualization.boxplot_interactive(
                    adata, 
                    annotation=on_anno_check(), 
                    layer=on_layer_check(), 
                    features=list(input.bp_features()),
                    showfliers=on_outlier_check(),
                    log_scale=input.bp_log_scale(),
                    orient=on_orient_check(),
                    figure_height=3, 
                    figure_width=4.8, 
                    figure_type="interactive"
                ).values()

                # Return the interactive Plotly figure object
                shared['df_boxplot'].set(df)
                shared['boxplot_fig'].set(fig)  # Store figure for HTML download
                print(type(fig))
                return fig

        return None


    @render.download(filename="boxplot_data.csv")
    def download_boxplot():
        df = shared['df_boxplot'].get()
        if df is not None:
            csv_string = df.to_csv(index=False)
            csv_bytes = csv_string.encode("utf-8")
            return csv_bytes, "text/csv"
        return None

    @render.download(filename="boxplot.html")
    def download_boxplot_html():
        fig = shared['boxplot_fig'].get()
        if fig is None:
            return None
        html_string = fig.to_html(include_plotlyjs='cdn')
        html_bytes = html_string.encode("utf-8")
        return html_bytes, "text/html"


    @render.ui
    @reactive.event(input.go_bp, ignore_none=True)
    def download_button_ui1():
        if shared['df_boxplot'].get() is not None:
            return ui.input_action_button(
                "show_download_modal_bp",
                "Download Data",
                class_="btn-warning"
            )
        return None

    @reactive.Effect
    @reactive.event(input.show_download_modal_bp)
    def show_download_modal():
        m = ui.modal(
            ui.div(
                ui.download_button("download_boxplot", "CSV", class_="btn-primary me-2"),
                ui.download_button("download_boxplot_html", "HTML", class_="btn-primary"),
                style="display: flex; gap: 10px; justify-content: center;"
            ),
            title="Select a Format:",
            easy_close=True,
            footer=None
        )
        ui.modal_show(m)


    @output
    @render_widget
    @reactive.event(input.go_bp, ignore_none=True)
    def boxplot_static():
        """
        This function produces a static (Plotly) boxplot image.
        """

         # Only run this function if both conditions are met

        if input.bp_output_type():
            return None

        else: 

            adata = ad.AnnData(
                X=shared['X_data'].get(), 
                obs=pd.DataFrame(shared['obs_data'].get()), 
                var=pd.DataFrame(shared['var_data'].get()), 
                layers=shared['layers_data'].get(), 
                dtype=shared['X_data'].get().dtype
            )

            # Proceed only if adata is valid
            if adata is not None and adata.var is not None:
                
                fig, df = spac.visualization.boxplot_interactive(
                    adata, 
                    annotation=on_anno_check(), 
                    layer=on_layer_check(), 
                    features=list(input.bp_features()),
                    showfliers=on_outlier_check(),
                    log_scale=input.bp_log_scale(),
                    orient=on_orient_check(),
                    figure_height=3, 
                    figure_width=4.8, 
                    figure_type="static"
                ).values()

                return fig

        return None

