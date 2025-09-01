from shiny import ui, render, reactive
from shinywidgets import render_widget
import anndata as ad
import pandas as pd
import spac.visualization
# Added...
import matplotlib.pyplot as plt


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
            # Added...
            font_size = input.bp_font_size()

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
                # Added...
                fig.update_layout(font=dict(size=font_size))
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


    @render.ui
    @reactive.event(input.go_bp, ignore_none=True)
    def download_button_ui1():
        if shared['df_boxplot'].get() is not None:
            return ui.download_button(
                "download_boxplot",
                "Download Data",
                class_="btn-warning"
            )
        return None


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
            # Added...
            font_size = input.bp_font_size()

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
                # Added...
                fig.update_layout(font=dict(size=font_size))
                return fig

        return None