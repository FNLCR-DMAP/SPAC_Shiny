from shiny import ui, render, reactive
import anndata as ad
import numpy as np
import pandas as pd
import spac.visualization


def features_server(input, output, session, shared):
    @output
    @render.plot
    @reactive.event(input.go_h1, ignore_none=True)
    def spac_Histogram_1():
        adata = ad.AnnData(
            X=shared['X_data'].get(), 
            obs=pd.DataFrame(shared['obs_data'].get()), 
            var=pd.DataFrame(shared['var_data'].get()), 
            layers=shared['layers_data'].get(), 
            dtype=shared['X_data'].get().dtype
        )
        btn_log_x = input.h1_log_x()
        btn_log_y = input.h1_log_y()

        if adata is not None:
            if input.h1_group_by_check() is not True:
                if input.h1_layer() != "Original":
                    fig1, ax, df = spac.visualization.histogram(
                        adata, 
                        feature=input.h1_feat(), 
                        layer=input.h1_layer(), 
                        x_log_scale=btn_log_x, 
                        y_log_scale=btn_log_y
                    ).values()
                    ax.tick_params(
                        axis='x', 
                        rotation=input.feat_slider(), 
                        labelsize=10
                    )
                    shared['df_histogram1'].set(df)
                    return fig1
                else:
                    fig1, ax, df  = spac.visualization.histogram(
                        adata, 
                        feature=input.h1_feat(), 
                        x_log_scale=btn_log_x, 
                        y_log_scale=btn_log_y
                    ).values()
                    ax.tick_params(
                        axis='x', 
                        rotation=input.feat_slider(), 
                        labelsize=10
                    )
                    shared['df_histogram1'].set(df)
                    return fig1

            if input.h1_group_by_check() is not False:
                if input.h1_layer() != "Original":
                    if input.h1_together_check() is  not False:
                        fig1, ax, df  = spac.visualization.histogram(
                            adata, 
                            feature=input.h1_feat(), 
                            layer=input.h1_layer(), 
                            group_by=input.h1_anno(), 
                            together=input.h1_together_check(), x_log_scale=btn_log_x, 
                            y_log_scale=btn_log_y, 
                            multiple=input.h1_together_drop()
                        ).values()
                        ax.tick_params(
                            axis='x', 
                            rotation=input.feat_slider(), 
                            labelsize=10
                        )
                        shared['df_histogram1'].set(df)
                        return fig1
                    else:
                        fig1, ax, df  = spac.visualization.histogram(
                            adata, 
                            feature=input.h1_feat(), 
                            layer=input.h1_layer(), 
                            group_by=input.h1_anno(), 
                            together=input.h1_together_check(), x_log_scale=btn_log_x, 
                            y_log_scale=btn_log_y
                        ).values()
                        axes = (
                            ax if isinstance(ax, (list, np.ndarray))
                            else [ax]
                        )
                        for ax in axes:
                            ax.tick_params(
                                axis='x', 
                                rotation=input.feat_slider(), 
                                labelsize=10
                            )
                            shared['df_histogram1'].set(df)
                        return fig1
                else:
                    if input.h1_together_check() is  not False:
                        fig1, ax, df  = spac.visualization.histogram(
                            adata, 
                            feature=input.h1_feat(), 
                            group_by=input.h1_anno(), 
                            together=input.h1_together_check(), x_log_scale=btn_log_x, 
                            y_log_scale=btn_log_y, 
                            multiple=input.h1_together_drop()
                        ).values()
                        ax.tick_params(
                            axis='x', 
                            rotation=input.feat_slider(), 
                            labelsize=10
                        )
                        shared['df_histogram1'].set(df)
                        return fig1
                    else:
                        fig1, ax, df  = spac.visualization.histogram(
                            adata, 
                            feature=input.h1_feat(), 
                            group_by=input.h1_anno(), 
                            together=input.h1_together_check(), x_log_scale=btn_log_x, 
                            y_log_scale=btn_log_y
                        ).values()
                        axes = (
                            ax if isinstance(ax, (list, np.ndarray)) 
                            else [ax]
                        )
                        for ax in axes:
                            ax.tick_params(
                                axis='x', 
                                rotation=input.feat_slider(), 
                                labelsize=10
                            )
                            shared['df_histogram1'].set(df)
                        return fig1
        return None

    histogram_ui_initialized = reactive.Value(False)


    @render.download(filename="features_histogram_data.csv")
    def download_histogram1_df():
        df = shared['df_histogram1'].get()
        if df is not None:
            csv_string = df.to_csv(index=False)
            csv_bytes = csv_string.encode("utf-8")
            return csv_bytes, "text/csv"
        return None


    @render.ui
    @reactive.event(input.go_h1, ignore_none=True)
    def download_histogram1_button_ui():
        if shared['df_histogram1'].get() is not None:
            return ui.download_button(
                "download_histogram1_df", 
                "Download Data", 
                class_="btn-warning"
            )
        return None


    @reactive.effect
    def histogram_reactivity():
        btn = input.h1_group_by_check()
        ui_initialized = histogram_ui_initialized.get()

        if btn and not ui_initialized:
            dropdown = ui.input_select(
                "h1_anno", 
                "Select an Annotation", 
                choices=shared['obs_names'].get()
            )
            ui.insert_ui(
                ui.div({"id": "inserted-dropdown"}, dropdown),
                selector="#main-h1_dropdown",
                where="beforeEnd",
            )

            together_check = ui.input_checkbox(
                "h1_together_check", 
                "Plot Together", 
                value=True
            )
            ui.insert_ui(
                ui.div({"id": "inserted-check"}, together_check),
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
            dropdown_together = ui.input_select(
                "h1_together_drop", 
                "Select Stack Type", 
                choices=['stack', 'layer', 'dodge', 'fill'], 
                selected='stack'
            )
            ui.insert_ui(
                ui.div(
                    {"id": "inserted-dropdown_together"}, 
                    dropdown_together
                ),
                selector="#main-h1_together_drop",
                where="beforeEnd",)      
        else:
            ui.remove_ui("#inserted-dropdown_together")
