from shiny import ui, render, reactive
import anndata as ad
import numpy as np
import pandas as pd
import spac.visualization


def feat_vs_anno_server(input, output, session, shared):
    def on_layer_check():
        return input.hm1_layer() if input.hm1_layer() != "Original" else None

    def on_dendro_check():
        '''
        Check if dendrogram is enabled and return the appropriate values.
        If dendrogram is enabled,
            return a tuple (annotation dendrogram, feature dendrogram).
        If dendrogram is disabled,
            return (None, None) to indicate that no dendrogram is available.
        '''
        return (
            (input.h2_anno_dendro(), input.h2_feat_dendro())
            if input.dendogram()
            else (None, None)
        )

    @output
    @render.plot
    @reactive.event(input.go_hm1, ignore_none=True)
    def spac_Heatmap():
        adata = ad.AnnData(
            X=shared['X_data'].get(),
            obs=pd.DataFrame(shared['obs_data'].get()),
            layers=shared['layers_data'].get(),
            dtype=shared['X_data'].get().dtype
        )
        if adata:
            vmin = input.min_select()
            vmax = input.max_select()
            cmap = input.hm1_cmap()  # Get the selected color map from the dropdown
            kwargs = {"vmin": vmin,"vmax": vmax,}

            cluster_annotations, cluster_features = on_dendro_check()
            df, fig, ax = spac.visualization.hierarchical_heatmap(
                adata,
                annotation=input.hm1_anno(),
                layer=on_layer_check(),
                z_score=None,
                cluster_annotations=cluster_annotations,
                cluster_feature=cluster_features,
                **kwargs
            )

            # Only update if a non-default color map is selected
            if cmap != "viridis":
                fig.ax_heatmap.collections[0].set_cmap(cmap)

            shared['df_heatmap'].set(df)

            # Rotate x-axis labels
            fig.ax_heatmap.set_xticklabels(
                fig.ax_heatmap.get_xticklabels(),
                rotation=input.hm_x_label_rotation(),  # degrees
                horizontalalignment='right'
            )
            # fig is a seaborn.matrix.ClusterGrid
            fig.fig.subplots_adjust(bottom=0.4)
            fig.fig.subplots_adjust(left=0.1)
            return fig

        return None

    heatmap_ui_initialized = reactive.Value(False)

    @reactive.effect
    def heatmap_reactivity():
        btn = input.dendogram()
        ui_initialized = heatmap_ui_initialized.get()

        if btn and not ui_initialized:
            annotation_check = ui.input_checkbox("h2_anno_dendro", "Annotation Cluster", value=False)
            ui.insert_ui(
                ui.div({"id": "inserted-check"}, annotation_check),
                selector="#main-hm1_check",
                where="beforeEnd",
            )

            feat_check = ui.input_checkbox("h2_feat_dendro", "Feature Cluster", value=False)
            ui.insert_ui(
                ui.div({"id": "inserted-check1"}, feat_check),
                selector="#main-hm2_check",
                where="beforeEnd",
            )
            heatmap_ui_initialized.set(True)

        elif not btn and ui_initialized:
            ui.remove_ui("#inserted-check")
            ui.remove_ui("#inserted-check1")
            heatmap_ui_initialized.set(False)


    @render.download(filename="heatmap_data.csv")
    def download_df():
        df = shared['df_heatmap'].get()
        if df:
            csv_string = df.to_csv(index=False)
            csv_bytes = csv_string.encode("utf-8")
            return csv_bytes, "text/csv"
        return None


    @render.ui
    @reactive.event(input.go_hm1, ignore_none=True)
    def download_button_ui():
        if shared['df_heatmap'].get() is not None:
            return ui.download_button("download_df", "Download Data", class_="btn-warning")
        return None


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

        min_num = ui.input_numeric(
            "min_select",
            "Minimum",
            min_val,
            min=min_val,
            max=max_val
        )
        ui.insert_ui(
            ui.div({"id": "inserted-min_num"}, min_num),
            selector="#main-min_num",
            where="beforeEnd",
        )

        max_num = ui.input_numeric(
            "max_select",
            "Maximum",
            max_val,
            min=min_val,
            max=max_val
        )
        ui.insert_ui(
            ui.div({"id": "inserted-max_num"}, max_num),
            selector="#main-max_num",
            where="beforeEnd",
        )
