from shiny import ui, render, reactive
import anndata as ad
import numpy as np
import pandas as pd
import spac.visualization


def feat_vs_anno_server(input, output, session, shared):
    rendering_state = reactive.Value(False)

    @reactive.effect
    @reactive.event(input.go_hm1)
    def handle_render_start():
        rendering_state.set(True)

    @reactive.effect
    @reactive.event(input.cancel_hm1)
    def handle_cancel_click():
        rendering_state.set(False)

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
            var=pd.DataFrame(shared['var_data'].get()),
            layers=shared['layers_data'].get(),
            dtype=shared['X_data'].get().dtype
        )
        if adata is None:
            return None

        vmin = input.min_select()
        vmax = input.max_select()
        cmap = input.hm1_cmap()
        fontsize = input.axis_label_fontsize()
        kwargs = {"vmin": vmin, "vmax": vmax}
        cluster_annotations, cluster_features = on_dendro_check()

        try:
            df, fig, ax = spac.visualization.hierarchical_heatmap(
                adata,
                annotation=input.hm1_anno(),
                layer=on_layer_check(),
                z_score=None,
                cluster_annotations=cluster_annotations,
                cluster_feature=cluster_features,
                **kwargs
            )
        except Exception as e:
            print("Heatmap generation failed:", e)
            return None

        if fig is None or not hasattr(fig, "ax_heatmap"):
            print("Invalid figure structure.")
            return None

        if cmap != "viridis":
            fig.ax_heatmap.collections[0].set_cmap(cmap)

        shared['df_heatmap'].set(df)

        #Rotate X and Y axis labels
        fig.ax_heatmap.set_xticklabels(
            fig.ax_heatmap.get_xticklabels(),
            rotation=input.hm_x_label_rotation(),
            horizontalalignment='right'
        )
        fig.ax_heatmap.set_yticklabels(
            fig.ax_heatmap.get_yticklabels(),
            rotation=input.hm_y_label_rotation(),
            verticalalignment='center'
        )
        
        # Abbreviate labels if enabled
        def abbreviate_labels(labels, limit):
            return [label.get_text()[:limit] if label.get_text() else "" for label in labels]

        if input.enable_abbreviation():
            limit = input.label_char_limit()
            abbreviated_xticks = abbreviate_labels(fig.ax_heatmap.get_xticklabels(), limit)
            fig.ax_heatmap.set_xticklabels(abbreviated_xticks, rotation=input.hm_x_label_rotation())
            abbreviated_yticks = abbreviate_labels(fig.ax_heatmap.get_yticklabels(), limit)
            fig.ax_heatmap.set_yticklabels(abbreviated_yticks, rotation=input.hm_y_label_rotation())

        for label in fig.ax_heatmap.get_xticklabels():
            label.set_fontsize(fontsize)
            label.set_fontfamily("DejaVu Sans")
        for label in fig.ax_heatmap.get_yticklabels():
            label.set_fontsize(fontsize)
            label.set_fontfamily("DejaVu Sans")

        fig.fig.subplots_adjust(bottom=0.4, left=0.1)
        return fig

    heatmap_ui_initialized = reactive.Value(False)

    @reactive.effect
    def heatmap_reactivity():
        btn = input.dendogram()
        ui_initialized = heatmap_ui_initialized.get()

        if btn and not ui_initialized:
            # Insert feature cluster first
            feat_check = ui.input_checkbox("h2_feat_dendro", "Feature Cluster", value=False)
            ui.insert_ui(
                ui.div({"id": "inserted-check1"}, feat_check),
                selector="#main-hm2_check",
                where="beforeEnd",
            )
            # Insert annotation cluster below
            annotation_check = ui.input_checkbox("h2_anno_dendro", "Annotation Cluster", value=False)
            ui.insert_ui(
                ui.div({"id": "inserted-check"}, annotation_check),
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
        if df is not None:
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

    @reactive.effect
    @reactive.event(input.enable_abbreviation)
    def toggle_label_char_limit_slider():
        if input.enable_abbreviation():
            slider = ui.input_slider(
                "label_char_limit",
                "Max Characters per Label",
                min=2,
                max=20,
                value=6
            )
            ui.insert_ui(
                ui.div({"id": "inserted-label-char-limit"}, slider),
                selector="#main-hm1_check",  # Or another appropriate container
                where="beforeEnd",
            )
        else:
            ui.remove_ui("#inserted-label-char-limit")
