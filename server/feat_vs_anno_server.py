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
            (input.hm1_anno_dendro(), input.hm1_feat_dendro())
            if input.hm1_dendogram()
            else (None, None)
        )

    @reactive.calc
    def get_adata():
        """Get the main AnnData object from shared state."""
        return ad.AnnData(
            X=shared['X_data'].get(),
            obs=pd.DataFrame(shared['obs_data'].get()),
            var=pd.DataFrame(shared['var_data'].get()),
            layers=shared['layers_data'].get(),
            dtype=shared['X_data'].get().dtype
        )
    
    @output
    @render.plot
    @reactive.event(input.go_hm1, ignore_none=True)
    def spac_Heatmap():
        adata = get_adata()
        if adata is None:
            return None

        vmin = input.hm1_min_select()
        vmax = input.hm1_max_select()
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

        cmap = input.hm1_cmap()
        if cmap != "viridis":
            fig.ax_heatmap.collections[0].set_cmap(cmap)

        shared['df_heatmap'].set(df)

        # Rotate X and Y axis labels
        fig.ax_heatmap.set_xticklabels(
            fig.ax_heatmap.get_xticklabels(),
            rotation=input.hm1_x_label_rotation(),
            horizontalalignment='right'
        )
        fig.ax_heatmap.set_yticklabels(
            fig.ax_heatmap.get_yticklabels(),
            rotation=input.hm1_y_label_rotation(),
            verticalalignment='center'
        )
        
        # Abbreviate labels if enabled
        def abbreviate_labels(labels, limit):
            return [label.get_text()[:limit] if label.get_text() else "" for label in labels]

        if input.hm1_enable_abbreviation():
            limit = input.hm1_label_char_limit()
            abbreviated_xticks = abbreviate_labels(fig.ax_heatmap.get_xticklabels(), limit)
            fig.ax_heatmap.set_xticklabels(abbreviated_xticks, rotation=input.hm1_x_label_rotation())
            abbreviated_yticks = abbreviate_labels(fig.ax_heatmap.get_yticklabels(), limit)
            fig.ax_heatmap.set_yticklabels(abbreviated_yticks, rotation=input.hm1_y_label_rotation())

        # Set font size for axis labels
        axis_fontsize = input.hm1_axis_label_fontsize()
        for label in fig.ax_heatmap.get_xticklabels():
            label.set_fontsize(axis_fontsize)
            label.set_fontfamily("DejaVu Sans")
        for label in fig.ax_heatmap.get_yticklabels():
            label.set_fontsize(axis_fontsize)
            label.set_fontfamily("DejaVu Sans")
        
        fig.fig.tight_layout(rect=[0.02, 0.02, 0.98, 0.98])  # Prevent the label to exceed the right border
        fig.fig.subplots_adjust(bottom=0.15, left=0)
        return fig

    @render.download(filename="heatmap_data.csv")
    def download_df_hm1():
        df = shared['df_heatmap'].get()
        if df is not None:
            csv_string = df.to_csv(index=False)
            csv_bytes = csv_string.encode("utf-8")
            return csv_bytes, "text/csv"
        return None

    @render.ui
    @reactive.event(input.go_hm1, ignore_none=True)
    def download_button_ui_hm1():
        if shared['df_heatmap'].get() is not None:
            return ui.download_button("download_df_hm1", "Download Data", class_="btn-warning")
        return None

    @reactive.effect
    @reactive.event(input.hm1_layer)
    def update_min_max():
        adata = get_adata()
        if input.hm1_layer() == "Original":
            layer_data = adata.X
        else:
            layer_data = adata.layers[input.hm1_layer()]
        mask = adata.obs[input.hm1_anno()].notna()
        layer_data = layer_data[mask]
        min_val = round(float(np.min(layer_data)), 2)
        max_val = round(float(np.max(layer_data)), 2)

        ui.remove_ui("#inserted-hm1_min_num")
        ui.remove_ui("#inserted-hm1_max_num")

        min_num = ui.input_numeric(
            "hm1_min_select", 
            "Minimum", 
            min_val, 
            min=min_val, 
            max=max_val
        )
        ui.insert_ui(
            ui.div({"id": "inserted-hm1_min_num"}, min_num),
            selector="#main-hm1_min_num",
            where="beforeEnd",
        )
        
        max_num = ui.input_numeric(
            "hm1_max_select", 
            "Maximum", 
            max_val, 
            min=min_val, 
            max=max_val
        )
        ui.insert_ui(
            ui.div({"id": "inserted-hm1_max_num"}, max_num),
            selector="#main-hm1_max_num",
            where="beforeEnd",
        )
