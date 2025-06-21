from shiny import ui, render, reactive
import spac.spatial_analysis
import spac.visualization


def nearest_neighbor_server(input, output, session, shared):
    @output
    @render.plot
    @reactive.event(input.go_nn, ignore_none=True)
    def spac_nearest_neighbor():
        adata = shared['adata_main'].get()
        annotation = input.nn_anno()
        label = str(input.nn_anno_label())
        if input.nn_plot_style() == 'numeric':
            plot_type = input.nn_plot_type_n()
        else:
            plot_type = input.nn_plot_type_d()
        if input.nn_stratify():
            stratify_by = input.nn_strat_select()
        else:
            stratify_by = None
        if annotation in adata.obs.columns:
            adata.obs[annotation] = adata.obs[annotation].astype(str)

        spac.spatial_analysis.calculate_nearest_neighbor(
            adata,
            annotation,
            spatial_associated_table=input.nn_spatial(),
            imageid=stratify_by,
            label='spatial_distance',
            verbose=True
        )
        if adata is not None:
            out = spac.visualization.visualize_nearest_neighbor(
                adata=adata,
                annotation=annotation,
                distance_from=label,
                method=input.nn_plot_style(),
                log=input.nn_log(),
                facet_plot=True,
                plot_type=plot_type,
                stratify_by=stratify_by
            )
            shared['df_nn'].set(out['data'])
            return out['fig']

    @render.download(filename="nearest_neighbor_data.csv")
    def download_df_nn():
        df = shared['df_nn'].get()
        if df is not None:
            csv_string = df.to_csv(index=False)
            csv_bytes = csv_string.encode("utf-8")
            return csv_bytes, "text/csv"
        return None

    @render.ui
    @reactive.event(input.go_nn, ignore_none=True)
    def download_button_ui_nn():
        if shared['df_nn'].get() is not None:
            return ui.download_button(
                "download_df_nn",
                "Download Data",
                class_="btn-warning"
            )
        return None
