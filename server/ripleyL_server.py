from shiny import ui, render, reactive
import anndata as ad
import numpy as np
import pandas as pd
import spac.spatial_analysis
import spac.visualization


def ripleyL_server(input, output, session, shared):
    @output
    @render.plot
    @reactive.event(input.go_rl, ignore_none=True)
    def spac_ripley_l_plot():
        adata = ad.AnnData(
            X=shared['X_data'].get(),
            var=pd.DataFrame(shared['var_data'].get()),
            obsm=shared['obsm_data'].get(),
            obs=shared['obs_data'].get()
        )
        annotation = input.rl_anno()
        if annotation in adata.obs.columns:
            adata.obs[annotation] = adata.obs[annotation].astype(str)
        phenotypes = list(map(str, input.rl_label()))
        region_anno = None
        n_simulations = 0
        seed = None
        region_labels = None
        distances = np.linspace(0, 500, 100).tolist()

        if input.region_check_rl():
            region_anno = input.region_select_rl()
            if region_anno in adata.obs.columns:
                adata.obs[region_anno] = adata.obs[region_anno].astype(str)
            region_labels = input.rl_region_labels()
        if input.sim_check_rl():
            n_simulations = input.num_sim_rl()
            seed = input.seed_rl()
        if input.slide_check_rl():
            slide_anno = input.slide_select_rl()
            slide_label = input.rl_slide_labels()
            subset = adata.obs[adata.obs[slide_anno] == slide_label]
            check_subset = subset[annotation].unique()
            if set(phenotypes).issubset(set(check_subset)):
                adata = adata[adata.obs[slide_anno] == slide_label].copy()
            else:
                return

        # Calculate Ripley’s L statistic for spatial data in adata
        spac.spatial_analysis.ripley_l(
            adata,
            annotation=annotation,
            phenotypes=phenotypes,
            regions=region_anno,
            distances=distances,
            n_simulations=n_simulations,
            seed=seed,
        )

        # Plot Ripley L Data
        fig, df = spac.visualization.plot_ripley_l(
            adata,
            phenotypes=phenotypes,
            regions=region_labels,
            sims=input.sim_check_rl(),
            return_df=True,
        )
        shared['df_ripley'].set(df)
        return fig

    @render.download(filename="ripley_plot_data.csv")
    def download_df_rl():
        df = shared['df_ripley'].get()
        if df is not None:
            csv_string = df.to_csv(index=False)
            csv_bytes = csv_string.encode("utf-8")
            return csv_bytes, "text/csv"
        return None

    @render.ui
    @reactive.event(input.go_rl, ignore_none=True)
    def download_button_ui_rl():
        if shared['df_ripley'].get() is not None:
            return ui.download_button(
                "download_df_rl",
                "Download Data",
                class_="btn-warning"
            )
        return None
