from shiny import ui, render, reactive
import anndata as ad
from typing import Tuple, Any

from utils.template_wrapper import (
    register_memory_object,
    unregister_memory_object,
)
from spac.templates.visualize_ripley_template import run_from_json


def ripleyL_server(input, output, session, shared):
    """
    Server logic for Ripley L visualization using the NIDAP-derived template.

    This implementation registers the in-memory AnnData object with the
    memory registry and delegates plotting to
    `spac.templates.visualize_ripley_template.run_from_json`.

    Parameters
    ----------
    input, output, session : shiny bindings
        Standard Shiny server arguments.
    shared : dict
        Shared reactive values (expects 'adata_main' and 'df_ripley').
    """

    @reactive.calc
    def get_adata() -> ad.AnnData:
        """Retrieve the main AnnData from shared state."""
        return shared['adata_main'].get()

    @output
    @render.plot
    @reactive.event(input.go_rl, ignore_none=True)
    def spac_ripley_l_plot():
        """Render the Ripley L plot by invoking the template.

        Returns
        -------
        matplotlib.figure.Figure or None
        """
        adata = get_adata()
        if adata is None:
            return None

        # Basic inputs: read selected pair in format 'CENTER -> NEIGHBOR'
        pair = input.rl_pair() or ""
        if not pair:
            return None
        try:
            center, neighbor = [p.strip() for p in pair.split("->", 1)]
        except Exception:
            return None

        # Regions
        plot_specific_regions = bool(input.region_check_rl())
        if plot_specific_regions:
            regions_labels = input.rl_region_labels()
        else:
            regions_labels = []

        # Simulations: controlled by 'show_sim_rl' checkbox in the UI
        plot_simulations = bool(input.show_sim_rl())

        # No slide stratification: operate on the full AnnData or the
        # region-subset above. Slide-specific stratification was removed.

        # Register adata in memory registry and call run_from_json
        try:
            virtual_path = register_memory_object(adata)

            params = {
                "Upstream_Analysis": virtual_path,
                "Center_Phenotype": center,
                "Neighbor_Phenotype": neighbor,
                "Plot_Specific_Regions": plot_specific_regions,
                "Regions_Labels": regions_labels,
                "Plot_Simulations": plot_simulations,
            }

            # Call template to get figure and dataframe in-memory
            figs_df: Tuple[Any, Any] = run_from_json(
                json_path=params, save_results=False, show_plot=False
            )
            if figs_df is None:
                return None

            fig, df = figs_df
            # Store dataframe for download
            shared['df_ripley'].set(df)

            return fig

        except Exception:
            import traceback
            traceback.print_exc()
            return None

        finally:
            try:
                unregister_memory_object(virtual_path)
            except Exception:
                # ignore cleanup errors
                pass

    @render.download(filename="ripley_plot_data.csv")
    def download_df_rl():
        df = shared['df_ripley'].get()
        if df:
            csv_string = df.to_csv(index=False)
            csv_bytes = csv_string.encode("utf-8")
            return csv_bytes, "text/csv"
        return None

    @render.ui
    @reactive.event(input.go_rl, ignore_none=True)
    def download_button_ui_rl():
        if shared['df_ripley'].get():
            return ui.download_button(
                "download_df_rl",
                "Download Data",
                class_="btn-warning"
            )
        return None
