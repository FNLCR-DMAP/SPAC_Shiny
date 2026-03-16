"""
Ripley L visualization server module for SPAC Shiny application.
"""

from shiny import render, reactive
import anndata as ad
import io
import tempfile
import multiprocessing
import matplotlib.pyplot as plt
from utils.plot_manager import PlotManager
from utils.template_wrapper import (
    register_memory_object,
    unregister_memory_object,
)
from spac.templates.visualize_ripley_template import run_from_json


def run_ripley_worker(queue, adata, center, neighbor, plot_specific_regions,
                      regions_labels, plot_simulations):
    try:
        plt.clf()
        plt.close('all')

        virtual_path = register_memory_object(adata)

        params = {
            "Upstream_Analysis": virtual_path,
            "Center_Phenotype": center,
            "Neighbor_Phenotype": neighbor,
            "Plot_Specific_Regions": plot_specific_regions,
            "Regions_Labels": regions_labels,
            "Plot_Simulations": plot_simulations,
        }

        try:
            figs_df = run_from_json(
                json_path=params,
                save_results=False,
                show_plot=False
            )
        finally:
            unregister_memory_object(virtual_path)

        if figs_df is None:
            queue.put("Error: No figure generated")
            return

        fig, df = figs_df

        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=150,
                    bbox_inches='tight', pad_inches=0.2)
        buf.seek(0)
        img_bytes = buf.read()
        plt.close(fig)

        queue.put((img_bytes, df))

    except Exception as e:
        import traceback
        traceback.print_exc()
        queue.put(f"Error: {str(e)}")


def ripleyL_server(input, output, session, shared):
    pm = PlotManager('rl', shared, plot_type='process', data_key='df_ripley')

    @reactive.calc
    def get_adata():
        return shared['adata_main'].get()

    @reactive.Effect
    @reactive.event(input.go_rl, ignore_none=True)
    def start_ripley_task():
        if pm.is_calculating.get():
            return

        adata = get_adata()
        if adata is None:
            return

        pair = input.rl_pair() or ""
        if not pair:
            return

        try:
            center, neighbor = [p.strip() for p in pair.split("->", 1)]
        except Exception:
            return

        plot_specific_regions = bool(input.region_check_rl())
        regions_labels = input.rl_region_labels() if plot_specific_regions else []
        plot_simulations = bool(input.show_sim_rl())

        pm.start_process(
            run_ripley_worker,
            args=(adata, center, neighbor, plot_specific_regions,
                  regions_labels, plot_simulations)
        )

    @reactive.Effect
    def check_status():
        def on_result(res):
            img_bytes, df = res
            pm.result.set(img_bytes)
            shared['df_ripley'].set(df)
        pm.check_process(on_result)

    @output
    @render.image
    def spac_ripley_l_plot():
        img_bytes = pm.result.get()
        if img_bytes is None:
            return None
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            tmp.write(img_bytes)
            tmp_path = tmp.name
        return {
            "src": tmp_path,
            "contentType": "image/png",
            "style": "max-width: 100%; height: auto;"
        }

    @render.ui
    def ripley_stop_button_ui():
        return pm.stop_button_ui('stop_rl')

    @render.download(filename="ripley_plot_data.csv")
    def download_df_rl():
        df = shared['df_ripley'].get()
        if df is not None:
            return df.to_csv(index=False).encode("utf-8"), "text/csv"
        return None

    @render.ui
    def download_button_ui_rl():
        return pm.download_button_ui('download_df_rl')

    @render.ui
    def download_ripley_plot_button_ui():
        return pm.plot_download_button_ui('download_ripley_plot')

    @render.download(filename="ripley_plot.png")
    def download_ripley_plot():
        return pm.create_plot_download_handler()()