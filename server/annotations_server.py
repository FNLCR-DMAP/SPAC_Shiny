from shiny import ui, render, reactive
import numpy as np
import io
import tempfile
import matplotlib.pyplot as plt
import multiprocessing
import spac.visualization
from utils.plot_manager import PlotManager


def run_histogram_worker(queue, adata, annotation, group_by, together, multiple, rotation):
    try:
        plt.clf()
        plt.close('all')

        res_dict = spac.visualization.histogram(
            adata,
            annotation=annotation,
            group_by=group_by,
            together=together,
            multiple=multiple
        )
        fig, ax, df = res_dict.values()

        fig.set_size_inches(12, 7)
        axes = ax if isinstance(ax, (list, np.ndarray)) else [ax]
        for a in axes:
            a.tick_params(axis='x', rotation=rotation, labelsize=10)

        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=140, bbox_inches=None, pad_inches=0.2)
        buf.seek(0)
        img_bytes = buf.read()
        plt.close(fig)

        queue.put((img_bytes, df))
    except Exception as e:
        queue.put(f"Error: {str(e)}")


def annotations_server(input, output, session, shared):
    pm = PlotManager('histogram', shared, plot_type='process', data_key='df_histogram2')

    @reactive.Effect
    @reactive.event(input.go_h2, ignore_none=True)
    def start_histogram_task():
        if pm.is_calculating.get():
            return

        adata = shared['adata_main'].get()
        if adata is None:
            return

        is_grouped = input.h2_group_by_check()
        annotation = input.h2_anno()
        group_by = input.h2_anno_1() if is_grouped else None
        together = input.h2_together_check() if is_grouped else False
        multiple = input.h2_together_drop() if (is_grouped and input.h2_together_check()) else "layer"
        rotation = input.anno_slider()

        pm.start_process(
            run_histogram_worker,
            args=(adata, annotation, group_by, together, multiple, rotation)
        )

    @reactive.Effect
    def check_status():
        def on_result(res):
            img_bytes, df = res
            pm.result.set(img_bytes)
            shared['df_histogram2'].set(df)
        pm.check_process(on_result)

    @output
    @render.image
    def spac_Histogram_2():
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
    def annotations_stop_button_ui():
        return pm.stop_button_ui('stop_h2')

    @render.ui
    def download_histogram_button_ui():
        return pm.download_button_ui('download_histogram2_df')

    @render.download(filename="annotation_histogram_data.csv")
    def download_histogram2_df():
        df = shared['df_histogram2'].get()
        if df is not None:
            return df.to_csv(index=False).encode("utf-8"), "text/csv"
        return None

    histogram2_ui_initialized = reactive.Value(False)
    
    @render.ui
    def download_histogram_plot_button_ui():
        return pm.plot_download_button_ui('download_histogram2_plot')

    @render.download(filename="annotation_plot.png")
    def download_histogram2_plot():
        return pm.create_plot_download_handler()()

    @reactive.effect
    def histogram_reactivity_2():
        btn = input.h2_group_by_check()
        ui_initialized = histogram2_ui_initialized.get()

        if btn and not ui_initialized:
            ui.insert_ui(
                ui.div(
                    {"id": "inserted-dropdown-1"},
                    ui.input_select(
                        "h2_anno_1",
                        "Select an Annotation",
                        choices=shared['obs_names'].get()
                    )
                ),
                selector="#main-h2_dropdown",
                where="beforeEnd",
            )
            ui.insert_ui(
                ui.div(
                    {"id": "inserted-check-1"},
                    ui.input_checkbox("h2_together_check", "Plot Together", value=True)
                ),
                selector="#main-h2_check",
                where="beforeEnd",
            )
            histogram2_ui_initialized.set(True)

        elif not btn and ui_initialized:
            ui.remove_ui("#inserted-dropdown-1")
            ui.remove_ui("#inserted-check-1")
            ui.remove_ui("#inserted-dropdown_together-1")
            histogram2_ui_initialized.set(False)

    @reactive.effect
    @reactive.event(input.h2_together_check)
    def update_stack_type_dropdown():
        if input.h2_together_check():
            ui.insert_ui(
                ui.div(
                    {"id": "inserted-dropdown_together-1"},
                    ui.input_select(
                        "h2_together_drop",
                        "Select Stack Type",
                        choices=['stack', 'layer', 'dodge', 'fill'],
                        selected='stack'
                    )
                ),
                selector="#main-h2_together_drop",
                where="beforeEnd"
            )
        else:
            ui.remove_ui("#inserted-dropdown_together-1")