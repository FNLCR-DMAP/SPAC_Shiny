from shiny import ui, render, reactive
import numpy as np
import spac.visualization
from utils.plot_utils import fig_to_png_bytes, png_bytes_to_figure


def annotations_server(input, output, session, shared):
    @output
    @render.plot
    @reactive.event(input.go_h2, ignore_none=True)
    def spac_Histogram_2():
        adata = shared['adata_main'].get()
        if adata is None:
            return None

        cache = shared['cache']
        version = shared['dataset_version'].get()

        is_grouped = input.h2_group_by_check()

        # Safely read dynamic inputs that may not exist yet
        try:
            group_by = input.h2_anno_1() if is_grouped else None
        except Exception:
            group_by = None

        try:
            together = input.h2_together_check() if is_grouped else False
        except Exception:
            together = False

        try:
            multiple = (
                input.h2_together_drop()
                if (is_grouped and together)
                else "layer"
            )
        except Exception:
            multiple = "layer"

        rotation = input.anno_slider()

        params = {
            'annotation': input.h2_anno(),
            'is_grouped': is_grouped,
            'group_by': group_by,
            'together': together,
            'multiple': multiple,
            'rotation': rotation,
        }

        def compute():
            kwargs = {'adata': adata, 'annotation': params['annotation']}
            if params['is_grouped']:
                kwargs['group_by'] = params['group_by']
                kwargs['together'] = params['together']
                if params['together']:
                    kwargs['multiple'] = params['multiple']
                else:
                    kwargs['multiple'] = "layer"

            fig, ax, df = spac.visualization.histogram(**kwargs).values()

            axes = ax if isinstance(ax, (list, np.ndarray)) else [ax]
            for a in axes:
                a.tick_params(
                    axis='x',
                    rotation=params['rotation'],
                    labelsize=10
                )
            return fig_to_png_bytes(fig), df

        img_bytes, df = cache.get_or_compute('histogram2', version, params, compute)

        if img_bytes is None:
            return None

        shared['df_histogram2'].set(df)
        return png_bytes_to_figure(img_bytes)


    @render.ui
    @reactive.event(input.go_h2, ignore_none=True)
    def download_histogram_button_ui():
        if shared['df_histogram2'].get() is not None:
            return ui.download_button(
                "download_histogram2_df",
                "Download Data",
                class_="btn-warning"
            )
        return None


    @render.download(filename="annotation_histogram_data.csv")
    def download_histogram2_df():
        df = shared['df_histogram2'].get()
        if df is not None:
            csv_string = df.to_csv(index=False)
            csv_bytes = csv_string.encode("utf-8")
            return csv_bytes, "text/csv"
        return None

    histogram2_ui_initialized = reactive.Value(False)

    @reactive.effect
    def histogram_reactivity_2():
        btn = input.h2_group_by_check()
        ui_initialized = histogram2_ui_initialized.get()

        if btn and not ui_initialized:
            dropdown = ui.input_select(
                "h2_anno_1",
                "Select an Annotation",
                choices=shared['obs_names'].get()
            )
            ui.insert_ui(
                ui.div({"id": "inserted-dropdown-1"}, dropdown),
                selector="#main-h2_dropdown",
                where="beforeEnd",
            )

            together_check = ui.input_checkbox(
                "h2_together_check",
                "Plot Together",
                value=True
            )
            ui.insert_ui(
                ui.div({"id": "inserted-check-1"}, together_check),
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
            dropdown_together = ui.input_select(
                "h2_together_drop",
                "Select Stack Type",
                choices=['stack', 'layer', 'dodge', 'fill'],
                selected='stack'
            )
            ui.insert_ui(
                ui.div({
                    "id": "inserted-dropdown_together-1"},
                    dropdown_together
                ),
                selector="#main-h2_together_drop",
                where="beforeEnd"
            )
        else:
            ui.remove_ui("#inserted-dropdown_together-1")
