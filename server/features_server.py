from shiny import ui, render, reactive
import anndata as ad
import numpy as np
import pandas as pd
import spac.visualization
from utils.plot_utils import fig_to_png_bytes, png_bytes_to_figure


def features_server(input, output, session, shared):
    def on_layer_check():
        return input.h1_layer() if input.h1_layer() != "Original" else None


    @output
    @render.plot
    @reactive.event(input.go_h1, ignore_none=True)
    def spac_Histogram_1():
        x_data = shared['X_data'].get()
        if x_data is None:
            return None

        cache = shared['cache']
        version = shared['dataset_version'].get()

        is_grouped = input.h1_group_by_check()

        try:
            group_by = input.h1_anno() if is_grouped else None
        except Exception:
            group_by = None

        try:
            together = input.h1_together_check() if is_grouped else False
        except Exception:
            together = False

        try:
            multiple = (
                input.h1_together_drop()
                if (is_grouped and together)
                else None
            )
        except Exception:
            multiple = None

        params = {
            'feature': input.h1_feat(),
            'layer': on_layer_check(),
            'log_x': input.h1_log_x(),
            'log_y': input.h1_log_y(),
            'is_grouped': is_grouped,
            'group_by': group_by,
            'together': together,
            'multiple': multiple,
            'rotation': input.feat_slider(),
        }

        def compute():
            adata = ad.AnnData(
                X=x_data,
                obs=pd.DataFrame(shared['obs_data'].get()),
                var=pd.DataFrame(shared['var_data'].get()),
                layers=shared['layers_data'].get(),
                dtype=x_data.dtype
            )

            kwargs = {
                "adata": adata,
                "feature": params['feature'],
                "layer": params['layer'],
                "x_log_scale": params['log_x'],
                "y_log_scale": params['log_y'],
            }

            if params['is_grouped']:
                kwargs["group_by"] = params['group_by']
                kwargs["together"] = params['together']
                if params['together'] and params['multiple']:
                    kwargs["multiple"] = params['multiple']

            fig1, ax, df = spac.visualization.histogram(**kwargs).values()

            axes = ax if isinstance(ax, (list, np.ndarray)) else [ax]
            for a in axes:
                a.tick_params(
                    axis='x',
                    rotation=params['rotation'],
                    labelsize=10
                )
            return fig_to_png_bytes(fig1), df

        img_bytes, df = cache.get_or_compute('histogram1', version, params, compute)

        if img_bytes is None:
            return None

        shared['df_histogram1'].set(df)
        return png_bytes_to_figure(img_bytes)

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
