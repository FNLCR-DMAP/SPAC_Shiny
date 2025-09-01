from shiny import ui, render, reactive
import numpy as np
import spac.visualization
# Added...
import matplotlib.pyplot as plt

def annotations_server(input, output, session, shared):
    @output
    @render.plot
    @reactive.event(input.go_h2, ignore_none=True)
    def spac_Histogram_2():
        adata = shared['adata_main'].get()
        if adata is None:
            return None

        # Added...
        # Note: This assumes your UI file has a slider with the id 'annotations_font_size'.
        # Please ensure this ID matches the one in your annotations_ui.py file.
        font_size = input.annotations_font_size()
        plt.rcParams.update({'font.size': font_size})
        rotation = input.anno_slider()

        # 1) If "Group By" is UNCHECKED, show a simple annotation histogram
        if not input.h2_group_by_check():
            fig, ax, df = spac.visualization.histogram(
                adata,
                annotation=input.h2_anno()
            ).values()
            shared['df_histogram2'].set(df)
            # Modified...
            ax.tick_params(axis='x', rotation=rotation, labelsize=font_size)
            return fig

        # 2) If "Group By" is CHECKED, we must always supply a
        #    valid multiple parameter
        else:
            # If user also checked "Plot Together", use their selected
            # stack type
            if input.h2_together_check():
                # e.g. 'stack', 'dodge', etc.
                multiple_param = input.h2_together_drop()

                together_flag = True
            else:
                # If grouping by but not "plot together", pick a default layout
                # or 'dodge' or any valid string
                multiple_param = "layer"
                together_flag = False

            fig, ax, df = spac.visualization.histogram(
                adata,
                annotation=input.h2_anno(),
                group_by=input.h2_anno_1(),
                together=together_flag,
                multiple=multiple_param
            ).values()
            shared['df_histogram2'].set(df)
            axes = ax if isinstance(ax, (list, np.ndarray)) else [ax]
            # Modified... (renamed loop variable to avoid shadowing)
            for current_ax in axes:
                # Modified...
                current_ax.tick_params(
                    axis='x',
                    rotation=rotation,
                    labelsize=font_size
                )
            return fig
        return None


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
        df = shared['df_human_histogram2'].get()
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