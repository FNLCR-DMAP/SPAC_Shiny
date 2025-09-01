from shiny import ui
from shinywidgets import output_widget


def boxplot_ui():
    # 4. BOXPLOTS PANEL --------------------------------------
    return ui.nav_panel(
        "Boxplot",
        ui.card(
            {"style": "width:100%;"},
            ui.row(
                ui.column(
                    3,
                    ui.input_select(
                        "bp_anno",
                        "Select an Annotation",
                        choices=[]
                    ),
                    ui.input_selectize(
                        "bp_features",
                        "Select Features",
                        multiple=True,
                        choices=[],
                        selected=[]
                    ),
                    ui.input_select(
                        "bp_layer",
                        "Select a Table",
                        choices=[],
                        selected="Original"
                    ),
                    ui.input_select(
                        "bp_outlier_check",
                        "Add Outliers",
                        choices={
                            "all": "All",
                            "downsample": "Downsampled",
                            "none": "None"
                        },
                        selected="none"
                    ),
                    ui.input_checkbox(
                        "bp_log_scale",
                        "Log Scale",
                        False
                    ),
                    ui.input_checkbox(
                        "bp_orient",
                        "Horizontal Orientation",
                        False
                    ),
                    ui.input_checkbox(
                        "bp_output_type",
                        "Enable Interactive Plot",
                        True
                    ),
                    # Added...
                    ui.input_slider(
                        "bp_font_size",
                        "Font Size",
                        min=5,
                        max=30,
                        value=12
                    ),
                    ui.input_action_button(
                        "go_bp",
                        "Render Plot",
                        class_="btn-success"
                    ),
                    ui.div(
                        {"style": "padding-top: 20px;"},
                        ui.output_ui("download_button_ui1")
                    )
                ),
                ui.column(
                    9,
                    ui.div(
                        {"style": "padding-bottom: 50px;"},
                        # Static plot conditional panel
                        # (when interactive unchecked)
                        ui.panel_conditional(
                            "input.bp_output_type === false",
                            output_widget(
                                "boxplot_static",
                                width="100%",
                                height="600px"
                            )
                        ),
                        # Interactive plot conditional panel
                        # (when interactive checked)
                        ui.panel_conditional(
                            "input.bp_output_type === true",
                            output_widget(
                                "spac_Boxplot",
                                width="100%",
                                height="600px"
                            )
                        )
                    )
                ),
            )
        ),
    )