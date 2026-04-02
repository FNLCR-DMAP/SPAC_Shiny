from shiny import ui
from shinywidgets import output_widget


def boxplot_ui():
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
                    ui.input_action_button(
                        "go_bp",
                        "Render Plot",
                        class_="btn-success",
                        style="width: 180px;"
                    ),
                    ui.div(
                        {"style": "padding-top: 10px;"},
                        ui.output_ui("boxplot_stop_button_ui")
                    ),
                    ui.div(
                        {"style": "padding-top: 10px;"},
                        ui.output_ui("download_button_ui1")
                    )
                ),
                ui.column(
                    9,
                    ui.div(
                        {"style": "padding-bottom: 50px;"},
                        ui.panel_conditional(
                            "input.bp_output_type === false",
                            output_widget(
                                "boxplot_static",
                                width="100%",
                                height="600px"
                            )
                        ),
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