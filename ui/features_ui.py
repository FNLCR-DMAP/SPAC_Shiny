from shiny import ui
from utils.accessibility import accessible_slider


def features_ui():
    return ui.nav_panel(
        "Features",
        ui.card(
            {"style": "width:100%;"},
            ui.column(
                12,
                ui.row(
                    ui.column(
                        2,
                        ui.input_select(
                            "h1_feat",
                            "Select a Feature",
                            choices=[]
                        ),
                        ui.input_select(
                            "h1_layer",
                            "Select a Table",
                            choices=[],
                            selected=["Original"]
                        ),
                        ui.input_checkbox(
                            "h1_group_by_check",
                            "Group By",
                            value=False
                        ),
                        ui.input_checkbox(
                            "h1_log_x",
                            "Log X-axis",
                            value=False
                        ),
                        ui.input_checkbox(
                            "h1_log_y",
                            "Log Y-axis",
                            value=False
                        ),
                        ui.div(id="main-h1_dropdown"),
                        ui.div(id="main-h1_check"),
                        ui.div(id="main-h1_together_drop"),
                        accessible_slider(
                            "feat_slider",
                            "Rotate X-axis Labels (degrees)",
                            min_val=0,
                            max_val=90,
                            value=0,
                            step=1
                        ),
                        ui.input_action_button(
                            "go_h1",
                            "Render Plot",
                            class_="btn-success",
                            style="width: 180px;"
                        ),
                        ui.div(
                            {"style": "padding-top: 10px;"},
                            ui.output_ui("features_stop_button_ui")
                        ),
                        ui.div(
                            {"style": "padding-top: 10px;"},
                            ui.output_ui("download_histogram1_button_ui")
                        ),
                        ui.div(
                            {"style": "padding-top: 10px;"},
                            ui.output_ui("download_features_plot_button_ui")
                        ),
                    ),
                    ui.column(
                        10,
                        ui.div(
                            {"style": "padding-bottom: 100px; overflow: hidden;"},
                            ui.output_image(
                                "spac_Histogram_1",
                                width="100%",
                                height="auto"
                            )
                        )
                    )
                )
            )
        )
    )