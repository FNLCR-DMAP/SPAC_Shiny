from shiny import ui
from utils.accessibility import accessible_slider


def annotations_ui():
    return ui.nav_panel(
        "Annotations",
        ui.card(
            {"style": "width:100%;"},
            ui.column(
                12,
                ui.row(
                    ui.column(
                        2,
                        ui.input_select(
                            "h2_anno",
                            "Select an Annotation",
                            choices=[]
                        ),
                        ui.input_checkbox(
                            "h2_group_by_check",
                            "Group By",
                            value=False
                        ),
                        ui.div(id="main-h2_dropdown"),
                        ui.div(id="main-h2_check"),
                        ui.div(id="main-h2_together_drop"),
                        accessible_slider(
                            "anno_slider",
                            "Rotate X-axis Labels (degrees)",
                            min_val=0,
                            max_val=90,
                            value=0,
                            step=1
                        ),
                        ui.input_action_button(
                            "go_h2",
                            "Render Plot",
                            class_="btn-success",
                            style="width: 180px;"
                        ),
                        ui.div(
                            {"style": "padding-top: 10px;"},
                            ui.output_ui("annotations_stop_button_ui")
                        ),
                        ui.div(
                            {"style": "padding-top: 10px;"},
                            ui.output_ui("download_histogram_button_ui")
                        ),
                        ui.div(
                            {"style": "padding-top: 10px;"},
                            ui.output_ui("download_histogram_plot_button_ui")
                        ),
                    ),
                    ui.column(
                        10,
                        ui.div(
                            {"style": "padding-bottom: 100px; overflow: hidden;"},
                            ui.output_image(
                                "spac_Histogram_2",
                                width="100%",
                                height="auto"
                            )
                        )
                    )
                )
            )
        )
    )