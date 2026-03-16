from shiny import ui
from utils.accessibility import accessible_slider


def feat_vs_anno_ui():
    return ui.nav_panel("Feat. Vs Anno.",
        ui.card({"style": "width:100%;"},
            ui.column(
                12,
                ui.row(
                    ui.column(
                        2,
                        ui.input_select(
                            "hm1_anno",
                            "Select an Annotation",
                            choices=[]
                        ),
                        ui.input_select(
                            "hm1_layer",
                            "Select a Table",
                            choices=[]
                        ),
                        ui.input_select(
                            "hm1_cmap",
                            "Select Color Map",
                            choices=[
                                "viridis", "plasma", "inferno", "magma",
                                "cividis", "coolwarm", "RdYlBu", "Spectral",
                                "PiYG", "PRGn"
                            ]
                        ),
                        accessible_slider(
                            "hm_x_label_rotation",
                            "Rotate X Axis Labels (degrees)",
                            min_val=0,
                            max_val=90,
                            value=50,
                            step=1
                        ),
                        ui.input_checkbox(
                            "dendogram",
                            "Include Dendrogram",
                            False
                        ),
                        ui.div(id="main-hm1_check"),
                        ui.div(id="main-hm2_check"),
                        ui.div(id="main-min_num"),
                        ui.div(id="main-max_num"),
                        ui.input_action_button(
                            "go_hm1",
                            "Render Plot",
                            class_="btn-success",
                            style="width: 180px;"
                        ),
                        ui.div(
                            {"style": "padding-top: 10px;"},
                            ui.output_ui("heatmap_stop_button_ui")
                        ),
                        ui.div(
                            {"style": "padding-top: 10px;"},
                            ui.output_ui("download_button_ui")
                        ),
                        ui.div(
                            {"style": "padding-top: 10px;"},
                            ui.output_ui("download_heatmap_plot_button_ui")
                        ),
                    ),
                    ui.column(
                        10,
                        ui.div(
                            {"style": "padding-bottom: 100px; overflow: hidden;"},
                            ui.output_image(
                                "spac_Heatmap",
                                width="100%",
                                height="auto"
                            )
                        )
                    )
                )
            )
        )
    )