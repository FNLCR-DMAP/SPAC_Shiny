from shiny import ui
from shinywidgets import output_widget

def anno_vs_anno_ui():
    # 6. ANNO. VS ANNO. (Sankey, Relational Heatmap) ---------
   return ui.nav_panel("Anno. Vs Anno.",
        ui.card({"style": "width:100%;"},
            ui.column(
                12,
                ui.row(
                    ui.column(
                        2,
                        ui.input_select(
                           "sk1_anno1",
                           "Select Source Annotation",
                           choices=[]
                        ),
                        ui.input_select(
                           "sk1_anno2",
                           "Select Target Annotation",
                           choices=[]
                        ),
                        # Added...
                        ui.input_slider(
                            "sankey_font_size",
                            "Font Size",
                            min=5,
                            max=30,
                            value=12
                        ),
                        ui.input_action_button(
                           "go_sk1",
                           "Render Plot",
                           class_="btn-success"
                        )
                    ),
                    ui.column(
                        10,
                        ui.div(
                            output_widget("spac_Sankey"),
                            style="width:100%; height:80vh;"
                        )
                    )
                )
            )
        ),
        ui.card({"style": "width:100%;"},
            ui.column(
                12,
                ui.row(
                    ui.column(
                        2,
                        ui.input_select(
                           "rhm_anno1",
                           "Select Source Annotation",
                           choices=[],
                           selected=[]
                        ),
                        ui.input_select(
                           "rhm_anno2",
                           "Select Target Annotation",
                           choices=[],
                           selected=[]
                        ),
                        # Added...
                        ui.input_slider(
                            "heatmap_font_size",
                            "Font Size",
                            min=5,
                            max=30,
                            value=12
                        ),
                        ui.input_action_button(
                           "go_rhm1",
                           "Render Plot",
                           class_="btn-success"
                        ),
                        ui.div(
                                {"style": "padding-top: 20px;"},
                                ui.output_ui("download_button_ui_1")
                            )
                    ),
                    ui.column(
                        10,
                        ui.div(
                            output_widget("spac_Relational"),
                            style="width:100%; height:80vh;"
                        )
                    )
                )
            )
        )
    )