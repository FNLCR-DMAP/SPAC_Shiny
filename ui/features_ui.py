from shiny import ui

def features_ui():
    # 3. FEATURES PANEL (Histogram) --------------------------
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
                        ui.input_slider(
                            "feat_slider",
                            "Rotate Axis",
                            min=0,
                            max=90,
                            value=0
                        ),
                        ui.input_action_button(
                            "go_h1",
                            "Render Plot",
                            class_="btn-success"
                        ),
                        ui.div(
                            {"style": "padding-top: 20px;"},
                            ui.output_ui("download_histogram1_button_ui")
                        ),
                    ),
                    ui.column(
                        10,
                        ui.div(
                            {"style": "padding-bottom: 100px;"},
                            ui.output_plot(
                                "spac_Histogram_1",
                                width="100%",
                                height="60vh"
                            )
                        )
                    )
                )
            )
        )
    )
