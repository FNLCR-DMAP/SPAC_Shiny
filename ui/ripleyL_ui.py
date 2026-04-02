from shiny import ui


def ripleyL_ui():
    return ui.nav_panel(
        "Ripley L",
        ui.card(
            {"style": "width:100%;"},
            ui.column(
                12,
                ui.row(
                    ui.column(
                        2,
                        ui.input_selectize(
                            "rl_pair",
                            "Select Phenotype Pair (center -> neighbor)",
                            multiple=False,
                            choices=[],
                            selected=None,
                        ),
                        ui.input_checkbox(
                            "region_check_rl",
                            "Stratify by Regions?",
                            False
                        ),
                        ui.panel_conditional(
                            "input.region_check_rl === true",
                            ui.input_select(
                                "region_select_rl",
                                "Select Region Annotation",
                                choices=[]
                            ),
                            ui.input_selectize(
                                "rl_region_labels",
                                "Select Regions",
                                multiple=True,
                                choices=[],
                                selected=[]
                            ),
                        ),
                        ui.input_checkbox(
                            "show_sim_rl",
                            "Show Simulations",
                            False,
                        ),
                        ui.input_action_button(
                            "go_rl",
                            "Render Plot",
                            class_="btn-success",
                            style="width: 180px;"
                        ),
                        ui.div(
                            {"style": "padding-top: 10px;"},
                            ui.output_ui("ripley_stop_button_ui")
                        ),
                        ui.div(
                            {"style": "padding-top: 10px;"},
                            ui.output_ui("download_button_ui_rl")
                        ),
                        ui.div(
                            {"style": "padding-top: 10px;"},
                            ui.output_ui("download_ripley_plot_button_ui")
                        ),
                    ),
                    ui.column(
                        10,
                        ui.div(
                            {"style": "padding-bottom: 100px; overflow: hidden;"},
                            ui.output_image(
                                "spac_ripley_l_plot",
                                width="100%",
                                height="auto"
                            )
                        )
                    )
                )
            )
        )
    )