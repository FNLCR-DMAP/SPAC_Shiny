from shiny import ui

def ripleyL_ui():
    # 11. Ripley L PANEL (Plot Ripley’s L statistic) --------
    return ui.nav_panel(
        "Ripley L",
        ui.card(
            {"style": "width:100%;"},
            ui.column(
                12,
                ui.row(
                    ui.column(
                        2,
                        ui.input_select(
                            "rl_anno",
                            "Select an Annotation",
                            choices=[]
                        ),
                        ui.input_selectize(
                            "rl_label",
                            "Select two Phenotypes",
                            multiple=True,
                            choices=[],
                            selected=[]
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
                            "slide_check_rl",
                            "Stratify by Slides?",
                            False
                        ),
                        ui.panel_conditional(
                            "input.slide_check_rl === true",
                            ui.input_select(
                                "slide_select_rl",
                                "Select Slide Annotation",
                                choices=[]
                            ),
                            ui.input_select(
                                "rl_slide_labels",
                                "Select Slides",
                                choices=[],
                                selected=[]
                            ),
                        ),
                        ui.input_checkbox(
                            "sim_check_rl",
                            "Simulations",
                            False
                        ),
                        ui.panel_conditional(
                            "input.sim_check_rl === true",
                            ui.input_slider(
                                "num_sim_rl",
                                "Number of simulations",
                                min=1,
                                max=100,
                                value=1
                            ),
                            ui.input_slider(
                                "seed_rl",
                                "Seed for Simulation Reproducibility",
                                min=0,
                                max=100,
                                value=42
                            )
                        ),
                        ui.input_action_button(
                            "go_rl",
                            "Render Plot",
                            class_="btn-success"
                        ),
                        ui.div(
                            {"style": "padding-top: 20px;"},
                            ui.output_ui("download_button_ui_rl")
                        ),
                    ),
                    ui.column(
                        10,
                        ui.div(
                            {"style": "padding-bottom: 100px;"},
                            ui.output_plot(
                                "spac_ripley_l_plot",
                                width="100%",
                                height="80vh"
                            )
                        ),
                        ui.output_text_verbatim("status_msg_rl")
                    )
                )
            )
        )
    )
