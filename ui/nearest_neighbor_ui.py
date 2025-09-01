from shiny import ui

def nearest_neighbor_ui():
    # 10. NEAREST NEIGHBORS PANEL ------------------------------------
    return ui.nav_panel(
        "Nearest Neigbors",
        ui.card(
            {"style": "width:100%;"},
            ui.column(
                12,
                ui.row(
                    ui.column(
                        2,
                        ui.input_select(
                            "nn_anno",
                            "Select an Annotation",
                            choices=[]
                        ),
                        ui.input_select(
                            "nn_anno_label",
                            "Select a Reference Phenotype",
                            choices=[]
                        ),
                        ui.input_select(
                            "nn_spatial",
                            "Select a Spatial Table",
                            choices=[]
                        ),
                        ui.input_checkbox(
                            "nn_stratify",
                            "Stratify by Annotation?",
                            False
                        ),
                        ui.panel_conditional(
                            "input.nn_stratify === true",
                            ui.input_select(
                                "nn_strat_select",
                                "Select Annotation",
                                choices=[]
                            ),
                        ),
                        ui.input_select(
                            "nn_plot_style",
                            "Select Plot Style",
                            choices=["numeric", "distribution"],
                            selected=["numeric"]
                        ),
                        ui.panel_conditional(
                            "input.nn_plot_style === 'numeric'",
                            ui.input_select(
                                "nn_plot_type_n",
                                "Select Plot Type",
                                choices=["box", "violin", "boxen"]
                            ),
                        ),
                        ui.panel_conditional(
                            "input.nn_plot_style === 'distribution'",
                            ui.input_select(
                                "nn_plot_type_d",
                                "Select Plot Type",
                                choices=["hist", "kde", "ecdf"]
                            ),
                        ),
                        ui.input_checkbox(
                            "nn_log",
                            "Apply Log to Distance Values",
                            value=False
                        ),
                        # ui.input_checkbox(
                        #     "nn_facet",
                        #     "Apply Facet Plots",
                        #     value=False
                        # ),
                        # Added...
                        ui.input_slider(
                            "nn_font_size",
                            "Font Size",
                            min=5,
                            max=30,
                            value=12
                        ),
                        ui.input_action_button(
                            "go_nn",
                            "Render Plot",
                            class_="btn_success"
                        ),
                        ui.div(
                            {"style": "padding-top: 20px;"},
                            ui.output_ui("download_button_ui_nn")
                        ),
                    ),
                    ui.column(
                        10,
                        ui.output_plot(
                            "spac_nearest_neighbor",
                            width="100%",
                            height="80vh"
                        )
                    )
                ),
            ),
        )
    )