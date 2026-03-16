from shiny import ui
def scatterplot_ui():
    return ui.nav_panel(
        "Scatterplot",
        ui.card(
            {"style": "width:100%;"},
            ui.column(
                12,
                ui.row(
                    ui.column(
                        2,
                        ui.input_select(
                            "scatter_layer",
                            "Select a Table",
                            choices=[],
                            selected="Original"
                        ),
                        ui.input_select(
                            "scatter_x",
                            "Select X Axis",
                            choices=[]
                        ),
                        ui.input_select(
                            "scatter_y",
                            "Select Y Axis",
                            choices=[]
                        ),
                        ui.input_checkbox(
                            "scatter_color_check",
                            "Color by Feature",
                            value=False
                        ),
                        ui.div(id="main-scatter_dropdown"),
                        ui.input_action_button(
                            "go_scatter",
                            "Render Plot",
                            class_="btn-success",
                            style="width: 180px;"
                        ),
                        ui.div(
                            {"style": "padding-top: 10px;"},
                            ui.output_ui("scatterplot_stop_button_ui")
                        ),
                        ui.div(
                            {"style": "padding-top: 10px;"},
                            ui.output_ui("download_scatter_plot_button_ui")
                        ),
                    ),
                    ui.column(
                        10,
                        ui.div(
                            {"style": "padding-bottom: 100px; overflow: hidden;"},
                            ui.output_image(
                                "spac_Scatter",
                                width="100%",
                                height="auto"
                            )
                        )
                    )
                )
            )
        )
    )