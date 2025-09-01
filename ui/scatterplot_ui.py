from shiny import ui

def scatterplot_ui():
    # 9. SCATTERPLOT PANEL ------------------------------------
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
                        # Added...
                        ui.input_slider(
                            "scatter_font_size",
                            "Font Size",
                            min=5,
                            max=30,
                            value=12
                        ),
                        ui.input_action_button(
                            "go_scatter",
                            "Render Plot",
                            class_="btn-success"
                        )
                    ),
                    ui.column(
                        10,
                        ui.output_plot(
                            "spac_Scatter",
                            width="100%",
                            height="80vh"
                        )
                    )
                )
            )
        )
    )