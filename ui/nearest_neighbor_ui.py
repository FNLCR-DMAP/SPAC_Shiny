"""
Nearest neighbor visualization UI module for SPAC Shiny application.
"""

from shiny import ui


def nearest_neighbor_ui():
    return ui.nav_panel(
        "Nearest Neighbors",
        ui.tags.style("""
            .nn-controls-panel {
                background-color: #f8f9fa;
                border-radius: 8px;
                padding: 10px;
            }
            .nn-plot-container {
                background-color: white;
                border-radius: 8px;
                border: 1px solid #dee2e6;
                padding: 10px;
            }
            .nn-controls-panel .form-group {
                margin-bottom: 12px;
            }
            .nn-controls-panel h4, .nn-controls-panel h5 {
                margin-bottom: 10px;
                margin-top: 15px;
            }
            .nn-controls-panel h4:first-child {
                margin-top: 0;
            }
            .accessible-tooltip:focus {
                outline: 2px solid #0056b3;
                background: #e9ecef;
            }
        """),
        ui.card(
            {"style": "width:100%;"},
            ui.column(
                12,
                ui.row(
                    ui.column(
                        3,
                        ui.div(
                            {
                                "class": "nn-controls-panel",
                                "style": (
                                    "height: 85vh; overflow-y: auto; "
                                    "padding-right: 10px; "
                                    "border-right: 1px solid #dee2e6;"
                                )
                            },
                            ui.div(
                                ui.h4("Core Parameters",
                                      class_="accessible-heading"),
                                ui.input_select(
                                    "nn_source_label",
                                    "Source Anchor Cell Label",
                                    choices=[]
                                ),
                                ui.input_selectize(
                                    "nn_target_label",
                                    "Target Cell Label",
                                    choices=[],
                                    multiple=True,
                                    options={
                                        "placeholder":
                                        "Select targets (empty for 'All')"
                                    }
                                ),
                                ui.input_select(
                                    "nn_image_id",
                                    ui.tags.span(
                                        "ImageID",
                                        ui.tags.span(
                                            "\u24D8",
                                            title=(
                                                "The annotation name used to distinguish "
                                                "different images or samples. "
                                                "If there's only one image, set to 'None.'"
                                            ),
                                            tabindex="0",
                                            class_="accessible-tooltip",
                                            style="margin-left:5px; cursor:help; color:#007bff;"
                                        )
                                    ),
                                    choices=["None"]
                                ),
                                ui.hr(),
                                ui.div(
                                    ui.input_checkbox(
                                        "nn_show_plot_config",
                                        "Show Plot Configuration",
                                        value=False
                                    ),
                                    ui.panel_conditional(
                                        "input.nn_show_plot_config",
                                        ui.input_select(
                                            "nn_plot_method",
                                            "Plot Method",
                                            choices=["numeric", "distribution"],
                                            selected="numeric"
                                        ),
                                        ui.panel_conditional(
                                            "input.nn_plot_method === 'numeric'",
                                            ui.input_select(
                                                "nn_plot_type_numeric",
                                                "Plot Type",
                                                choices=["boxen", "box", "violin",
                                                         "strip", "swarm"],
                                                selected="boxen"
                                            ),
                                        ),
                                        ui.panel_conditional(
                                            "input.nn_plot_method === 'distribution'",
                                            ui.input_select(
                                                "nn_plot_type_distribution",
                                                "Plot Type",
                                                choices=["kde", "hist", "ecdf"],
                                                selected="kde"
                                            ),
                                        ),
                                        ui.input_checkbox(
                                            "nn_log_scale",
                                            "Log Scale",
                                            value=False
                                        ),
                                        ui.input_checkbox(
                                            "nn_facet_plot",
                                            "Facet Plot",
                                            value=False
                                        ),
                                        ui.output_ui("nn_color_mapping_ui"),
                                    )
                                ),
                                ui.hr(),
                                ui.div(
                                    ui.input_checkbox(
                                        "nn_show_figure_config",
                                        "Show Figure Configuration",
                                        value=False
                                    ),
                                    ui.panel_conditional(
                                        "input.nn_show_figure_config",
                                        ui.row(
                                            ui.column(
                                                6,
                                                ui.input_numeric(
                                                    "nn_figure_width",
                                                    "Width",
                                                    value=10,
                                                    min=4,
                                                    max=20,
                                                    step=1
                                                ),
                                            ),
                                            ui.column(
                                                6,
                                                ui.input_numeric(
                                                    "nn_figure_height",
                                                    "Height",
                                                    value=6,
                                                    min=3,
                                                    max=15,
                                                    step=1
                                                ),
                                            )
                                        ),
                                        ui.row(
                                            ui.column(
                                                6,
                                                ui.input_numeric(
                                                    "nn_font_size",
                                                    "Font Size",
                                                    value=11,
                                                    min=8,
                                                    max=20
                                                ),
                                            ),
                                            ui.column(
                                                6,
                                                ui.input_numeric(
                                                    "nn_figure_dpi",
                                                    "DPI",
                                                    value=150,
                                                    min=72,
                                                    max=600,
                                                    step=25
                                                ),
                                            )
                                        ),
                                    )
                                ),
                                ui.div(
                                    ui.input_checkbox(
                                        "nn_show_axis_settings",
                                        "Show Axis Settings",
                                        value=False
                                    ),
                                    ui.panel_conditional(
                                        "input.nn_show_axis_settings",
                                        ui.input_numeric(
                                            "nn_x_axis_rotation",
                                            "X Axis Label Rotation (degrees)",
                                            value=0,
                                            min=-90,
                                            max=90
                                        ),
                                        ui.input_checkbox(
                                            "nn_shared_x_title",
                                            "Shared X Axis Title",
                                            value=True
                                        ),
                                        ui.input_numeric(
                                            "nn_x_title_fontsize",
                                            "X Axis Title Font Size",
                                            value=12,
                                            min=8,
                                            max=20
                                        ),
                                    )
                                ),
                                ui.br(),
                                ui.input_action_button(
                                    "go_nn_viz",
                                    "Generate Visualization",
                                    class_="btn-success w-100"
                                ),
                                ui.div(
                                    {"style": "padding-top: 10px;"},
                                    ui.output_ui("nn_stop_button_ui")
                                ),
                                ui.div(
                                    {"style": "padding-top: 10px;"},
                                    ui.output_ui("download_button_ui_nn")
                                ),
                                ui.div(
                                    {"style": "padding-top: 10px;"},
                                    ui.output_ui("download_nn_plot_button_ui")
                                ),
                            )
                        )
                    ),
                    ui.column(
                        9,
                        ui.div(
                            {
                                "class": "nn-plot-container",
                                "style": (
                                    "height: 85vh; overflow: auto; "
                                    "padding-left: 15px;"
                                )
                            },
                            ui.output_image(
                                "nn_visualization_plot",
                                width="100%",
                                height="auto"
                            )
                        )
                    )
                ),
            ),
        )
    )