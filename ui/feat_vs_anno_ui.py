"""
Feature vs Annotation heatmap visualization UI module for SPAC Shiny
application.

This module defines the user interface components for generating and
customizing heatmaps using the hierarchical_heatmap function.
"""

from shiny import ui


def feat_vs_anno_ui():
    """
    Create the feature vs annotation heatmap visualization UI.

    Returns
    -------
    shiny.ui.NavPanel
        UI components for the feature vs annotation heatmap feature
    """
    return ui.nav_panel(
        "Feat. Vs Anno.",
        # Custom CSS for improved layout
        ui.tags.style("""
            .hm1-controls-panel {
                background-color: #f8f9fa;
                border-radius: 8px;
                padding: 10px;
            }
            .hm1-plot-container {
                background-color: white;
                border-radius: 8px;
                border: 1px solid #dee2e6;
                padding: 10px;
            }
            .hm1-controls-panel .form-group {
                margin-bottom: 12px;
            }
            .hm1-controls-panel h4, .hm1-controls-panel h5 {
                margin-bottom: 10px;
                margin-top: 15px;
            }
            .hm1-controls-panel h4:first-child {
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
                                "class": "hm1-controls-panel",
                                "style": (
                                    "height: 85vh; overflow-y: auto; "
                                    "padding-right: 10px; "
                                    "border-right: 1px solid #dee2e6;"
                                )
                            },
                            ui.h4("Core Parameters",
                                class_="accessible-heading"),

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

                            ui.hr(),

                            # Plot configuration in expandable section
                            ui.div(
                                ui.input_checkbox(
                                    "hm1_show_plot_config",
                                    "Show Plot Configuration",
                                    value=False
                                ),
                                ui.panel_conditional(
                                    "input.hm1_show_plot_config",
                                    ui.input_select(
                                        "hm1_cmap",
                                        "Select Color Map",
                                        choices=[
                                            "viridis", "plasma", "inferno",
                                            "magma", "cividis", "coolwarm",
                                            "RdYlBu", "Spectral", "PiYG",
                                            "PRGn"
                                        ],
                                        selected="viridis"
                                    ),
                                    ui.input_checkbox(
                                        "hm1_dendogram",
                                        "Include Dendrogram",
                                        False
                                    ),
                                    ui.panel_conditional(
                                        "input.hm1_dendogram",
                                        ui.input_checkbox(
                                            "hm1_feat_dendro",
                                            "Feature Cluster",
                                            value=False
                                        ),
                                        ui.input_checkbox(
                                            "hm1_anno_dendro",
                                            "Annotation Cluster",
                                            value=False
                                        ),
                                    ),
                                    ui.div(id="main-hm1_min_num"),
                                    ui.div(id="main-hm1_max_num"),
                                ),
                            ),

                            ui.hr(),

                            # Axis settings in expandable section
                            ui.div(
                                ui.input_checkbox(
                                    "hm1_show_axis_settings",
                                    "Show Axis Settings",
                                    value=False
                                ),
                                ui.panel_conditional(
                                    "input.hm1_show_axis_settings",
                                    ui.input_numeric(
                                        "hm1_x_label_rotation",
                                        "Rotate X Axis Labels (degrees)",
                                        min=0,
                                        max=90,
                                        value=50
                                    ),
                                    ui.input_numeric(
                                        "hm1_y_label_rotation",
                                        "Rotate Y Axis Labels (degrees)",
                                        min=0,
                                        max=90,
                                        value=25
                                    ),
                                    ui.input_numeric(
                                        "hm1_axis_label_fontsize",
                                        "Axis Label Font Size",
                                        min=3,
                                        max=24,
                                        value=10
                                    ),
                                    ui.input_checkbox(
                                        "hm1_enable_abbreviation",
                                        "Abbreviate Axis Labels",
                                        value=False
                                    ),
                                    ui.panel_conditional(
                                        "input.hm1_enable_abbreviation",
                                        ui.input_numeric(
                                            "hm1_label_char_limit",
                                            "Axis Label Character Limit",
                                            min=3,
                                            max=20,
                                            value=6,
                                        ),
                                    ),
                                ),
                            ),

                            ui.br(),
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
                                ui.output_ui("download_button_ui_hm1")
                            ),
                            ui.div(
                                {"style": "padding-top: 10px;"},
                                ui.output_ui("download_heatmap_plot_button_ui")
                            ),
                        ),
                    ),
                    ui.column(
                        9,
                        ui.div(
                            {
                                "class": "hm1-plot-container",
                                "style": (
                                    "height: 85vh; overflow: auto; "
                                    "padding-left: 15px;"
                                ),
                            },
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
