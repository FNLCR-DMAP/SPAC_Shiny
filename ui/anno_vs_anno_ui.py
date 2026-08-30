from shiny import ui
from shinywidgets import output_widget


def anno_vs_anno_ui():
    # 6. ANNO. VS ANNO. (Sankey, Relational Heatmap) ---------
    return ui.nav_panel(
        "Anno. Vs Anno.",
        # Overview Section
        ui.card(
            {"style": "width:100%; margin-bottom: 20px; "
             "background-color: #f8f9fa;"},
            ui.card_header(
                ui.h3(
                    "Annotation Relationship Analysis",
                    style="margin: 0; color: #2c3e50; text-align: center;"
                )
            ),
            ui.p(
                "This page provides two complementary visualizations for "
                "exploring relationships between annotations:",
                style="text-align: center; margin-bottom: 15px; "
                      "font-size: 1.1em;"
            ),
            ui.row(
                    ui.column(
                        6,
                        ui.div(
                            {"style": "text-align: center; padding: 15px; "
                             "border: 2px solid #28a745; "
                             "border-radius: 8px; margin: 5px;"},
                            ui.h5(
                                "📊 Sankey Diagram",
                                style="color: #155724; margin-bottom: 10px; "
                                      "font-weight: 600;"
                            ),
                            ui.p(
                                "Interactive flow visualization showing "
                                "connections between annotation categories",
                                style="margin: 0; color: #343a40;"
                            )
                        )
                    ),
                    ui.column(
                        6,
                        ui.div(
                            {"style": "text-align: center; padding: 15px; "
                             "border: 2px solid #ffc107; "
                             "border-radius: 8px; margin: 5px;"},
                            ui.h5(
                                "🔥 Relational Heatmap",
                                style="color: #721c24; margin-bottom: 10px; "
                                      "font-weight: 600;"
                            ),
                            ui.p(
                                "Quantitative analysis with interactive "
                                "heatmap and downloadable data",
                                style="margin: 0; color: #343a40;"
                            )
                        )
                    )
                )
        ),
        # Sankey Plot Section
        ui.card(
            {"style": "width:100%;"},
            ui.card_header(
                ui.h4(
                    "Sankey Diagram",
                    style="margin: 0; color: #2c3e50;"
                ),
                ui.p(
                    "Visualize the flow and relationships between source "
                    "and target annotations using an interactive Sankey "
                    "diagram.",
                    style="margin: 5px 0 0 0; color: #6c757d; "
                          "font-size: 0.9em;"
                )
            ),
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
                        ui.input_action_button(
                            "go_sk1",
                            "Generate Sankey Plot",
                            class_="btn-success"
                        ),
                        ui.div(
                            {"style": "padding-top: 20px;"},
                            ui.output_ui("download_button_ui_sankey")
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
        # Relational Heatmap Section
        ui.card(
            {"style": "width:100%;"},
            ui.card_header(
                ui.h4(
                    "Relational Heatmap",
                    style="margin: 0; color: #2c3e50;"
                ),
                ui.p(
                    "Explore quantitative relationships between annotations "
                    "with an interactive heatmap and downloadable data.",
                    style="margin: 5px 0 0 0; color: #6c757d; "
                          "font-size: 0.9em;"
                )
            ),
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
                        ui.input_action_button(
                            "go_rhm1",
                            "Generate Heatmap",
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
