from shiny import ui


def data_input_ui():
    # 1. DATA INPUT PANEL -----------------------------------
    return ui.nav_panel("Data Input",
        # Add custom CSS to increase height of upload message/progress bar
        ui.tags.head(ui.tags.style("""
            .shiny-file-input-progress {
                height: 30px !important; 
                line-height: 30px !important;
            }
            .progress-bar {
                height: 30px !important;
                line-height: 30px !important;
                font-size: 16px !important;
            }
            /* Ensure text doesn't get cut off */
            .progress-bar span {
                white-space: nowrap;
                overflow: visible;
            }
            /* Style for the metric output text - larger but not bold */
            .metric-output {
                font-size: 18px;
            }
        """)),
        ui.card({"style": "width:100%;"},
            ui.column(
                12,
                ui.row(
                    ui.column(
                        6,
                        ui.card(
                            ui.div({
                                "style": "font-weight: bold; font-size: 30px;"
                            },
                            ui.p("SPAC Interactive Dashboard")),
                            ui.div({"style": "margin-bottom: 15px;"},
                                ui.input_file(
                                    "input_file", "Choose a file to upload:", 
                                    multiple=False, 
                                    width="100%"
                                ),
                                # ADDED: Helper text to show supported file types
                                ui.p(
                                    "Supported types: .pickle, .h5ad",
                                    style="font-size: 0.9em; font-style: italic; color: #6c757d;"
                                )
                            ),
                            ui.row(
                                ui.column(
                                    6,
                                    ui.card(
                                        ui.div({"class": "metric-output"},
                                            ui.output_text("print_rows")
                                        ),
                                        height="auto", class_="p-2 mb-2"
                                    ),
                                    ui.card(
                                        ui.div({"class": "metric-output"},
                                            ui.output_text("print_columns")
                                        ),
                                        height="auto", class_="p-2 mb-2"
                                    ),
                                    ui.card(
                                        ui.div({"class": "metric-output"},
                                            ui.output_text("print_obs_names")
                                        ),
                                        height="auto", class_="p-2 mb-2"
                                    )
                                ),
                                ui.column(
                                    6,
                                    ui.card(
                                        ui.div({"class": "metric-output"},
                                            ui.output_text("print_obsm_names")
                                        ),
                                        height="auto", class_="p-2 mb-2"
                                    ),
                                    ui.card(
                                        ui.div({"class": "metric-output"},
                                            ui.output_text(
                                                "print_layers_names"
                                            )
                                        ),
                                        height="auto", class_="p-2 mb-2"
                                    ),
                                    ui.card(
                                        ui.div({"class": "metric-output"},
                                            ui.output_text("print_uns_names")
                                        ),
                                        height="auto", class_="p-2 mb-2"
                                    )
                                )
                            ),
                            class_="mb-3"
                        ),
                        #SPAC TERMINOLOGY 
                        ui.card({"style": "width:100%;"},   
                            ui.column(
                                12,
                                ui.h2("SPAC Terminology"),
                                ui.p(
                                    "SPAC uses general terminology to " 
                                    "simplify technical terms from the "
                                    "AnnData object for less technical "
                                    "users. Here is a quick guide:"
                                ),
                                ui.tags.ul(
                                    ui.tags.li([
                                        ui.tags.b("Cells:"), 
                                        " Rows in the X matrix of AnnData."
                                    ]),
                                    ui.tags.li([
                                        ui.tags.b("Features:"), 
                                        " Columns in the X matrix of AnnData,"
                                        " representing gene expression or"
                                        " antibody intensity."
                                    ]),
                                    ui.tags.li([
                                        ui.tags.b("Tables:"), 
                                        " Originally called layers in "
                                        " AnnData, these represent transformed"
                                        " features."
                                    ]),
                                    ui.tags.li([
                                        ui.tags.b("Associated Tables:"), 
                                        " Corresponds to .obsm in AnnData and"
                                        " can store spatial coordinates, UMAP"
                                        " embeddings, etc."
                                    ]),
                                    ui.tags.li([
                                        ui.tags.b("Annotation:"), 
                                        " Corresponds to .obs in AnnData and" " can store cell phenotypes,"
                                        " experiment names, slide IDs, etc."
                                    ])
                                ),
                                ui.p(
                                    "For more in-depth explanations, "
                                    "visit our ",
                                    ui.a(
                                        "GitHub page", 
                                        href=
                                        "https://github.com/FNLCR-DMAP/spac_datamine/blob/main/CONTRIBUTING.md", 
                                        target="_blank"
                                    ),
                                    ".")
                            )
                        )
                    ),
                    ui.column(
                        6,
                        ui.card(
                            ui.input_checkbox(
                                "subset_select_check", 
                                "Subset Annotation", 
                                False
                            ),
                            ui.div(id="main-subset_anno_dropdown"),
                            ui.div(id="main-subset_label_dropdown"),
                            ui.input_action_button(
                                "go_subset", 
                                "Subset Data", 
                                class_="btn-success"
                            ),
                            ui.input_action_button(
                                "restore_data", 
                                "Restore Original Data", 
                                class_="btn-warning"
                            ),
                            ui.div({"class": "metric-output"},
                                ui.output_text("print_subset_history")
                            ),
                            class_="mb-3"
                        )
                    ),
                    ui.card(
                        {"style": "width:100%;"},
                        ui.h4("Annotation Summary with Top 10 Labels"),
                        ui.output_ui("annotation_labels_display")
                    )
                )
            )
        )
    )