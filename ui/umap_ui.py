from shiny import ui

def umap_ui():
    # 8. UMAP PANEL ------------------------------------------
    return ui.nav_panel("UMAP",
        ui.card({"style": "width:100%;"},
            ui.column(
                12,
                ui.row(
                    ui.column(
                        6,
                        ui.input_radio_buttons(
                            "umap_rb", 
                            "Choose one:", 
                            ["Annotation", "Feature"]
                        ),
                        ui.input_select(
                            "plottype", 
                            "Select a plot type", 
                            choices=["umap", "pca", "tsne"]
                        ),
                        ui.div(id="main-ump_rb_dropdown_anno"),
                        ui.div(id="main-ump_rb_dropdown_feat"),
                        ui.div(id="main-ump_table_dropdown_feat"),
                        ui.input_slider(
                            "umap_slider_1", 
                            "Point Size", 
                            min=.5, 
                            max=10, 
                            value=3
                        ),
                        # Added...
                        ui.input_slider(
                            "umap_font_size_1",
                            "Font Size",
                            min=5,
                            max=30,
                            value=12
                        ),
                        ui.input_action_button(
                            "go_umap1", 
                            "Render Plot", 
                            class_="btn-success"
                        ),
                        ui.output_plot(
                            "spac_UMAP", 
                            width="100%", 
                            height="80vh"
                        )
                    ),
                    ui.column(
                        6,
                        ui.input_radio_buttons(
                            "umap_rb2", 
                            "Choose one:", 
                            ["Annotation", "Feature"]
                        ),
                        ui.input_select(
                            "plottype2", 
                            "Select a plot type", 
                            choices=["umap", "pca", "tsne"]
                        ),
                        ui.div(id="main-ump_rb_dropdown_anno2"),
                        ui.div(id="main-ump_rb_dropdown_feat2"),
                        ui.div(id="main-ump_table_dropdown_feat2"),
                        ui.input_slider(
                            "umap_slider_2", 
                            "Point Size", 
                            min=.5, 
                            max=10,
                            value=3
                        ),
                        # Added missing font size slider for the second plot
                        ui.input_slider(
                            "umap_font_size_2",
                            "Font Size",
                            min=5,
                            max=30,
                            value=12
                        ),
                        ui.input_action_button(
                            "go_umap2", 
                            "Render Plot", 
                            class_="btn-success"
                        ),
                        ui.output_plot(
                            "spac_UMAP2", 
                            width="100%", 
                            height="80vh"
                        )
                    )
                )
            )
        )
    )