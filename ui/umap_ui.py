from shiny import ui
from utils.accessibility import accessible_slider


def umap_ui():
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
                        accessible_slider(
                            "umap_slider_1",
                            "Point Size",
                            min_val=0.5,
                            max_val=10,
                            value=3,
                            step=0.1
                        ),
                        ui.input_action_button(
                            "go_umap1",
                            "Render Plot",
                            class_="btn-success",
                            style="width: 180px;"
                        ),
                        ui.div(
                            {"style": "padding-top: 10px;"},
                            ui.output_ui("umap1_stop_button_ui")
                        ),
                        ui.div(
                            {"style": "padding-top: 20px; padding-bottom: 20px;"},
                            ui.output_ui("download_umap1_plot_button_ui")
                        ),
                        ui.div(
                            {"style": "padding-top: 10px;"},
                            ui.output_image(
                                "spac_UMAP",
                                width="100%",
                                height="auto"
                            )
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
                        accessible_slider(
                            "umap_slider_2",
                            "Point Size",
                            min_val=0.5,
                            max_val=10,
                            value=3,
                            step=0.1
                        ),
                        ui.input_action_button(
                            "go_umap2",
                            "Render Plot",
                            class_="btn-success",
                            style="width: 180px;"
                        ),
                        ui.div(
                            {"style": "padding-top: 10px;"},
                            ui.output_ui("umap2_stop_button_ui")
                        ),
                        ui.div(
                            {"style": "padding-top: 20px; padding-bottom: 20px;"},
                            ui.output_ui("download_umap2_plot_button_ui")
                        ),
                        ui.div(
                            {"style": "padding-top: 10px;"},
                            ui.output_image(
                                "spac_UMAP2",
                                width="100%",
                                height="auto"
                            )
                        )
                    )
                )
            )
        )
    )