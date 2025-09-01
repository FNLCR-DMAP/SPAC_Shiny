from shiny import ui
from shinywidgets import output_widget

# TODO: Stratify by region + slide at same time breaks
# TODO: Also, straitfy by slide annotation not same as current annotation breaks
# TODO: Also, selecting feature then trying to stratify breaks
def spatial_ui():
    # 7. SPATIAL PANEL ---------------------------------------
    return ui.nav_panel(
        "Spatial",
        ui.card(
            {"style": "width:100%;"},
            ui.column(
                12,
                ui.row(
                    ui.column(
                        2,
                        ui.input_radio_buttons(
                            "spatial_rb",
                            "Color by:",
                            ["Annotation", "Feature"]
                        ),
                        ui.div(id="main-spatial_dropdown_anno"),
                        ui.div(id="main-spatial_dropdown_feat"),
                        ui.div(id="main-spatial_table_dropdown_feat"),
                        ui.input_slider(
                            "spatial_slider",
                            "Point Size",
                            min=2,
                            max=10,
                            value=3
                        ),
                        # Added...
                        ui.input_slider(
                            "spatial_font_size",
                            "Font Size",
                            min=5,
                            max=30,
                            value=12
                        ),
                        ui.input_checkbox(
                            "slide_select_check",
                            "Stratify by Slide",
                            False
                        ),
                        ui.div(id="main-slide_dropdown"),
                        ui.div(id="main-label_dropdown"),
                        ui.input_checkbox(
                            "region_select_check",
                            "Stratify by Region",
                            False
                        ),
                        ui.div(id="main-region_dropdown"),
                        ui.div(id="main-region_label_select_dropdown"),
                        ui.input_action_button(
                            "go_sp1",
                            "Render Plot",
                            class_="btn-success"
                        )
                    ),
                    ui.column(
                        10,
                        ui.div(
                            {"style": "padding-bottom: 20px;"},
                            output_widget("spac_Spatial"),
                            style="width:100%; height:80vh;"
                        )
                    )
                )
            )
        )
    )