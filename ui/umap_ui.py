from shiny import ui
from utils.accessibility import accessible_slider


def _figure_legend_panel(suffix: str):
    """
    Collapsible figure and legend controls. suffix '' for column 1, '2' for column 2.
    """
    show_id = f"umap_show_figure_config{suffix}"
    cond = f"input.{show_id}"
    mode_is_annotation = f"input.umap_rb{suffix} == 'Annotation'"
    return ui.div(
        ui.input_checkbox(
            show_id,
            "Show figure & legend settings",
            value=False,
        ),
        ui.panel_conditional(
            cond,
            ui.row(
                ui.column(
                    6,
                    ui.input_numeric(
                        f"umap_fig_width{suffix}",
                        "Figure width",
                        value=10,
                        min=4,
                        max=30,
                    ),
                ),
                ui.column(
                    6,
                    ui.input_numeric(
                        f"umap_fig_height{suffix}",
                        "Figure height",
                        value=8,
                        min=4,
                        max=30,
                    ),
                ),
            ),
            ui.row(
                ui.column(
                    6,
                    ui.input_numeric(
                        f"umap_fig_dpi{suffix}",
                        "DPI",
                        value=110,
                        min=72,
                        max=600,
                        step=20,
                    ),
                ),
                ui.column(
                    6,
                    ui.input_numeric(
                        f"umap_font_size{suffix}",
                        "Font size",
                        value=8,
                        min=6,
                        max=24,
                    ),
                ),
            ),
            ui.row(
                ui.column(
                    6,
                    ui.input_numeric(
                        f"umap_legend_font_size{suffix}",
                        "Legend font size",
                        value=8,
                        min=6,
                        max=32,
                    ),
                ),
                ui.column(
                    6,
                    ui.panel_conditional(
                        mode_is_annotation,
                        ui.input_numeric(
                            f"umap_legend_marker_scale{suffix}",
                            "Legend marker scale",
                            value=1.5,
                            min=0.5,
                            max=20,
                            step=0.5,
                        ),
                    ),
                ),
            ),
        ),
    )


def _feature_value_range_panel(suffix: str):
    """Min/max for feature coloring; inputs always present for stable server reads."""
    mode_is_feature = f"input.umap_rb{suffix} == 'Feature'"
    return ui.panel_conditional(
        mode_is_feature,
        ui.div(
            ui.input_checkbox(
                f"umap_use_val_range{suffix}",
                "Set feature color min/max",
                value=False,
            ),
            ui.panel_conditional(
                f"input.umap_use_val_range{suffix}",
                ui.row(
                    ui.column(
                        6,
                        ui.input_numeric(
                            f"umap_val_min{suffix}",
                            "Value min",
                            value=0.0,
                        ),
                    ),
                    ui.column(
                        6,
                        ui.input_numeric(
                            f"umap_val_max{suffix}",
                            "Value max",
                            value=1.0,
                        ),
                    ),
                ),
            ),
        ),
    )


def umap_ui():
    return ui.nav_panel(
        "UMAP",
        ui.card(
            {"style": "width:100%;"},
            ui.column(
                12,
                ui.row(
                    ui.column(
                        6,
                        ui.input_radio_buttons(
                            "umap_rb",
                            "Choose one:",
                            ["Annotation", "Feature"],
                        ),
                        ui.input_select(
                            "plottype",
                            "Select a plot type",
                            choices=["umap", "pca", "tsne"],
                        ),
                        ui.div(id="main-ump_rb_dropdown_anno"),
                        ui.div(id="main-ump_rb_dropdown_feat"),
                        ui.div(id="main-ump_table_dropdown_feat"),
                        accessible_slider(
                            "umap_slider_1",
                            "Point Size",
                            min_val=1,
                            max_val=10,
                            value=3,
                            step=1,
                        ),
                        _feature_value_range_panel(""),
                        _figure_legend_panel(""),
                        ui.input_action_button(
                            "go_umap1",
                            "Render Plot",
                            class_="btn-success",
                        ),
                        ui.output_plot(
                            "spac_UMAP",
                            width="100%",
                            height="80vh",
                        ),
                    ),
                    ui.column(
                        6,
                        ui.input_radio_buttons(
                            "umap_rb2",
                            "Choose one:",
                            ["Annotation", "Feature"],
                        ),
                        ui.input_select(
                            "plottype2",
                            "Select a plot type",
                            choices=["umap", "pca", "tsne"],
                        ),
                        ui.div(id="main-ump_rb_dropdown_anno2"),
                        ui.div(id="main-ump_rb_dropdown_feat2"),
                        ui.div(id="main-ump_table_dropdown_feat2"),
                        accessible_slider(
                            "umap_slider_2",
                            "Point Size",
                            min_val=1,
                            max_val=10,
                            value=3,
                            step=1,
                        ),
                        _feature_value_range_panel("2"),
                        _figure_legend_panel("2"),
                        ui.input_action_button(
                            "go_umap2",
                            "Render Plot",
                            class_="btn-success",
                        ),
                        ui.output_plot(
                            "spac_UMAP2",
                            width="100%",
                            height="80vh",
                        ),
                    ),
                ),
            ),
        ),
    )
