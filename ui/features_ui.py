"""
Features histogram UI module for SPAC Shiny application.

This module provides the user interface components for visualizing feature
distributions using the histogram template functionality.
"""

from shiny import ui
from utils.accessibility import accessible_slider


def features_ui():
    """
    Create the Features histogram panel UI.

    Returns
    -------
    shiny.ui.NavPanel
        UI nav panel for the features histogram visualization
    """
    return ui.nav_panel(
        "Features",
        ui.card(
            {"style": "width:100%;"},
            ui.column(
                12,
                ui.row(

                    # ── Left Controls Column ──────────────────────────────
                    ui.column(
                        2,

                        # ── Core Parameters ───────────────────────────────
                        ui.input_select(
                            "h1_feat",
                            "Select a Feature",
                            choices=[]
                        ),
                        ui.input_select(
                            "h1_layer",
                            "Select a Table",
                            choices=[],
                            selected=["Original"]
                        ),
                        ui.input_checkbox(
                            "h1_group_by_check",
                            "Group By",
                            value=False
                        ),
                             # ▼▼▼ NEW: shown only when Group By is checked
                        ui.panel_conditional(
                            "input.h1_group_by_check",
                            ui.div(id="main-h1_dropdown"),
                            ui.div(
                                {"style": "margin-top: 6px;"},
                                ui.input_checkbox(
                                    "h1_facet",
                                    "Facet Plot",
                                    value=False
                                ),
                            ),
                            ui.div(id="main-h1_check"),
                            ui.div(id="main-h1_together_drop"),
                        ),
                        ui.input_checkbox(
                            "h1_log_x",
                            "Log X-axis",
                            value=False
                        ),
                        ui.input_checkbox(
                            "h1_log_y",
                            "Log Y-axis",
                            value=False
                        ),

                        # Dynamic Group By UI injection targets
                        ui.div(id="main-h1_dropdown"),
                        ui.div(id="main-h1_check"),
                        ui.div(id="main-h1_together_drop"),

                        ui.hr(),

                        # ── Bins Settings ─────────────────────────────────
                        ui.div(
                            ui.input_checkbox(
                                "h1_show_bins",
                                "Show Bins Settings",
                                value=False
                            ),
                            ui.panel_conditional(
                                "input.h1_show_bins",
                                ui.input_radio_buttons(
                                    "h1_bins_type",
                                    "Bins Input Type",
                                    choices={
                                        "auto":   "Auto — Determined automatically",
                                        "number": "Number — Set number of bins",
                                        #"list":   "List — Set custom bin edges",
                                    },
                                    selected="auto"
                                ),
                                ui.panel_conditional(
                                    "input.h1_bins_type === 'number'",
                                    ui.input_numeric(
                                        "h1_bins_number",
                                        "Number of Bins",
                                        value=10,
                                        min=1,
                                        step=1
                                    ),
                                ),
                                ui.panel_conditional(
                                    "input.h1_bins_type === 'list'",
                                    ui.input_text(
                                        "h1_bins_list",
                                        "Bin Edges (comma-separated)",
                                        placeholder="e.g. 0, 1, 2, 3"
                                    ),
                                    ui.tags.small(
                                        {"style": "color: #6c757d;"},
                                        "Creates bins: [0,1), [1,2), [2,3]"
                                    ),
                                ),
                            ),
                        ),

                        ui.hr(),

                        # ── Stat Settings ─────────────────────────────────
                        ui.div(
                            ui.input_checkbox(
                                "h1_show_stat",
                                "Show Stat Settings",
                                value=False
                            ),
                            ui.panel_conditional(
                                "input.h1_show_stat",
                                ui.input_select(
                                    "h1_stat",
                                    "Statistical Transformation",
                                    choices={
                                        "count":       "Count — Number of observations per bin",
                                        "frequency":   "Frequency — Count divided by bin width",
                                        "density":     "Density — Total area normalized to 1",
                                        "probability": "Probability — Bar height as observation probability",
                                    },
                                    selected="count"
                                ),
                            ),
                        ),

                        ui.hr(),

                        # ── Element Settings ──────────────────────────────
                        ui.div(
                            ui.input_checkbox(
                                "h1_show_element",
                                "Show Element Settings",
                                value=False
                            ),
                            ui.panel_conditional(
                                "input.h1_show_element",
                                ui.input_select(
                                    "h1_element",
                                    "Visual Representation of Bins",
                                    choices={
                                        "bars": "Bars — Standard bar-style histogram (default)",
                                        "step": "Step — Step line plot without bars",
                                        "poly": "Poly — Polygon with x-axis as bottom edge",
                                    },
                                    selected="bars"
                                ),
                            ),
                        ),

                        ui.hr(),

                        # ── Axis Settings ─────────────────────────────────
                        ui.div(
                            ui.input_checkbox(
                                "h1_show_axis_settings",
                                "Show Axis Settings",
                                value=False
                            ),
                            ui.panel_conditional(
                                "input.h1_show_axis_settings",
                                accessible_slider(
                                    "feat_slider",
                                    "Rotate X-axis Labels (degrees)",
                                    min_val=0,
                                    max_val=90,
                                    value=0,
                                    step=1
                                ),
                                ui.input_numeric(
                                    "feat_slider_num",
                                    "Or Type a Value (degrees)",
                                    value=0,
                                    min=0,
                                    max=90,
                                    step=1
                                ),
                            ),
                        ),

                        ui.hr(),

                        # ── Figure Configuration ──────────────────────────
                        ui.div(
                            ui.input_checkbox(
                                "h1_show_figure_config",
                                "Show Figure Configuration",
                                value=False
                            ),
                            ui.panel_conditional(
                                "input.h1_show_figure_config",
                                ui.row(
                                    ui.column(
                                        6,
                                        ui.input_numeric(
                                            "h1_figure_width",
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
                                            "h1_figure_height",
                                            "Height",
                                            value=6,
                                            min=3,
                                            max=15,
                                            step=1
                                        ),
                                    ),
                                ),
                                ui.row(
                                    ui.column(
                                        6,
                                        ui.input_numeric(
                                            "h1_font_size",
                                            "Font Size",
                                            value=11,
                                            min=8,
                                            max=20,
                                            step=1
                                        ),
                                    ),
                                    ui.column(
                                        6,
                                        ui.input_numeric(
                                            "h1_figure_dpi",
                                            "DPI",
                                            value=150,
                                            min=72,
                                            max=600,
                                            step=25
                                        ),
                                    ),
                                ),
                            ),
                        ),

                        ui.hr(),

                        # ── Action Button & Download ──────────────────────
                        ui.input_action_button(
                            "go_h1",
                            "Render Plot",
                            class_="btn-success"
                        ),
                        ui.div(
                            {"style": "padding-top: 20px;"},
                            ui.output_ui("download_histogram1_button_ui")
                        ),
                    ),

                    # ── Right Plot Column ─────────────────────────────────
                    ui.column(
                        10,
                        ui.div(
                            {"style": "padding-bottom: 100px;"},
                            ui.output_plot(
                                "spac_Histogram_1",
                                width="100%",
                                # height="60vh"
                                height="150vh",
                            )
                        )
                    )
                )
            )
        )
    )