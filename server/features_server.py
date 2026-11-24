"""
Feature histogram visualization module for SPAC Shiny application.

This module handles the server-side logic for generating feature histograms
through the standard SPAC Histogram template workflow.
"""

import tempfile

import matplotlib.pyplot as plt
from shiny import ui, render, reactive

# Import plot display utilities.
from utils.plot_utils import fig_to_png_bytes
# Import template wrapper utilities.
from utils.template_wrapper import (
    register_memory_object,
    unregister_memory_object,
)
# Import the specific template
from spac.templates.histogram_template import run_from_json


def features_server(input, output, session, shared):
    """
    Server logic for the feature histogram visualization.

    Parameters
    ----------
    input : shiny.session.Inputs
        Shiny input object
    output : shiny.session.Outputs
        Shiny output object
    session : shiny.session.Session
        Shiny session object
    shared : dict
        Shared reactive values across server modules
    """

    @reactive.calc
    def get_adata():
        """Get the main AnnData object from shared state."""
        return shared['adata_main'].get()

    @reactive.calc
    def get_features_inputs():
        """Return normalized semantic state from the Features controls."""

        group_by_enabled = input.h1_group_by_check()
        group_by = input.h1_anno() if group_by_enabled else None
        together = input.h1_together_check() if group_by_enabled else False
        multiple = input.h1_together_drop() if together else "stack"

        return {
            "feature": input.h1_feat(),
            "layer": (
                input.h1_layer() if input.h1_layer() != "Original" else None
            ),
            "x_log_scale": input.h1_log_x(),
            "y_log_scale": input.h1_log_y(),
            "group_by": group_by,
            "together": together,
            "multiple": multiple,
            "x_axis_label_rotation": input.feat_slider(),
        }

    def build_template_params(input_values, virtual_path):
        """Convert normalized semantic values to template parameters.

        The template supplies defaults for controls not exposed by the UI.
        """

        return {
            "Upstream_Analysis": virtual_path,
            "Feature": input_values["feature"] or "None",
            "Table_": input_values["layer"] or "Original",
            "Group_by": input_values["group_by"] or "None",
            "Together": input_values["together"],
            "Take_X_Log": input_values["x_log_scale"],
            "Take_Y_log": input_values["y_log_scale"],
            "Multiple": input_values["multiple"],
            "X_Axis_Label_Rotation": input_values["x_axis_label_rotation"],
            "Plot_By": "Feature",  # Features tab invariant
            "Facet": False,  # Facet UI deferred
        }

    @output
    @render.image(delete_file=True)
    @reactive.event(input.go_h1, ignore_none=True)
    def spac_Histogram_1():
        """Generate a feature histogram through the SPAC template."""

        adata = get_adata()
        if adata is None:
            return None
        features_inputs = get_features_inputs()

        virtual_path = register_memory_object(adata)
        try:
            params = build_template_params(features_inputs, virtual_path)
            fig, df = run_from_json(
                json_path=params,
                save_to_disk=False,
                show_plot=False,
            )
        finally:
            unregister_memory_object(virtual_path)

        try:
            png_bytes = fig_to_png_bytes(fig)
        finally:
            plt.close(fig)
        with tempfile.NamedTemporaryFile(
            prefix="spac_feature_histogram_",
            suffix=".png",
            delete=False,
        ) as temp_file:
            temp_file.write(png_bytes)
            temp_path = temp_file.name

        shared['df_histogram1'].set(df)
        return {
            "src": temp_path,
            "width": "100%",
            "height": "100%",
            "style": "object-fit: contain;",
            "alt": "Feature histogram",
        }

    @render.download(filename="features_histogram_data.csv")
    def download_histogram1_df():
        df = shared['df_histogram1'].get()
        if df is not None:
            csv_string = df.to_csv(index=False)
            csv_bytes = csv_string.encode("utf-8")
            return csv_bytes, "text/csv"
        return None

    @render.ui
    @reactive.event(input.go_h1, ignore_none=True)
    def download_histogram1_button_ui():
        if shared['df_histogram1'].get() is not None:
            return ui.download_button(
                "download_histogram1_df", 
                "Download Data", 
                class_="btn-warning"
            )
        return None
