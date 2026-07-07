"""
SPAC Shiny App - Effect Updates Server Components
This module contains server-side logic for updating UI elements
based on user input and shared data.
"""

from shiny import ui, render, reactive
from utils.data_processing import get_annotation_label_counts
import anndata as ad
import spac.data_utils


def effect_update_server(input, output, session, shared):
    @reactive.Effect
    def update_select_input_feat():
        choices = shared['var_names'].get()
        ui.update_select("h1_feat", choices=choices)
        ui.update_select("umap_feat", choices=choices)
        if choices:
            ui.update_select("bp_features", choices=choices)
            ui.update_selectize("bp_features", selected=choices[:2])

    @reactive.Effect
    def update_select_input_anno():
        choices = shared['obs_names'].get()
        ui.update_select("bp_anno", choices=choices)
        if choices:
            new_choices = choices + ["No Annotation"]
            ui.update_select("bp_anno", choices=new_choices)
        ui.update_select("h2_anno", choices=choices)
        ui.update_select("hm1_anno", choices=choices)

        if choices and len(choices) > 1:
            ui.update_select("sk1_anno1", choices=choices)
            ui.update_selectize("sk1_anno1", selected=choices[0])
            ui.update_select("sk1_anno2", choices=choices)
            ui.update_selectize("sk1_anno2", selected=choices[1])
            ui.update_select("rhm_anno1", choices=choices)
            ui.update_selectize("rhm_anno1", selected=choices[0])
            ui.update_select("rhm_anno2", choices=choices)
            ui.update_selectize("rhm_anno2", selected=choices[1])
        ui.update_select("spatial_anno", choices=choices)
        ui.update_select("nn_image_id", choices=["None"] + (choices or []))
        # rl_anno removed; no annotation selector needed for Ripley
        ui.update_select("region_select_rl", choices=choices)
        return

    @reactive.Effect
    def update_nearest_neighbor_choices():
        """Update nearest neighbor choices from spatial_distance columns."""
        # Get spatial_distance columns from shared state
        phenotype_choices = shared['spatial_distance_columns'].get()

        if phenotype_choices and len(phenotype_choices) > 0:
            # Update source label dropdown
            ui.update_select(
                "nn_source_label",
                choices=phenotype_choices,
            )

            # Update target label dropdown
            ui.update_selectize(
                "nn_target_label",
                choices=phenotype_choices,
            )

            # Set default source label to first available
            ui.update_select(
                "nn_source_label",
                selected=phenotype_choices[0],
            )
        else:
            # Clear choices if no spatial_distance data available
            ui.update_select("nn_source_label", choices=[])
            ui.update_selectize("nn_target_label", choices=[])
        return

    @reactive.Effect
    def update_select_input_layer():
        if shared['layers_names'].get():
            new_choices = shared['layers_names'].get() + ["Original"]
            ui.update_select("h1_layer", choices=new_choices)
            ui.update_select("bp_layer", choices=new_choices)
            ui.update_select("hm1_layer", choices=new_choices)
            ui.update_select("scatter_layer", choices=new_choices)
        return

    @reactive.Effect
    def update_rl_pairs():
        """Populate available Ripley phenotype pairs from
        adata.uns['ripley_l'].

        Choices are formatted as 'CENTER -> NEIGHBOR'.
        """
        adata = shared['adata_main'].get()
        if adata is None:
            ui.update_selectize("rl_pair", choices=[])
            return

        from utils.data_processing import get_rl_pairs

        choices = get_rl_pairs(adata)
        ui.update_selectize("rl_pair", choices=choices)
        if choices:
            ui.update_selectize("rl_pair", selected=choices[0])
        return

    @reactive.Effect
    def update_rl_region_and_slide_labels():
        """Populate region and slide label selectize controls based on the
        currently selected region/slide annotation columns.
        """
        adata = shared['adata_main'].get()
        if adata is None:
            ui.update_selectize("rl_region_labels", choices=[])
            return
        # Use utility to get label counts per annotation (handles None safely)
        label_counts = get_annotation_label_counts(adata)
        region_name = input.region_select_rl()

        if label_counts is None:
            return

        if region_name and region_name in label_counts:
            # Use helper to fetch sorted labels for the region annotation
            from utils.data_processing import get_annotation_top_labels

            top = get_annotation_top_labels(adata, region_name, top_n=None)
            region_label_list = [str(x[0]) for x in top]
            ui.update_selectize(
                "rl_region_labels",
                selected=region_label_list[:1] if region_label_list else [],
                choices=region_label_list,
            )

        # Slide stratification removed; no rl_slide_labels population needed

        return

    from ui.annotation_tables_ui import build_annotation_tables

    @output
    @render.ui
    def annotation_labels_display():
        """Retrieve and display top annotation labels in compact tables.

        Markup delegated to ui.annotation_tables_ui.build_annotation_tables.
        """
        adata = shared['adata_main'].get()
        annotation_counts = get_annotation_label_counts(adata)

        if not annotation_counts:
            return ui.tags.div("No annotations or data found.")

        return build_annotation_tables(annotation_counts)

    @reactive.Calc
    @render.text
    def print_obsm_names():
        obsm = shared['obsm_names'].get()
        if not obsm:
            return "Associated Tables: None"
        if obsm:
            if len(obsm) > 1:
                obsm_str = ", ".join(obsm)
            else:
                obsm_str = obsm[0] if obsm else ""
            return "Associated Tables: " + obsm_str
        else:
            return "Empty"

    # Initialize a flag to track dropdown creation
    subset_ui_initialized = reactive.Value(False)

    @reactive.Effect
    def subset_reactivity():
        # Get input values
        _ = ad.AnnData(obs=shared['obs_data'].get())
        annotations = shared['obs_names'].get()
        btn = input.subset_select_check()

        # Access the flag value
        ui_initialized = subset_ui_initialized.get()

        if btn and not ui_initialized:
            # Insert annotation dropdown
            anno_dropdown = ui.input_select(
                "subset_anno_select",
                "Select Annotation to subset",
                choices=annotations,
            )
            ui.insert_ui(
                ui.div({"id": "inserted-subset_anno_dropdown"}, anno_dropdown),
                selector="#main-subset_anno_dropdown",
                where="beforeEnd",
            )

            # Insert label dropdown
            label_dropdown = ui.input_selectize(
                "subset_label_select",
                "Select a Label",
                choices=[],
                selected=[],
                multiple=True,
            )
            ui.insert_ui(
                ui.div(
                    {"id": "inserted-subset_label_dropdown"},
                    label_dropdown,
                ),
                selector="#main-subset_label_dropdown",
                where="beforeEnd",
            )

            # Mark the dropdowns as initialized
            subset_ui_initialized.set(True)

        elif not btn and ui_initialized:
            # Remove dropdowns if the checkbox is unchecked
            ui.remove_ui("#inserted-subset_anno_dropdown")
            ui.remove_ui("#inserted-subset_label_dropdown")

            # Reset the flag
            subset_ui_initialized.set(False)

    # Create a reactive trigger for label updates
    label_update_trigger = reactive.Value(0)

    @reactive.effect
    def update_subset_labels():
        # This effect depends on both the selected annotation and the trigger
        _ = label_update_trigger.get()  # Add trigger dependency
        adata = shared['adata_main'].get()
        selected_anno = input.subset_anno_select()

        if adata and selected_anno:
            labels = adata.obs[selected_anno].unique().tolist()
            print(f"Updating labels for {selected_anno}: {labels}")
            ui.update_selectize("subset_label_select", choices=labels)

    @reactive.effect
    @reactive.event(input.go_subset, ignore_none=True)
    def subset_stratification():
        adata = shared['adata_main'].get()
        if adata:
            annotation = input.subset_anno_select()
            labels = list(input.subset_label_select())

            # Perform the subsetting
            adata_subset = spac.data_utils.select_values(
                adata, annotation=annotation, values=labels
            )

            # Update the main data
            shared['adata_main'].set(adata_subset)

            # Increment the label update trigger to force update
        label_update_trigger.set(label_update_trigger.get() + 1)

    # Reactive variable to store subset history as a simple string
    subset_history = reactive.Value("")

    @reactive.effect
    @reactive.event(input.go_subset, ignore_none=True)
    def track_subset():
        """Append the current annotation and selected labels to history."""
        annotation = input.subset_anno_select()
        labels = input.subset_label_select()

        # If both annotation and labels are valid, add them to history
        if annotation and labels:
            # Create a string like "Annotation:label1,label2"
            new_entry = f"{annotation}: {','.join(labels)}"

            # Update the subset history string by appending the new entry
            current_history = subset_history.get()
            if current_history:
                # If there's already history, append with a separator
                subset_history.set(f"{current_history} -> {new_entry}")
            else:  # Otherwise, just set the first entry
                subset_history.set(new_entry)

    @output
    @render.text
    def print_subset_history():
        """
        Render the subset history as plain text for the UI.
        """
        history = subset_history.get()
        return history if history else "No subsets have been made yet."

    # Reactive variable to store the master copy of the inputted adata object
    adata_master = reactive.Value(None)

    @reactive.Effect
    def store_master_copy():
        """Store a master copy of the adata object when first loaded."""
        adata = shared['adata_main'].get()
        if adata and adata_master.get() is None:
            # Make a copy of the adata object and store it as the master copy
            adata_master.set(adata.copy())

    # Add a button to the UI for restoring the data
    ui.input_action_button(
        "restore_data",
        "Restore Original Data",
        class_="btn-warning",
    )

    @reactive.effect
    @reactive.event(input.restore_data, ignore_none=True)
    def restore_to_master():
        """
        Restore adata_main to the master copy stored in adata_master.
        """
        master_data = adata_master.get()
        if master_data:
            # Set adata_main to a copy of the master data to ensure
            # independence
            shared['adata_main'].set(master_data.copy())

            # Clear the subset history since the data is restored
            subset_history.set("")

