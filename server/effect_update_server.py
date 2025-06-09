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
        ui.update_select("bp_features", choices=choices)


    @reactive.Effect
    def update_select_input_anno():
        choices = shared['obs_names'].get()
        ui.update_select("bp_anno", choices=choices)
        ui.update_select("h2_anno", choices=choices)
        ui.update_select("hm1_anno", choices=choices)
        ui.update_select("sk1_anno1", choices=choices)
        ui.update_select("sk1_anno2", choices=choices)
        ui.update_select("rhm_anno1", choices=choices)
        ui.update_select("rhm_anno2", choices=choices)
        ui.update_select("spatial_anno", choices=choices)
        ui.update_select("nn_anno", choices=choices)   
        ui.update_select("nn_strat_select", choices=choices)
        ui.update_select("rl_anno", choices=choices)
        ui.update_select("region_select_rl", choices=choices)
        ui.update_select("slide_select_rl", choices=choices)
        return

    @reactive.Effect
    def update_select_input_layer():
        if shared['layers_names'].get() is not None:
            new_choices = shared['layers_names'].get() + ["Original"]
            ui.update_select("h1_layer", choices=new_choices)
            ui.update_select("bp_layer", choices=new_choices)
            ui.update_select("hm1_layer", choices=new_choices)
            ui.update_select("scatter_layer", choices=new_choices)
        return
        
    @reactive.Effect
    def update_select_input_anno_bp():
        if shared['obs_names'].get() is not None:
            new_choices = shared['obs_names'].get() + ["No Annotation"]
            ui.update_select("bp_anno", choices=new_choices)


    @reactive.Effect
    def update_boxplot_selectize():
        selected_names=shared['var_names'].get()
        if selected_names is not None:
            ui.update_selectize("bp_features", selected=selected_names[:2])
            return

    @reactive.Effect
    def update_relational_select():
        selected_names=shared['obs_names'].get()
        if selected_names is not None and len(selected_names) > 1:
            ui.update_selectize("rhm_anno1", selected=selected_names[0])
            ui.update_selectize("rhm_anno2", selected=selected_names[1])
        return
    
    @reactive.Effect
    def update_rl_selectize():
        adata = shared['adata_main'].get()  # Reactive AnnData
        label_counts =  get_annotation_label_counts(adata)
        anno_name = input.rl_anno() # Selected annotation column
        region_name = input.region_select_rl() # Selected region column
        slide_name = input.slide_select_rl() # Selected slide column

        if anno_name is None or anno_name not in label_counts:
            return
        if region_name is None or region_name not in label_counts:
            return
        if slide_name is None or slide_name not in label_counts:
            return
        pheno_label_list = list(map(str, label_counts[anno_name].keys()))
        region_label_list =list(map(str, label_counts[region_name].keys())) 
        slide_label_list =list(map(str, label_counts[slide_name].keys()))

        if label_counts is not None:
            ui.update_selectize("rl_label", selected=pheno_label_list[:2], choices=pheno_label_list)
            ui.update_selectize("rl_region_labels", selected=region_label_list[0], choices=region_label_list)
            ui.update_selectize("rl_slide_labels", choices=slide_label_list)
            return

    @output
    @render.ui
    def annotation_labels_display():
        """
        1) Retrieve ALL annotations (via get_annotation_label_counts).
        2) Within each annotation, keep only the top 5 labels (sorted by count).
        3) Display them in separate cards, each card listing the top 5 labels.
        """
        adata = shared['adata_main'].get()  # your reactive AnnData
        annotation_counts = get_annotation_label_counts(adata)

        # If no data, show a simple message
        if not annotation_counts:
            return ui.tags.div("No annotations or data found.")

        # Build a list of annotation cards
        container = []
        for annotation_name, label_counts_dict in annotation_counts.items():
            # Sort labels by count (descending) and take top 5
            sorted_label_counts = sorted(
                label_counts_dict.items(),
                key=lambda x: x[1],
                reverse=True
            )[:10]

            # Build bullet points of "Label (count cells)"
            list_items = [
                ui.tags.li(f"{label} ({count} cells)")
                for label, count in sorted_label_counts
            ]

            # Wrap it in a card for this annotation
            annotation_card = ui.card(
                ui.h5(annotation_name),
                ui.tags.ul(*list_items),
                style="margin-bottom: 15px;"
            )
            container.append(annotation_card)

        # Return them as one TagList so each annotation is its own card
        return ui.TagList(*container)

    @reactive.effect
    def update_select_label_nn():
        adata = ad.AnnData(obs=shared['obs_data'].get())
        if input.nn_anno():
            selected_anno = input.nn_anno()
            labels = adata.obs[selected_anno].unique().tolist()
            ui.update_select("nn_anno_label", choices=labels)

    @reactive.effect
    def update_select_df_nn():
        df_names = shared['obsm_names'].get()
        ui.update_select("nn_spatial", choices=df_names)
            
    region_ui_initialized = reactive.Value(False)

    @reactive.Calc
    @render.text
    def print_obsm_names():
        obsm = shared['obsm_names'].get()
        if not obsm:
            return "Associated Tables: None"
        if obsm is not None:
            if len(obsm) > 1:
                obsm_str = ", ".join(obsm)
            else:
                obsm_str = obsm[0] if obsm else ""
            return "Associated Tables: " + obsm_str
        else:
            return "Empty"
        return


    # Initialize a flag to track dropdown creation
    subset_ui_initialized = reactive.Value(False)

    @reactive.effect
    def subset_reactivity():
        # Get input values
        adata = ad.AnnData(obs=shared['obs_data'].get())
        annotations = shared['obs_names'].get()
        btn = input.subset_select_check()

        # Access the flag value
        ui_initialized = subset_ui_initialized.get()

        if btn and not ui_initialized:
            # Insert annotation dropdown
            anno_dropdown = ui.input_select(
                "subset_anno_select", "Select Annotation to subset", choices=annotations
            )
            ui.insert_ui(
                ui.div({"id": "inserted-subset_anno_dropdown"}, anno_dropdown),
                selector="#main-subset_anno_dropdown",
                where="beforeEnd",
            )

            # Insert label dropdown
            label_dropdown = ui.input_selectize(
                "subset_label_select", "Select a Label", choices=[], selected=[], multiple=True
            )
            ui.insert_ui(
                ui.div({"id": "inserted-subset_label_dropdown"}, label_dropdown),
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
        trigger = label_update_trigger.get()  # Add trigger dependency
        adata = shared['adata_main'].get()
        selected_anno = input.subset_anno_select()

        if adata is not None and selected_anno:
            labels = adata.obs[selected_anno].unique().tolist()
            print(f"Updating labels for {selected_anno}: {labels}")
            ui.update_selectize("subset_label_select", choices=labels)

    @reactive.effect
    @reactive.event(input.go_subset, ignore_none=True)
    def subset_stratification():
        adata = shared['adata_main'].get()
        if adata is not None:
            annotation = input.subset_anno_select()
            labels = list(input.subset_label_select())

            # Perform the subsetting
            adata_subset = spac.data_utils.select_values(adata, annotation=annotation, values=labels)

            # Update the main data
            shared['adata_main'].set(adata_subset)

            # Increment the label update trigger to force update
            label_update_trigger.set(label_update_trigger.get() + 1)


    # Reactive variable to store subset history as a simple string
    subset_history = reactive.Value("")

    @reactive.effect
    @reactive.event(input.go_subset, ignore_none=True)
    def track_subset():
        """
        Append the current annotation and selected labels to the subset history.
        """
        annotation = input.subset_anno_select()
        labels = input.subset_label_select()

        # If both annotation and labels are valid, add them to history
        if annotation and labels:
            # Create a string like "Annotation:label1,label2"
            new_entry = f"{annotation}: {','.join(labels)}"

            # Update the subset history string by appending the new entry
            current_history = subset_history.get()
            if current_history:  # If there's already history, append with a separator
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
        """
        Store a master copy of the adata object when it is first loaded into adata_main.
        """
        adata = shared['adata_main'].get()
        if adata is not None and adata_master.get() is None:
            # Make a copy of the adata object and store it as the master copy
            adata_master.set(adata.copy())

    # Add a button to the UI for restoring the data
    ui.input_action_button("restore_data", "Restore Original Data", class_="btn-warning")

    @reactive.effect
    @reactive.event(input.restore_data, ignore_none=True)
    def restore_to_master():
        """
        Restore adata_main to the master copy stored in adata_master.
        """
        master_data = adata_master.get()
        if master_data is not None:
            # Set adata_main to a copy of the master data to ensure independence
            shared['adata_main'].set(master_data.copy())

            # Clear the subset history since the data is restored
            subset_history.set("")

