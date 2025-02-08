from shiny import App, Inputs, Outputs, Session, render, ui, reactive
from shinywidgets import output_widget, render_widget
import pickle
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path as path
import spac
import spac.visualization
import spac.spatial_analysis

app_ui = ui.page_fluid(

    ui.navset_card_tab(

        # 1. DATA INPUT PANEL -----------------------------------
        ui.nav_panel("Data Input",
            ui.card({"style": "width:100%;"},
                ui.column(12,
                    ui.row(
                        ui.column(6,
                            ui.div({"style": "font-weight: bold; font-size: 30px;"},
                                   ui.p("SPAC Interactive Dashboard")),
                            ui.input_file("input_file", "Choose a file to upload:", multiple=False),
                            ui.output_text("print_rows"),
                            ui.output_text("print_columns"),
                            ui.output_text("print_obs_names"),
                            ui.output_text("print_obsm_names"),
                            ui.output_text("print_layers_names"),
                            ui.output_text("print_uns_names")
                        ),
                        ui.column(6,
                            ui.input_checkbox("subset_select_check", "Subset Annotation", False),
                            ui.div(id="main-subset_anno_dropdown"),
                            ui.div(id="main-subset_label_dropdown"),
                            ui.input_action_button("go_subset", "Subset Data", class_="btn-success"),
                            ui.input_action_button("restore_data", "Restore Original Data", class_="btn-warning"),
                            ui.output_text("print_subset_history")
                        )
                    )
                )
            )
        ),

        # 2. ANNOTATIONS PANEL (Histogram of annotations) --------
        ui.nav_panel("Annotations",
            ui.card({"style": "width:100%;"},
                ui.column(12,
                    ui.row(
                        ui.column(2,
                            ui.input_select("h2_anno", "Select an Annotation", choices=[]),
                            ui.input_checkbox("h2_group_by_check", "Group By", value=False),
                            ui.div(id="main-h2_dropdown"),
                            ui.div(id="main-h2_check"),
                            ui.div(id="main-h2_together_drop"),
                            ui.input_action_button("go_h2", "Render Plot", class_="btn-success"),
                        ),
                        ui.column(10,
                            ui.div(
                            {"style": "padding-bottom: 100px;"},
                            ui.output_plot("spac_Histogram_2", width="100%", height="80vh")
                            )
                        )
                    )
                )
            )
        ),


        # 3. FEATURES PANEL (Histogram) --------------------------
        ui.nav_panel("Features",
            ui.card({"style": "width:100%;"},
                ui.column(12,
                    ui.row(
                        ui.column(2,
                            ui.input_select("h1_feat", "Select a Feature", choices=[]),
                            ui.input_select("h1_layer", "Select a Table", choices=[], selected=["Original"]),
                            ui.input_checkbox("h1_group_by_check", "Group By", value=False),
                            ui.input_checkbox("h1_log_x", "Log X-axis", value=False),
                            ui.input_checkbox("h1_log_y", "Log Y-axis", value=False),
                            ui.div(id="main-h1_dropdown"),
                            ui.div(id="main-h1_check"),
                            ui.div(id="main-h1_together_drop"),
                            ui.input_action_button("go_h1", "Render Plot", class_="btn-success")
                        ),
                        ui.column(10,
                            ui.div(
                            {"style": "padding-bottom: 100px;"},
                            ui.output_plot("spac_Histogram_1", width="100%", height="60vh")
                            )
                        )
                    )
                )
            )
        ),

        # 4. BOXPLOTS PANEL --------------------------------------
        ui.nav_panel("Boxplots",
            ui.card({"style": "width:100%;"},
                ui.column(12,
                    ui.row(
                        ui.column(6,
                            ui.card({"style": "width:100%;"},
                                ui.column(12,
                                    ui.input_select("bp1_anno", "Select an Annotation", choices=[]),
                                    ui.input_select("bp1_layer", "Select a Table", choices=[], selected="Original"),
                                    ui.input_selectize("bp1_features", "Select Features", multiple=True, choices=[], selected=[]),
                                    ui.input_checkbox("bp1_outlier_check", "Add Outliers", False),
                                    ui.input_action_button("go_bp1", "Render Plot", class_="btn-success"),
                                    ui.div(
                                    {"style": "padding-bottom: 100px;"},
                                    ui.output_plot("spac_Boxplot_1", width="100%", height="60vh")
                                    )
                                )
                            )
                        ),
                        ui.column(6,
                            ui.card({"style": "width:100%;"},
                                ui.column(12,
                                    ui.input_select("bp2_anno", "Select an Annotation", choices=[]),
                                    ui.input_select("bp2_layer", "Select a Table", choices=[], selected="Original"),
                                    ui.input_selectize("bp2_features", "Select Features", multiple=True, choices=[], selected=[]),
                                    ui.input_checkbox("bp2_outlier_check", "Add Outliers", False),
                                    ui.input_action_button("go_bp2", "Render Plot", class_="btn-success"),
                                    ui.div(
                                    {"style": "padding-bottom: 100px;"},
                                    ui.output_plot("spac_Boxplot_2", width="100%", height="60vh")
                                    )
                                )
                            )
                        )
                    )
                )
            )
        ),


        # 5. FEAT. VS ANNO. (Heatmap) ----------------------------
        ui.nav_panel("Feat. Vs Anno.",
            ui.card({"style": "width:100%;"},
                ui.column(12,
                    ui.row(
                        ui.column(2,
                            ui.input_select("hm1_anno", "Select an Annotation", choices=[]),
                            ui.input_select("hm1_layer", "Select a Table", choices=[]),
                            ui.input_checkbox("dendogram", "Include Dendrogram", False),
                            ui.div(id="main-hm1_check"),
                            ui.div(id="main-hm2_check"),
                            ui.div(id="main-min_num"),
                            ui.div(id="main-max_num"),
                            ui.input_action_button("go_hm1", "Render Plot", class_="btn-success"),
                            ui.div(
                                    {"style": "padding-top: 20px;"},
                                    ui.output_ui("download_button_ui")
                                )
                        ),
                        ui.column(10,
                            ui.div(
                                    {"style": "padding-bottom: 100px;"},
                            ui.output_plot("spac_Heatmap", width="100%", height="60vh")
                            )
                        )
                    )
                )
            )
        ),

        # 6. ANNO. VS ANNO. (Sankey, Relational Heatmap) ---------
        ui.nav_panel("Anno. Vs Anno.",
            ui.card({"style": "width:100%;"},
                ui.column(12,
                    ui.row(
                        ui.column(2,
                            ui.input_select("sk1_anno1", "Select Source Annotation", choices=[]),
                            ui.input_select("sk1_anno2", "Select Target Annotation", choices=[]),
                            ui.input_action_button("go_sk1", "Render Plot", class_="btn-success")
                        ),
                        ui.column(10,
                            ui.div(
                                output_widget("spac_Sankey"),
                                style="width:100%; height:80vh;"
                            )
                        )
                    )
                )
            ),
            ui.card({"style": "width:100%;"},
                ui.column(12,
                    ui.row(
                        ui.column(2,
                            ui.input_select("rhm_anno1", "Select Source Annotation", choices=[], selected=[]),
                            ui.input_select("rhm_anno2", "Select Target Annotation", choices=[], selected=[]),
                            ui.input_action_button("go_rhm1", "Render Plot", class_="btn-success"),
                            ui.div(
                                    {"style": "padding-top: 20px;"},
                                    ui.output_ui("download_button_ui_1")
                                )
                        ),
                        ui.column(10,
                            ui.div(
                                output_widget("spac_Relational"),
                                style="width:100%; height:80vh;"
                            )
                        )
                    )
                )
            )
        ),

        # 7. SPATIAL PANEL ---------------------------------------
        ui.nav_panel("Spatial",
            ui.card({"style": "width:100%;"},
                ui.column(12,
                    ui.row(
                        ui.column(2,
                            ui.input_select("spatial_anno", "Select an Annotation", choices=[]),
                            ui.input_slider("spatial_slider", "Point Size", min=2, max=10, value=3),
                            ui.input_checkbox("slide_select_check", "Stratify by Slide", False),
                            ui.div(id="main-slide_dropdown"),
                            ui.div(id="main-label_dropdown"),
                            ui.input_checkbox("region_select_check", "Stratify by Region", False),
                            ui.div(id="main-region_dropdown"),
                            ui.div(id="main-region_label_select_dropdown"),
                            ui.input_action_button("go_sp1", "Render Plot", class_="btn-success")
                        ),
                        ui.column(10,
                            ui.div(
                                {"style": "padding-bottom: 20px;"},
                                output_widget("spac_Spatial"),
                                style="width:100%; height:80vh;"
                            )
                        )
                    )
                )
            )
        ),

        # 8. UMAP PANEL ------------------------------------------
        ui.nav_panel("UMAP",
            ui.card({"style": "width:100%;"},
                ui.column(12,
                    ui.row(
                        ui.column(6,
                            ui.input_radio_buttons("umap_rb", "Choose one:", ["Annotation", "Feature"]),
                            ui.input_select("plottype", "Select a plot type", choices=["umap", "pca", "tsne"]),
                            ui.div(id="main-ump_rb_dropdown_anno"),
                            ui.div(id="main-ump_rb_dropdown_feat"),
                            ui.div(id="main-ump_table_dropdown_feat"),
                            ui.input_slider("umap_slider_1", "Point Size", min=.5, max=10, value=3),
                            ui.input_action_button("go_umap1", "Render Plot", class_="btn-success"),
                            ui.output_plot("spac_UMAP", width="100%", height="80vh")
                        ),
                        ui.column(6,
                            ui.input_radio_buttons("umap_rb2", "Choose one:", ["Annotation", "Feature"]),
                            ui.input_select("plottype2", "Select a plot type", choices=["umap", "pca", "tsne"]),
                            ui.div(id="main-ump_rb_dropdown_anno2"),
                            ui.div(id="main-ump_rb_dropdown_feat2"),
                            ui.div(id="main-ump_table_dropdown_feat2"),
                            ui.input_slider("umap_slider_2", "Point Size", min=.5, max=10, value=3),
                            ui.input_action_button("go_umap2", "Render Plot", class_="btn-success"),
                            ui.output_plot("spac_UMAP2", width="100%", height="80vh")
                        )
                    )
                )
            )
        ),

        # 9. SCATTERPLOT PANEL ------------------------------------
        ui.nav_panel("Scatterplot",
            ui.card({"style": "width:100%;"},
                ui.column(12,
                    ui.row(
                        ui.column(2,
                            ui.input_select("scatter_layer", "Select a Table", choices=[], selected="Original"),
                            ui.input_select("scatter_x", "Select X Axis", choices=[]),
                            ui.input_select("scatter_y", "Select Y Axis", choices=[]),
                            ui.input_checkbox("scatter_color_check", "Color by Feature", value=False),
                            ui.div(id="main-scatter_dropdown"),
                            ui.input_action_button("go_scatter", "Render Plot", class_="btn-success")
                        ),
                        ui.column(10,
                            ui.output_plot("spac_Scatter", width="100%", height="80vh")
                        )
                    )
                )
            )
        )
    )
)


def server(input, output, session):



    # Define a reactive variable to track if data is loaded
    data_loaded = reactive.Value(False)

    @reactive.Effect
    def adata_filter():
        print("Calling Data")
        file_info = input.input_file()
        if not file_info:
            data_loaded.set(False)  # Set to False if no file is uploaded
            return
        else:
            file_path = file_info[0]['datapath']
            with open(file_path, 'rb') as file:
                if file_path.endswith('.pickle'):
                    adata_main.set(pickle.load(file))
                elif file_path.endswith('.h5ad'):
                    adata_main.set(ad.read_h5ad(file_path))
                else:
                    adata_main.set(ad.read(file_path))
            data_loaded.set(True)  # Set to True if a file is successfully uploaded





    # Create a reactive variable for the main data
    adata_main = reactive.Value(None)

    # Create reactive variables for parts of the anndata object
    slide_annotation = reactive.Value(None)
    X_data = reactive.Value(None)
    obs_data = reactive.Value(None) #AKA Annotations
    obsm_data = reactive.Value(None)
    layers_data = reactive.Value(None)
    var_data = reactive.Value(None) #AKA Features
    uns_data = reactive.Value(None)
    shape_data = reactive.Value(None)
    obs_names= reactive.Value(None)
    obsm_names = reactive.Value(None)
    layers_names = reactive.Value(None)
    var_names = reactive.Value(None)
    uns_names = reactive.Value(None)
    df_heatmap = reactive.Value(None)
    df_relational = reactive.Value(None)

    @reactive.Effect
    def update_parts():
        print("Updating Parts")
        adata = adata_main.get()
        if adata is not None:

            if hasattr(adata, 'X'):
                X_data.set(adata.X)
            else:
                X_data.set(None)

            if hasattr(adata, 'obs'):
                obs_data.set(adata.obs)
            else:
                obs_data.set(None)

            if hasattr(adata, 'obsm'):
                obsm_data.set(adata.obsm)
            else:
                obsm_data.set(None)

            if hasattr(adata, 'layers'):
                layers_data.set(adata.layers)
            else:
                layers_data.set(None)

            if hasattr(adata, 'var'):
                var_data.set(adata.var)
            else:
                var_data.set(None)

            if hasattr(adata, 'uns'):
                uns_data.set(adata.uns)
            else:
                uns_data.set(None)

            shape_data.set(adata.shape)

            if hasattr(adata, 'obs'):
                obs_names.set(list(adata.obs.keys()))
            else:
                obs_names.set(None)

            if hasattr(adata, 'obsm'):
                obsm_names.set(list(adata.obsm.keys()))
            else:
                obsm_names.set(None)

            if hasattr(adata, 'layers'):
                layers_names.set(list(adata.layers.keys()))
            else:
                layers_names.set(None)

            if hasattr(adata, 'var'):
                var_names.set(list(adata.var.index.tolist()))
            else:
                var_names.set(None)

            if hasattr(adata, 'uns'):
                uns_names.set(list(adata.uns.keys()))
            else:
                uns_names.set(None)
        else:
            obs_data.set(None)
            obsm_data.set(None)
            layers_data.set(None)
            var_data.set(None)
            uns_data.set(None)
            shape_data.set(None)
            obs_names.set(None)
            obsm_names.set(None)
            layers_names.set(None)
            var_names.set(None)
            uns_names.set(None)


    @reactive.Calc
    @render.text
    def print_obs_names():
        obs = obs_names.get()
        if obs is not None:
            if len(obs) > 1:
                obs_str = ", ".join(obs)
            else:
                obs_str = obs[0] if obs else ""
            return "Annotations: " + obs_str
        return

    @reactive.Calc
    @render.text
    def print_obsm_names():
        obsm = obsm_names.get()
        if obsm is not None:
            if len(obsm) > 1:
                obsm_str = ", ".join(obsm)
            else:
                obsm_str = obsm[0] if obsm else ""
            return "Associated Tables: " + obsm_str
        return

    @reactive.Calc
    @render.text
    def print_layers_names():
        layers = layers_names.get()
        if layers is not None:
            if len(layers) > 1:
                layers_str = ", ".join(layers)
            elif len(layers) > 1:
                layers_str = layers[0] if layers else ""
            return "Tables: " + layers_str
        return

    @reactive.Calc
    @render.text
    def print_uns_names():
        uns = uns_names.get()
        if uns is not None:
            if len(uns) > 1:
                uns_str = ", ".join(uns)
            else:
                uns_str = uns[0] if uns else ""
            return "Unstructured Data: " + uns_str
        return

    @reactive.Calc
    @render.text
    def print_rows():
        shape = shape_data.get()
        if shape is not None:
            return "# of Cells: " + str(shape[0])
        return

    @reactive.Calc
    @render.text
    def print_columns():
        shape = shape_data.get()
        if shape is not None:
            return "# of Features: " + str(shape[1])
        return



    @reactive.Effect
    def update_select_input_feat():
        choices = var_names.get()
        ui.update_select("h1_feat", choices=choices)
        ui.update_select("umap_feat", choices=choices)
        ui.update_select("bp1_features", choices=choices)
        ui.update_select("bp2_features", choices=choices)


    @reactive.Effect
    def update_select_input_anno():
        choices = obs_names.get()
        ui.update_select("bp1_anno", choices=choices)
        ui.update_select("bp2_anno", choices=choices)
        ui.update_select("h2_anno", choices=choices)
        ui.update_select("hm1_anno", choices=choices)
        ui.update_select("sk1_anno1", choices=choices)
        ui.update_select("sk1_anno2", choices=choices)
        ui.update_select("rhm_anno1", choices=choices)
        ui.update_select("rhm_anno2", choices=choices)
        ui.update_select("spatial_anno", choices=choices)

        return

    @reactive.Effect
    def update_select_input_layer():
        if layers_names.get() is not None:
            new_choices = layers_names.get() + ["Original"]
            ui.update_select("h1_layer", choices=new_choices)
            ui.update_select("bp1_layer", choices=new_choices)
            ui.update_select("bp2_layer", choices=new_choices)
            ui.update_select("hm1_layer", choices=new_choices)
            ui.update_select("scatter_layer", choices=new_choices)
        return
    @reactive.Effect
    def update_select_input_anno_bp():
        if obs_names.get() is not None:
            new_choices = obs_names.get() + ["No Annotation"]
            ui.update_select("bp1_anno", choices=new_choices)
            ui.update_select("bp2_anno", choices=new_choices)

    @reactive.Effect
    def update_select_input_layer_scatter():
        choices = get_scatterplot_names()
        ui.update_select("scatter_x", choices=choices)
        ui.update_select("scatter_y", choices=choices)
        return


    @reactive.Effect
    def update_boxplot_selectize():
        selected_names=var_names.get()
        if selected_names is not None:
            ui.update_selectize("bp1_features", selected=selected_names[:2])
            ui.update_selectize("bp2_features", selected=selected_names[:2])
            return
    @reactive.Effect
    def update_relational_select():
        selected_names=obs_names.get()
        if selected_names is not None and len(selected_names) > 1:
            ui.update_selectize("rhm_anno1", selected=selected_names[0])
            ui.update_selectize("rhm_anno2", selected=selected_names[1])
        return



    # Initialize a flag to track dropdown creation
    subset_ui_initialized = reactive.Value(False)

    @reactive.effect
    def subset_reactivity():
        # Get input values
        adata = ad.AnnData(obs=obs_data.get())
        annotations = obs_names.get()
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
        adata = adata_main.get()
        selected_anno = input.subset_anno_select()

        if adata is not None and selected_anno:
            labels = adata.obs[selected_anno].unique().tolist()
            print(f"Updating labels for {selected_anno}: {labels}")
            ui.update_selectize("subset_label_select", choices=labels)

    @reactive.effect
    @reactive.event(input.go_subset, ignore_none=True)
    def subset_stratification():
        adata = adata_main.get()
        if adata is not None:
            annotation = input.subset_anno_select()
            labels = input.subset_label_select()

            # Perform the subsetting
            adata_subset = adata[adata.obs[annotation].isin(labels)].copy()

            # Update the main data
            adata_main.set(adata_subset)

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
        adata = adata_main.get()
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
            adata_main.set(master_data.copy())

            # Clear the subset history since the data is restored
            subset_history.set("")





    @output
    @render.plot
    @reactive.event(input.go_h1, ignore_none=True)
    def spac_Histogram_1():
        adata = ad.AnnData(X=X_data.get(), obs=pd.DataFrame(obs_data.get()), var=pd.DataFrame(var_data.get()), layers=layers_data.get(), dtype=X_data.get().dtype)
        btn_log_x = input.h1_log_x()
        btn_log_y = input.h1_log_y()
        if adata is not None:
            if input.h1_group_by_check() is not True:
                if input.h1_layer() != "Original":
                    fig1 = spac.visualization.histogram(adata, feature=input.h1_feat(), layer=input.h1_layer(), log_scale=(btn_log_x, btn_log_y))
                else:
                    fig1 = spac.visualization.histogram(adata, feature=input.h1_feat(), log_scale=(btn_log_x, btn_log_y))

            if input.h1_group_by_check() is not False:
                if input.h1_layer() != "Original":
                    if input.h1_together_check() is  not False:
                        fig1 = spac.visualization.histogram(adata, feature=input.h1_feat(), layer=input.h1_layer(), group_by=input.h1_anno(), together=input.h1_together_check(), log_scale=(btn_log_x, btn_log_y), multiple=input.h1_together_drop())
                    else:
                        fig1 = spac.visualization.histogram(adata, feature=input.h1_feat(), layer=input.h1_layer(), group_by=input.h1_anno(), together=input.h1_together_check(), log_scale=(btn_log_x, btn_log_y))
                else:
                    if input.h1_together_check() is  not False:
                        fig1 = spac.visualization.histogram(adata, feature=input.h1_feat(), group_by=input.h1_anno(), together=input.h1_together_check(), log_scale=(btn_log_x, btn_log_y), multiple=input.h1_together_drop())
                    else:
                        fig1 = spac.visualization.histogram(adata, feature=input.h1_feat(), group_by=input.h1_anno(), together=input.h1_together_check(), log_scale=(btn_log_x, btn_log_y))
            fig = fig1[0]
            axes = fig1[1] if isinstance(fig1[1], list) else [fig1[1]]
            for ax in axes:
                ax.tick_params(axis='x', rotation=90, labelsize=10)
            return fig1
        return None

    histogram_ui_initialized = reactive.Value(False)

    @reactive.effect
    def histogram_reactivity():
        btn = input.h1_group_by_check()
        ui_initialized = histogram_ui_initialized.get()

        if btn and not ui_initialized:
            dropdown = ui.input_select("h1_anno", "Select an Annotation", choices=obs_names.get())
            ui.insert_ui(
                ui.div({"id": "inserted-dropdown"}, dropdown),
                selector="#main-h1_dropdown",
                where="beforeEnd",
            )

            together_check = ui.input_checkbox("h1_together_check", "Plot Together", value=False)
            ui.insert_ui(
                ui.div({"id": "inserted-check"}, together_check),
                selector="#main-h1_check",
                where="beforeEnd",
            )

            histogram_ui_initialized.set(True)

        elif not btn and ui_initialized:
            ui.remove_ui("#inserted-dropdown")
            ui.remove_ui("#inserted-check")
            ui.remove_ui("#inserted-dropdown_together")
            histogram_ui_initialized.set(False)

    @reactive.effect
    @reactive.event(input.h1_together_check)
    def update_stack_type_dropdown():
        if input.h1_together_check():
            dropdown_together = ui.input_select("h1_together_drop", "Select Stack Type", 
                                                choices=['stack', 'layer', 'dodge', 'fill'], 
                                                selected='stack')
            ui.insert_ui(
                ui.div({"id": "inserted-dropdown_together"}, dropdown_together),
                selector="#main-h1_together_drop",
                where="beforeEnd",)      
        else:
            ui.remove_ui("#inserted-dropdown_together")

    @output
    @render.plot
    @reactive.event(input.go_bp1, ignore_none=True)
    def spac_Boxplot_1():
        adata = ad.AnnData(X=X_data.get(), obs=pd.DataFrame(obs_data.get()), var=pd.DataFrame(var_data.get()), layers=layers_data.get(), dtype=X_data.get().dtype)
        if adata is not None and adata.var is not None:
            if input.bp1_layer() != "Original" and input.bp1_anno() != "No Annotation":
                fig, ax, df = spac.visualization.boxplot(adata, annotation=input.bp1_anno(), layer=input.bp1_layer(), features=list(input.bp1_features()),showfliers=input.bp1_outlier_check())
                return ax.legend(loc='upper right', bbox_to_anchor=(1, 1))
            if input.bp1_layer() == "Original" and input.bp1_anno() != "No Annotation":
                fig, ax, df = spac.visualization.boxplot(adata, annotation=input.bp1_anno(), features=list(input.bp1_features()),showfliers=input.bp1_outlier_check())
                return ax.legend(loc='upper right', bbox_to_anchor=(1, 1))
            if input.bp1_layer() != "Original" and input.bp1_anno() == "No Annotation":
                fig, ax, df = spac.visualization.boxplot(adata, layer=input.bp1_layer(), features=list(input.bp1_features()),showfliers=input.bp1_outlier_check())
                return ax.legend(loc='upper right', bbox_to_anchor=(1, 1))
            if input.bp1_layer() == "Original" and input.bp1_anno() == "No Annotation":
                fig, ax, df = spac.visualization.boxplot(adata, features=list(input.bp1_features()),showfliers=input.bp1_outlier_check())
                return ax.legend(loc='upper right', bbox_to_anchor=(1, 1))
        return None

    @output
    @render.plot
    @reactive.event(input.go_bp2, ignore_none=True)
    def spac_Boxplot_2():
        adata = ad.AnnData(X=X_data.get(), obs=pd.DataFrame(obs_data.get()), var=pd.DataFrame(var_data.get()), layers=layers_data.get(), dtype=X_data.get().dtype)
        if adata is not None and adata.var is not None:
            if input.bp2_layer() != "Original" and input.bp2_anno() != "No Annotation":
                fig, ax, df = spac.visualization.boxplot(adata, annotation=input.bp2_anno(), layer=input.bp2_layer(), features=list(input.bp2_features()),showfliers=input.bp2_outlier_check())
                return ax.legend(loc='upper right', bbox_to_anchor=(1, 1))
            if input.bp2_layer() == "Original" and input.bp2_anno() != "No Annotation":
                fig, ax, df = spac.visualization.boxplot(adata, annotation=input.bp2_anno(), features=list(input.bp2_features()),showfliers=input.bp2_outlier_check())
                return ax.legend(loc='upper right', bbox_to_anchor=(1, 1))
            if input.bp2_layer() != "Original" and input.bp2_anno() == "No Annotation":
                fig, ax, df = spac.visualization.boxplot(adata, layer=input.bp2_layer(), features=list(input.bp2_features()),showfliers=input.bp2_outlier_check())
                return ax.legend(loc='upper right', bbox_to_anchor=(1, 1))
            if input.bp2_layer() == "Original" and input.bp2_anno() == "No Annotation":
                fig, ax, df = spac.visualization.boxplot(adata, features=list(input.bp2_features()),showfliers=input.bp2_outlier_check())
                return ax.legend(loc='upper right', bbox_to_anchor=(1, 1))
        return None

    @output
    @render.plot
    @reactive.event(input.go_h2, ignore_none=True)
    def spac_Histogram_2():
        adata = adata_main.get()
        if adata is None:
            return None

        # 1) If "Group By" is UNCHECKED, show a simple annotation histogram
        if not input.h2_group_by_check():
            fig = spac.visualization.histogram(
                adata,
                annotation=input.h2_anno()
            )
            fig[1].tick_params(axis='x', rotation=90, labelsize=10)
            return fig

        # 2) If "Group By" is CHECKED, we must always supply a valid multiple parameter
        else:
            # If user also checked "Plot Together", use their selected stack type
            if input.h2_together_check():
                multiple_param = input.h2_together_drop()  # e.g. 'stack', 'dodge', etc.
                together_flag = True
            else:
                # If grouping by but not "plot together", pick a default layout
                multiple_param = "layer"  # or 'dodge' or any valid string
                together_flag = False

            fig = spac.visualization.histogram(
                adata,
                annotation=input.h2_anno(),
                group_by=input.h2_anno_1(),
                together=together_flag,
                multiple=multiple_param
            )
            return fig

        return None

    histogram2_ui_initialized = reactive.Value(False)

    @reactive.effect
    def histogram_reactivity_2():
        btn = input.h2_group_by_check()
        ui_initialized = histogram2_ui_initialized.get()

        if btn and not ui_initialized:
            dropdown = ui.input_select("h2_anno_1", "Select an Annotation", choices=obs_names.get())
            ui.insert_ui(
                ui.div({"id": "inserted-dropdown-1"}, dropdown),
                selector="#main-h2_dropdown",
                where="beforeEnd",
            )

            together_check = ui.input_checkbox("h2_together_check", "Plot Together", value=False)
            ui.insert_ui(
                ui.div({"id": "inserted-check-1"}, together_check),
                selector="#main-h2_check",
                where="beforeEnd",
            )
            histogram2_ui_initialized.set(True)

        elif not btn and ui_initialized:
            ui.remove_ui("#inserted-dropdown-1")
            ui.remove_ui("#inserted-check-1")
            ui.remove_ui("#inserted-dropdown_together-1")
            histogram2_ui_initialized.set(False)
    
    @reactive.effect
    @reactive.event(input.h2_together_check)
    def update_stack_type_dropdown():
        if input.h2_together_check():
            dropdown_together = ui.input_select("h2_together_drop", "Select Stack Type", 
                                                choices=['stack', 'layer', 'dodge', 'fill'], 
                                                selected='stack')
            ui.insert_ui(
                ui.div({"id": "inserted-dropdown_together-1"}, dropdown_together),
                selector="#main-h2_together_drop",
                where="beforeEnd",)      
        else:
            ui.remove_ui("#inserted-dropdown_together-1")

    @output
    @render.plot
    @reactive.event(input.go_hm1, ignore_none=True)
    def spac_Heatmap():
        adata = ad.AnnData(X=X_data.get(), obs=pd.DataFrame(obs_data.get()), var=pd.DataFrame(var_data.get()), layers=layers_data.get(), dtype=X_data.get().dtype)
        if adata is not None:
            vmin = input.min_select()
            vmax = input.max_select()    
            if input.dendogram() is not True:
                if input.hm1_layer() != "Original":
                    df, fig, ax = spac.visualization.hierarchical_heatmap(adata, annotation=input.hm1_anno(), layer=input.hm1_layer(), z_score=None, vmin=vmin, vmax=vmax)
                else:
                    df, fig, ax = spac.visualization.hierarchical_heatmap(adata, annotation=input.hm1_anno(), layer=None, z_score=None, vmin=vmin, vmax=vmax)
            elif input.dendogram() is not False:
                cluster_annotations = input.h2_anno_dendro()
                cluster_features = input.h2_feat_dendro()
                if input.hm1_layer() != "Original":
                    df, fig, ax = spac.visualization.hierarchical_heatmap(adata, annotation=input.hm1_anno(), layer=input.hm1_layer(), z_score=None, cluster_annotations=cluster_annotations, cluster_feature=cluster_features, vmin=vmin, vmax=vmax)
                else:
                    df, fig, ax = spac.visualization.hierarchical_heatmap(adata, annotation=input.hm1_anno(), layer=None, z_score=None, cluster_annotations=cluster_annotations, cluster_feature=cluster_features, vmin=vmin, vmax=vmax)
            df_heatmap.set(df)
            return fig

        return None

    heatmap_ui_initialized = reactive.Value(False)

    @reactive.effect
    def heatmap_reactivity():
        btn = input.dendogram()
        ui_initialized = heatmap_ui_initialized.get()

        if btn and not ui_initialized:
            annotation_check = ui.input_checkbox("h2_anno_dendro", "Annotation Cluster", value=False)
            ui.insert_ui(
                ui.div({"id": "inserted-check"}, annotation_check),
                selector="#main-hm1_check",
                where="beforeEnd",
            )

            feat_check = ui.input_checkbox("h2_feat_dendro", "Feature Cluster", value=False)
            ui.insert_ui(
                ui.div({"id": "inserted-check1"}, feat_check),
                selector="#main-hm2_check",
                where="beforeEnd",
            )
            heatmap_ui_initialized.set(True)

        elif not btn and ui_initialized:
            ui.remove_ui("#inserted-check")
            ui.remove_ui("#inserted-check1")
            heatmap_ui_initialized.set(False)

    @session.download(filename="heatmap_data.csv")
    def download_df():
        df = df_heatmap.get()
        if df is not None:
            csv_string = df.to_csv(index=False)
            csv_bytes = csv_string.encode("utf-8")
            return csv_bytes, "text/csv"
        return None

    @render.ui
    @reactive.event(input.go_hm1, ignore_none=True)
    def download_button_ui():
        if df_heatmap.get() is not None:
            return ui.download_button("download_df", "Download Data", class_="btn-warning")
        return None

    @reactive.effect
    @reactive.event(input.hm1_layer)
    def update_min_max():
        adata = ad.AnnData(X=X_data.get(), obs=pd.DataFrame(obs_data.get()), var=pd.DataFrame(var_data.get()), layers=layers_data.get())
        if input.hm1_layer() == "Original":
            layer_data = adata.X
        else:
            layer_data = adata.layers[input.hm1_layer()]
        mask = adata.obs[input.hm1_anno()].notna()
        layer_data = layer_data[mask]
        min_val = round(float(np.min(layer_data)), 2)
        max_val = round(float(np.max(layer_data)), 2)

        ui.remove_ui("#inserted-min_num")
        ui.remove_ui("#inserted-max_num")

        min_num = ui.input_numeric("min_select", "Minimum", min_val, min=min_val, max=max_val)
        ui.insert_ui(
            ui.div({"id": "inserted-min_num"}, min_num),
            selector="#main-min_num",
            where="beforeEnd",
        )
        
        max_num = ui.input_numeric("max_select", "Maximum", max_val, min=min_val, max=max_val)
        ui.insert_ui(
            ui.div({"id": "inserted-max_num"}, max_num),
            selector="#main-max_num",
            where="beforeEnd",
        )

    @output
    @render_widget
    @reactive.event(input.go_sk1, ignore_none=True)
    def spac_Sankey():
        adata = ad.AnnData(X=X_data.get(), obs=pd.DataFrame(obs_data.get()), layers=layers_data.get(), dtype=X_data.get().dtype)
        if adata is not None:
            fig = spac.visualization.sankey_plot(adata, source_annotation=input.sk1_anno1(), target_annotation=input.sk1_anno2())
            return fig
        return None

    @output
    @render_widget
    @reactive.event(input.go_rhm1, ignore_none=True)
    def spac_Relational():
        adata = ad.AnnData(X=X_data.get(), obs=pd.DataFrame(obs_data.get()))
        if adata is not None:
            result = spac.visualization.relational_heatmap(adata, source_annotation=input.rhm_anno1(), target_annotation=input.rhm_anno2())
            df_relational.set(result['data'])
            return result['figure']
        return None

    @session.download(filename="relational_data.csv")
    def download_df_1():
        df = df_relational.get()
        if df is not None:
            csv_string = df.to_csv(index=False)
            csv_bytes = csv_string.encode("utf-8")
            return csv_bytes, "text/csv"
        return None

    @render.ui
    @reactive.event(input.go_rhm1, ignore_none=True)
    def download_button_ui_1():
        if df_relational.get() is not None:
            return ui.download_button("download_df_1", "Download Data", class_="btn-warning")
        return None

    @output
    @render.plot
    @reactive.event(input.go_umap1, ignore_none=True)
    def spac_UMAP():
        adata = ad.AnnData(X=X_data.get(), var=pd.DataFrame(var_data.get()), obsm=obsm_data.get(), obs=obs_data.get(), dtype=X_data.get().dtype, layers=layers_data.get())
        point_size=input.umap_slider_1()
        if adata is not None:
            if input.umap_rb() == "Feature":
                if input.umap_layer() == "Original":
                    layer = None
                else:
                    layer = input.umap_layer()
                out = spac.visualization.dimensionality_reduction_plot(adata, method=input.plottype(), feature=input.umap_rb_feat(), layer=layer, point_size=point_size)
                return out
            elif input.umap_rb() == "Annotation":
                out1 = spac.visualization.dimensionality_reduction_plot(adata, method=input.plottype(), annotation=input.umap_rb_anno(), point_size=point_size)
                return out1
        return None

    # Track the UI state
    umap_annotation_initialized = reactive.Value(False)
    umap_feature_initialized = reactive.Value(False)

    @reactive.effect
    def umap_reactivity():
        flipper = data_loaded.get()
        if flipper is not False:
            btn = input.umap_rb()

            if btn == "Annotation":
                if not umap_annotation_initialized.get():
                    # Create the Annotation dropdown
                    dropdown = ui.input_select(
                        "umap_rb_anno", "Select an Annotation", choices=obs_names.get()
                    )
                    ui.insert_ui(
                        ui.div({"id": "inserted-rbdropdown_anno"}, dropdown),
                        selector="#main-ump_rb_dropdown_anno",
                        where="beforeEnd",
                    )
                    # Update the state
                    umap_annotation_initialized.set(True)
                # Remove the Feature dropdown and table
                if umap_feature_initialized.get():
                    ui.remove_ui("#inserted-rbdropdown_feat")
                    ui.remove_ui("#inserted-umap_table")
                    umap_feature_initialized.set(False)

            elif btn == "Feature":
                if not umap_feature_initialized.get():
                    # Create the Feature dropdown
                    dropdown1 = ui.input_select(
                        "umap_rb_feat", "Select a Feature", choices=var_names.get()
                    )
                    ui.insert_ui(
                        ui.div({"id": "inserted-rbdropdown_feat"}, dropdown1),
                        selector="#main-ump_rb_dropdown_feat",
                        where="beforeEnd",
                    )

                    # Create the Table dropdown
                    new_choices = layers_names.get() + ["Original"]
                    table_umap = ui.input_select(
                        "umap_layer", "Select a Table", choices=new_choices, selected=["Original"]
                    )
                    ui.insert_ui(
                        ui.div({"id": "inserted-umap_table"}, table_umap),
                        selector="#main-ump_table_dropdown_feat",
                        where="beforeEnd",
                    )
                    # Update the state
                    umap_feature_initialized.set(True)
                # Remove the Annotation dropdown
                if umap_annotation_initialized.get():
                    ui.remove_ui("#inserted-rbdropdown_anno")
                    umap_annotation_initialized.set(False)

            elif btn == "None":
                # Remove all dropdowns and reset states
                if umap_annotation_initialized.get():
                    ui.remove_ui("#inserted-rbdropdown_anno")
                    umap_annotation_initialized.set(False)
                if umap_feature_initialized.get():
                    ui.remove_ui("#inserted-rbdropdown_feat")
                    ui.remove_ui("#inserted-umap_table")
                    umap_feature_initialized.set(False)



    @output
    @render.plot
    @reactive.event(input.go_umap2, ignore_none=True)
    def spac_UMAP2():
        adata = ad.AnnData(X=X_data.get(), var=pd.DataFrame(var_data.get()), obsm=obsm_data.get(), obs=obs_data.get(), dtype=X_data.get().dtype, layers=layers_data.get())
        point_size_2=input.umap_slider_2()
        if adata is not None:
            if input.umap_rb2() == "Feature":
                if input.umap_layer2() == "Original":
                    layer2 = None
                else:
                    layer2 = input.umap_layer2()
                out = spac.visualization.dimensionality_reduction_plot(adata, method=input.plottype2(), feature=input.umap_rb_feat2(), layer=layer2, point_size=point_size_2)
                return out
            elif input.umap_rb2() == "Annotation":
                out1 = spac.visualization.dimensionality_reduction_plot(adata, method=input.plottype2(), annotation=input.umap_rb_anno2(), point_size=point_size_2)
                return out1
        return None

    # Track the UI state
    umap2_annotation_initialized = reactive.Value(False)
    umap2_feature_initialized = reactive.Value(False)

    @reactive.effect
    def umap_reactivity2():
        flipper = data_loaded.get()
        if flipper is not False:
            btn = input.umap_rb2()

            if btn == "Annotation":
                if not umap2_annotation_initialized.get():
                    dropdown = ui.input_select(
                        "umap_rb_anno2", "Select an Annotation", choices=obs_names.get()
                    )
                    ui.insert_ui(
                        ui.div({"id": "inserted-rbdropdown_anno2"}, dropdown),
                        selector="#main-ump_rb_dropdown_anno2",
                        where="beforeEnd",
                    )
                    umap2_annotation_initialized.set(True)
                if umap2_feature_initialized.get():
                    ui.remove_ui("#inserted-rbdropdown_feat2")
                    ui.remove_ui("#inserted-umap_table2")
                    umap2_feature_initialized.set(False)

            elif btn == "Feature":
                if not umap2_feature_initialized.get():
                    dropdown1 = ui.input_select(
                        "umap_rb_feat2", "Select a Feature", choices=var_names.get()
                    )
                    ui.insert_ui(
                        ui.div({"id": "inserted-rbdropdown_feat2"}, dropdown1),
                        selector="#main-ump_rb_dropdown_feat2",
                        where="beforeEnd",
                    )

                    new_choices = layers_names.get() + ["Original"]
                    table_umap_1 = ui.input_select(
                        "umap_layer2", "Select a Table", choices=new_choices, selected=["Original"]
                    )
                    ui.insert_ui(
                        ui.div({"id": "inserted-umap_table2"}, table_umap_1),
                        selector="#main-ump_table_dropdown_feat2",
                        where="beforeEnd",
                    )
                    umap2_feature_initialized.set(True)
                if umap2_annotation_initialized.get():
                    ui.remove_ui("#inserted-rbdropdown_anno2")
                    umap2_annotation_initialized.set(False)

            elif btn == "None":
                if umap2_annotation_initialized.get():
                    ui.remove_ui("#inserted-rbdropdown_anno2")
                    umap2_annotation_initialized.set(False)
                if umap2_feature_initialized.get():
                    ui.remove_ui("#inserted-rbdropdown_feat2")
                    ui.remove_ui("#inserted-umap_table2")
                    umap2_feature_initialized.set(False)


    slide_ui_initialized = reactive.Value(False)

    @reactive.effect
    def slide_reactivity():
        btn = input.slide_select_check()
        ui_initialized = slide_ui_initialized.get()

        if btn and not ui_initialized:
            dropdown_slide = ui.input_select("slide_select_drop", "Select the Slide Annotation", choices=obs_names.get())
            ui.insert_ui(
                ui.div({"id": "inserted-slide_dropdown"}, dropdown_slide),
                selector="#main-slide_dropdown",
                where="beforeEnd",
            )

            dropdown_label = ui.input_select("slide_select_label", "Select a Slide", choices=[])
            ui.insert_ui(
                ui.div({"id": "inserted-label_dropdown"}, dropdown_label),
                selector="#main-label_dropdown",
                where="beforeEnd",
            )
            slide_ui_initialized.set(True)

        elif not btn and ui_initialized:
            ui.remove_ui("#inserted-slide_dropdown")
            ui.remove_ui("#inserted-label_dropdown")
            slide_ui_initialized.set(False)

    @reactive.effect
    def update_slide_select_drop():
        adata = ad.AnnData(obs=obs_data.get())
        if input.slide_select_drop():
            selected_anno = input.slide_select_drop()
            labels = adata.obs[selected_anno].unique().tolist()
            ui.update_select("slide_select_label", choices=labels)

    region_ui_initialized = reactive.Value(False)

    @reactive.effect
    def region_reactivity():
        btn = input.region_select_check()
        ui_initialized = region_ui_initialized.get()

        if btn and not ui_initialized:
            dropdown_region = ui.input_select("region_select_drop", "Select the Region Annotation", choices=obs_names.get())
            ui.insert_ui(
                ui.div({"id": "inserted-region_dropdown"}, dropdown_region),
                selector="#main-region_dropdown",
                where="beforeEnd",
            )

            dropdown_label = ui.input_select("region_label_select", "Select a Region", choices=[])
            ui.insert_ui(
                ui.div({"id": "inserted-region_label_select_dropdown"}, dropdown_label),
                selector="#main-region_label_select_dropdown",
                where="beforeEnd",
            )
            region_ui_initialized.set(True)

        elif not btn and ui_initialized:
            ui.remove_ui("#inserted-region_dropdown")
            ui.remove_ui("#inserted-region_label_select_dropdown")
            region_ui_initialized.set(False)

    @reactive.effect
    def update_region_select_drop():
        adata = ad.AnnData(obs=obs_data.get())
        if input.region_select_drop():
            selected_anno = input.region_select_drop()
            labels = adata.obs[selected_anno].unique().tolist()
            ui.update_select("region_label_select", choices=labels)




    @output
    @render_widget
    @reactive.event(input.go_sp1, ignore_none=True)
    def spac_Spatial():
        adata = ad.AnnData(X=X_data.get(), obs=pd.DataFrame(obs_data.get()), obsm=obsm_data.get(), dtype=X_data.get().dtype)
        slide_check = input.slide_select_check()
        region_check = input.region_select_check()
        if adata is not None:
            if slide_check is False and region_check is False:
                out = spac.visualization.interative_spatial_plot(adata, annotations=input.spatial_anno(), figure_width=4, figure_height=4, dot_size=input.spatial_slider())
                out[0]['image_object'].update_xaxes(showticklabels=True, ticks="outside", tickwidth=2, ticklen=10)
                out[0]['image_object'].update_yaxes(showticklabels=True, ticks="outside", tickwidth=2, ticklen=10)
                return out[0]['image_object']
            if slide_check is True and region_check is False:

                adata_subset = adata[adata.obs[input.slide_select_drop()] == input.slide_select_label()].copy()
                out = spac.visualization.interative_spatial_plot(adata_subset, annotations=input.spatial_anno(), figure_width=4, figure_height=4, dot_size=input.spatial_slider())
                out[0]['image_object'].update_xaxes(showticklabels=True, ticks="outside", tickwidth=2, ticklen=10)
                out[0]['image_object'].update_yaxes(showticklabels=True, ticks="outside", tickwidth=2, ticklen=10)
                return out[0]['image_object']
            if slide_check is True and region_check is True:

                adata_subset = adata[(adata.obs[input.slide_select_drop()] == input.slide_select_label()) & (adata.obs[input.region_select_drop()] == input.region_label_select())].copy()
                out = spac.visualization.interative_spatial_plot(adata_subset, annotations=input.spatial_anno(), figure_width=4, figure_height=4, dot_size=input.spatial_slider())
                out[0]['image_object'].update_xaxes(showticklabels=True, ticks="outside", tickwidth=2, ticklen=10)
                out[0]['image_object'].update_yaxes(showticklabels=True, ticks="outside", tickwidth=2, ticklen=10)
                return out[0]['image_object']
            if slide_check is False and region_check is True:

                adata_subset = adata[(adata.obs[input.region_select_drop()] == input.region_label_select())].copy()
                out = spac.visualization.interative_spatial_plot(adata_subset, annotations=input.spatial_anno(), figure_width=4, figure_height=4, dot_size=input.spatial_slider())
                out[0]['image_object'].update_xaxes(showticklabels=True, ticks="outside", tickwidth=2, ticklen=10)
                out[0]['image_object'].update_yaxes(showticklabels=True, ticks="outside", tickwidth=2, ticklen=10)
                return out[0]['image_object']

        return None

    #@output
    #@render.plot
    #def spac_Neighborhood():
        #adata = ad.AnnData(X=X_data.get(), obs=pd.DataFrame(obs_data.get()))
        #if adata is not None:
            #out = spac.spatial_analysis.spatial_interaction(adata, annotation=input.neighbor_anno(), analysis_method=input.anno_method())
            #return out
        #return None

    @reactive.Calc
    def get_scatterplot_names():
        if obsm_names.get() is not None and var_names.get() is not None:
            obsm_list = obsm_names.get()
            var_list = var_names.get()
            obsm_dict = {item: item for item in obsm_list}
            features_dict = {item: item for item in var_list}
            dict = {"Annotated Tables" : obsm_dict, "Features" : features_dict}

            return dict
        return []

    @reactive.Calc
    def get_scatterplot_coordinates_x():
        adata = ad.AnnData(X=X_data.get(), var=pd.DataFrame(var_data.get()), obsm=obsm_data.get(), layers=layers_data.get())
        obsm = obsm_names.get()
        features = var_names.get()
        layer_selection = input.scatter_layer()
        selection = input.scatter_x()

        if selection in obsm:
            coords = adata.obsm[selection]
            x_coords = coords[:, 0]  # Extract the first column for x-coordinates
            return x_coords
        elif selection in features and layer_selection == "Original":
            column_index = adata.var_names.get_loc(selection)
            x_coords = adata.X[:, column_index]  # Extract the column corresponding to the feature
            return x_coords
        elif selection in features and layer_selection != "Original":
            column_index = adata.var_names.get_loc(selection)
            new_layer = adata.layers[layer_selection]
            x_coords = new_layer[:, column_index]  # Extract the column corresponding to the feature
            return x_coords

        return None



    @reactive.Calc
    def get_scatterplot_coordinates_y():
        adata = adata_main.get()
        obsm = obsm_names.get()
        features = var_names.get()
        layer_selection = input.scatter_layer()
        selection = input.scatter_y()

        if selection in obsm:
            coords = adata.obsm[selection]
            y_coords = coords[:, 1]  # Extract the second column for y-coordinates
            return y_coords
        elif selection in features and layer_selection == "Original":
            column_index = adata.var_names.get_loc(selection)
            y_coords = adata.X[:, column_index]  # Extract the column corresponding to the feature
            return y_coords
        elif selection in features and layer_selection != "Original":
            column_index = adata.var_names.get_loc(selection)
            new_layer = adata.layers[layer_selection]
            y_coords = new_layer[:, column_index]  # Extract the column corresponding to the feature
            return y_coords

        return None


    # Track the UI state for scatterplot dropdowns
    scatter_ui_initialized = reactive.Value(False)

    @reactive.effect
    def scatter_reactivity():
        btn = input.scatter_color_check()
        if btn and not scatter_ui_initialized.get():
            # Insert the color selection dropdown if not already initialized
            dropdown = ui.input_select("scatter_color", "Select Feature", choices=var_names.get())
            ui.insert_ui(
                ui.div({"id": "inserted-scatter_dropdown"}, dropdown),
                selector="#main-scatter_dropdown",
                where="beforeEnd",
            )
            scatter_ui_initialized.set(True)
        elif not btn and scatter_ui_initialized.get():
            # Remove the color selection dropdown if it exists
            ui.remove_ui("#inserted-scatter_dropdown")
            scatter_ui_initialized.set(False)

    @reactive.Calc
    def get_color_values():
        adata = ad.AnnData(X=X_data.get(), var=pd.DataFrame(var_data.get()))
        column_index = adata.var_names.get_loc(input.scatter_color())
        color_values = adata.X[:, column_index]
        return color_values



    @output
    @render.plot
    @reactive.event(input.go_scatter, ignore_none=True)
    def spac_Scatter():
        x_points = get_scatterplot_coordinates_x()
        y_points = get_scatterplot_coordinates_y()
        btn = input.scatter_color_check()
        if btn is False:
            fig, ax = spac.visualization.visualize_2D_scatter(x_points,y_points)
            return ax
        elif btn is True:
            fig1, ax1 = spac.visualization.visualize_2D_scatter(x_points,y_points, labels=get_color_values())
            return ax1





app = App(app_ui, server)


