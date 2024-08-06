# Import necessary modules
from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from shinywidgets import output_widget, render_widget
import pickle
import anndata as ad
import pandas as pd
import spac
import spac.visualization
import spac.spatial_analysis


# Define the UI for the application
app_ui = ui.page_fluid(
    ui.panel_title(ui.h2("SPAC Interactive Dashboard", class_="text-center")),
    ui.navset_tab_card(
        ui.nav_panel(
            "Data Input",
            ui.input_file("input_file", "Choose a file to upload:", multiple=False),
            ui.output_text("print_rows"),
            ui.output_text("print_columns"),
            ui.output_text("print_obs_names"),
            ui.output_text("print_obsm_names"),
            ui.output_text("print_layers_names"),
            ui.output_text("print_uns_names")
        ),
        ui.nav_panel(
            "Features",
            ui.card(
                ui.row(
                    ui.column(
                        2,
                        ui.input_select("h1_feat", "Select a Feature", choices=[]),
                        ui.input_select(
                            "h1_layer", "Select a Table", choices=[], 
                            selected=["Original"]
                        ),
                        ui.input_slider(
                            "h1_slider", "Select Number of Bins", min=1, max=200, 
                            value=20
                        ),
                        ui.input_checkbox("h1_group_by_check", "Group By", value=False),
                        ui.div(id="main-h1_dropdown"),
                        ui.div(id="main-h1_check"),
                        ui.input_action_button(
                            "go_h1", "Render Plot", class_="btn-success"
                        )
                    ),
                    ui.column(10, ui.output_plot("spac_Histogram_1"))
                )
            )
        ),
        ui.nav_panel(
            "Boxplots",
            ui.row(
                ui.column(
                    6,
                    ui.card(
                        ui.column(
                            12,
                            ui.input_select(
                                "bp1_anno", "Select an Annotation", choices=[]
                            ),
                            ui.input_select(
                                "bp1_layer", "Select a Table", choices=[], 
                                selected="Original"
                            ),
                            ui.input_selectize(
                                "bp1_features", "Select Features", multiple=True, 
                                choices=[], selected=[]
                            ),
                            ui.input_action_button(
                                "go_bp1", "Render Plot", class_="btn-success"
                            ),
                            ui.output_plot("spac_Boxplot_1")
                        )
                    )
                ),
                ui.column(
                    6,
                    ui.card(
                        ui.column(
                            12,
                            ui.input_select(
                                "bp2_anno", "Select an Annotation", choices=[]
                            ),
                            ui.input_select(
                                "bp2_layer", "Select a Table", choices=[], 
                                selected="Original"
                            ),
                            ui.input_selectize(
                                "bp2_features", "Select Features", multiple=True, 
                                choices=[], selected=[]
                            ),
                            ui.input_action_button(
                                "go_bp2", "Render Plot", class_="btn-success"
                            ),
                            ui.output_plot("spac_Boxplot_2")
                        )
                    )
                )
            )
        ),
        ui.nav_panel(
            "Annotations",
            ui.card(
                ui.row(
                    ui.column(
                        2,
                        ui.input_select("h2_anno", "Select an Annotation", choices=[]),
                        ui.input_action_button(
                            "go_h2", "Render Plot", class_="btn-success"
                        ),
                    ),
                    ui.column(10, ui.output_plot("spac_Histogram_2"))
                )
            )
        ),
        ui.nav_panel(
            "Feat. Vs Anno.",
            ui.card(
                ui.row(
                    ui.column(
                        2,
                        ui.input_select("hm1_anno", "Select an Annotation", choices=[]),
                        ui.input_select("hm1_layer", "Select a Table", choices=[]),
                        ui.input_checkbox("dendogram", "Include Dendogram", False),
                        ui.input_action_button(
                            "go_hm1", "Render Plot", class_="btn-success"
                        )
                    ),
                    ui.column(10, ui.output_plot("spac_Heatmap"))
                )
            )
        ),
        ui.nav_panel(
            "Anno. Vs Anno.",
            ui.card(
                ui.row(
                    ui.column(
                        2,
                        ui.input_select(
                            "sk1_anno1", "Select Source Annotation", choices=[]
                        ),
                        ui.input_select(
                            "sk1_anno2", "Select Target Annotation", choices=[]
                        ),
                        ui.input_action_button(
                            "go_sk1", "Render Plot", class_="btn-success"
                        )
                    ),
                    ui.column(10, output_widget("spac_Sankey"))
                )
            ),
            ui.card(
                ui.row(
                    ui.column(
                        2,
                        ui.input_select(
                            "rhm_anno1", "Select Source Annotation", choices=[], 
                            selected=[]
                        ),
                        ui.input_select(
                            "rhm_anno2", "Select Target Annotation", choices=[], 
                            selected=[]
                        ),
                        ui.input_action_button(
                            "go_rhm1", "Render Plot", class_="btn-success"
                        )
                    ),
                    ui.column(10, output_widget("spac_Relational"))
                )
            )
        ),
        ui.nav_panel(
            "Spatial",
            ui.card(
                ui.row(
                    ui.column(
                        2,
                        ui.input_select("spatial_anno", "Select an Object", choices=[]),
                        ui.input_slider(
                            "spatial_slider", "Point Size", min=2, max=10, value=3
                        ),
                        ui.input_action_button(
                            "go_sp1", "Render Plot", class_="btn-success"
                        )
                    ),
                    ui.column(10, output_widget("spac_Spatial"))
                )
            )
        ),
        ui.nav_panel(
            "UMAP",
            ui.card(
                ui.row(
                    ui.column(
                        6,
                        ui.input_radio_buttons(
                            "umap_rb", "Choose one:", ["Annotation", "Feature"]
                        ),
                        ui.input_select(
                            "plottype", "Select a plot type", 
                            choices=["umap", "pca", "tsne"]
                        ),
                        ui.div(id="main-ump_rb_dropdown_anno"),
                        ui.div(id="main-ump_rb_dropdown_feat"),
                        ui.input_action_button(
                            "go_umap1", "Render Plot", class_="btn-success"
                        ),
                        ui.output_plot("spac_UMAP")
                    ),
                    ui.column(
                        6,
                        ui.input_radio_buttons(
                            "umap_rb2", "Choose one:", ["Annotation", "Feature"]
                        ),
                        ui.input_select(
                            "plottype2", "Select a plot type", 
                            choices=["umap", "pca", "tsne"]
                        ),
                        ui.div(id="main-ump_rb_dropdown_anno2"),
                        ui.div(id="main-ump_rb_dropdown_feat2"),
                        ui.input_action_button(
                            "go_umap2", "Render Plot", class_="btn-success"
                        ),
                        ui.output_plot("spac_UMAP2")
                    )
                )
            )
        ),
        ui.nav_panel(
            "Scatterplot",
            ui.card(
                ui.row(
                    ui.column(
                        2,
                        ui.input_select(
                            "scatter_layer", "Select a Table", choices=[], 
                            selected="Original"
                        ),
                        ui.input_select("scatter_x", "Select X Axis", choices=[]),
                        ui.input_select("scatter_y", "Select Y Axis", choices=[]),
                        ui.input_checkbox(
                            "scatter_color_check", "Color by Feature", value=False
                        ),
                        ui.div(id="main-scatter_dropdown"),
                        ui.input_action_button(
                            "go_scatter", "Render Plot", class_="btn-success"
                        )
                    ),
                    ui.column(10, ui.output_plot("spac_Scatter"))
                )
            )
        )
    )
)

# Define the server logic for the application
def server(input: Inputs, output: Outputs, session: Session):
    # Define a reactive variable to track if data is loaded
    data_loaded = reactive.Value(False)

    @reactive.Effect
    def adata_filter():
        print("Calling Data")
        file_info = input.input_file()
        if not file_info:
            data_loaded.set(False)  # Set to False if no file is uploaded
            return
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
    X_data = reactive.Value(None)
    obs_data = reactive.Value(None)  # AKA Annotations
    obsm_data = reactive.Value(None)
    layers_data = reactive.Value(None)
    var_data = reactive.Value(None)  # AKA Features
    uns_data = reactive.Value(None)
    shape_data = reactive.Value(None)
    obs_names = reactive.Value(None)
    obsm_names = reactive.Value(None)
    layers_names = reactive.Value(None)
    var_names = reactive.Value(None)
    uns_names = reactive.Value(None)

    @reactive.Effect
    def update_parts():
        print("Updating Parts")
        adata = adata_main.get()
        if adata is not None:
            X_data.set(getattr(adata, 'X', None))
            obs_data.set(getattr(adata, 'obs', None))
            obsm_data.set(getattr(adata, 'obsm', None))
            layers_data.set(getattr(adata, 'layers', None))
            var_data.set(getattr(adata, 'var', None))
            uns_data.set(getattr(adata, 'uns', None))
            shape_data.set(adata.shape)

            obs_names.set(
                list(adata.obs.keys()) if hasattr(adata, 'obs') else None
            )
            obsm_names.set(
                list(adata.obsm.keys()) if hasattr(adata, 'obsm') else None
            )
            layers_names.set(
                list(adata.layers.keys()) if hasattr(adata, 'layers') else None
            )
            var_names.set(
                list(adata.var.index.tolist()) if hasattr(adata, 'var') else None
            )
            uns_names.set(
                list(adata.uns.keys()) if hasattr(adata, 'uns') else None
            )
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
        if obs:
            obs_str = ", ".join(obs)
            return f"Obs: {obs_str}"
        return ""

    @reactive.Calc
    @render.text
    def print_obsm_names():
        obsm = obsm_names.get()
        if obsm:
            obsm_str = ", ".join(obsm)
            return f"Obsm: {obsm_str}"
        return ""

    @reactive.Calc
    @render.text
    def print_layers_names():
        layers = layers_names.get()
        if layers:
            layers_str = ", ".join(layers)
            return f"Layers: {layers_str}"
        return ""

    @reactive.Calc
    @render.text
    def print_uns_names():
        uns = uns_names.get()
        if uns:
            uns_str = ", ".join(uns)
            return f"Uns: {uns_str}"
        return ""

    @reactive.Calc
    @render.text
    def print_rows():
        shape = shape_data.get()
        if shape:
            return f"# of Rows: {shape[0]}"
        return ""

    @reactive.Calc
    @render.text
    def print_columns():
        shape = shape_data.get()
        if shape:
            return f"# of Columns: {shape[1]}"
        return ""

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

    @reactive.Effect
    def update_select_input_layer():
        layers = layers_names.get()
        if layers is not None:
            new_choices = layers + ["Original"]
            ui.update_select("h1_layer", choices=new_choices)
            ui.update_select("bp1_layer", choices=new_choices)
            ui.update_select("bp2_layer", choices=new_choices)
            ui.update_select("hm1_layer", choices=new_choices)
            ui.update_select("scatter_layer", choices=new_choices)

    @reactive.Effect
    def update_select_input_anno_bp():
        annotations = obs_names.get()
        if annotations:
            new_choices = annotations + ["No Annotation"]
            ui.update_select("bp1_anno", choices=new_choices)
            ui.update_select("bp2_anno", choices=new_choices)

    @reactive.Effect
    def update_select_input_layer_scatter():
        choices = get_scatterplot_names()
        ui.update_select("scatter_x", choices=choices)
        ui.update_select("scatter_y", choices=choices)

    @reactive.Effect
    def update_boxplot_selectize():
        selected_names = var_names.get()
        if selected_names:
            ui.update_selectize("bp1_features", selected=selected_names[:2])
            ui.update_selectize("bp2_features", selected=selected_names[:2])

    @reactive.Effect
    def update_relational_select():
        selected_names = obs_names.get()
        if selected_names and len(selected_names) > 1:
            ui.update_selectize("rhm_anno1", selected=selected_names[0])
            ui.update_selectize("rhm_anno2", selected=selected_names[1])

    @output
    @render.plot
    @reactive.event(input.go_h1, ignore_none=True)
    def spac_Histogram_1():
        adata = ad.AnnData(
            X=X_data.get(), obs=pd.DataFrame(obs_data.get()), 
            var=pd.DataFrame(var_data.get()), layers=layers_data.get()
        )
        if adata:
            if not input.h1_group_by_check():
                fig1 = spac.visualization.histogram(
                    adata, feature=input.h1_feat(), bins=input.h1_slider(), 
                    layer=input.h1_layer() if input.h1_layer() != "Original" else None
                )
                return fig1
            if input.h1_group_by_check():
                fig = spac.visualization.histogram(
                    adata, feature=input.h1_feat(), bins=input.h1_slider(), 
                    group_by=input.h1_anno(), together=input.h1_together_check()
                )
                return fig
        return None

    @reactive.Effect
    def histogram_reactivity():
        if input.h1_group_by_check():
            dropdown = ui.input_select(
                "h1_anno", "Select an Annotation", choices=obs_names.get()
            )
            together_check = ui.input_checkbox(
                "h1_together_check", "Plot Together", value=False
            )
            ui.insert_ui(
                ui.div({"id": "inserted-dropdown"}, dropdown), 
                selector="#main-h1_dropdown", where="beforeEnd"
            )
            ui.insert_ui(
                ui.div({"id": "inserted-check"}, together_check), 
                selector="#main-h1_check", where="beforeEnd"
            )
        else:
            ui.remove_ui("#inserted-dropdown")
            ui.remove_ui("#inserted-check")

    @output
    @render.plot
    @reactive.event(input.go_bp1, ignore_none=True)
    def spac_Boxplot_1():
        adata = ad.AnnData(
            X=X_data.get(), obs=pd.DataFrame(obs_data.get()), 
            var=pd.DataFrame(var_data.get()), layers=layers_data.get()
        )
        if adata and adata.var is not None:
            annotation = input.bp1_anno()
            layer = input.bp1_layer()
            features = list(input.bp1_features())
            fig, ax = spac.visualization.boxplot(
                adata, annotation=annotation if annotation != "No Annotation" else None, 
                layer=layer if layer != "Original" else None, features=features
            )
            return ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
        return None

    @output
    @render.plot
    @reactive.event(input.go_bp2, ignore_none=True)
    def spac_Boxplot_2():
        adata = ad.AnnData(
            X=X_data.get(), obs=pd.DataFrame(obs_data.get()), 
            var=pd.DataFrame(var_data.get()), layers=layers_data.get()
        )
        if adata and adata.var is not None:
            annotation = input.bp2_anno()
            layer = input.bp2_layer()
            features = list(input.bp2_features())
            fig, ax = spac.visualization.boxplot(
                adata, annotation=annotation if annotation != "No Annotation" else None, 
                layer=layer if layer != "Original" else None, features=features
            )
            return ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
        return None

    @output
    @render.plot
    @reactive.event(input.go_h2, ignore_none=True)
    def spac_Histogram_2():
        adata = ad.AnnData(X=X_data.get(), obs=pd.DataFrame(obs_data.get()))
        if adata:
            fig = spac.visualization.histogram(adata, annotation=input.h2_anno())
            return fig
        return None

    @output
    @render.plot
    @reactive.event(input.go_hm1, ignore_none=True)
    def spac_Heatmap():
        adata = ad.AnnData(
            X=X_data.get(), obs=pd.DataFrame(obs_data.get()), 
            var=pd.DataFrame(var_data.get()), layers=layers_data.get()
        )
        if adata:
            df, fig, ax = spac.visualization.hierarchical_heatmap(
                adata, annotation=input.hm1_anno(), 
                layer=input.hm1_layer() if input.hm1_layer() != "Original" else None, 
                z_score=None, cluster_annotations=input.dendogram()
            )
            return fig
        return None

    @output
    @render_widget
    @reactive.event(input.go_sk1, ignore_none=True)
    def spac_Sankey():
        adata = ad.AnnData(
            X=X_data.get(), obs=pd.DataFrame(obs_data.get()), 
            layers=layers_data.get()
        )
        if adata:
            fig = spac.visualization.sankey_plot(
                adata, source_annotation=input.sk1_anno1(), 
                target_annotation=input.sk1_anno2()
            )
            return fig
        return None

    @output
    @render_widget
    @reactive.event(input.go_rhm1, ignore_none=True)
    def spac_Relational():
        adata = ad.AnnData(X=X_data.get(), obs=pd.DataFrame(obs_data.get()))
        if adata:
            fig = spac.visualization.relational_heatmap(
                adata, source_annotation=input.rhm_anno1(), 
                target_annotation=input.rhm_anno2()
            )
            return fig
        return None

    @output
    @render.plot
    @reactive.event(input.go_umap1, ignore_none=True)
    def spac_UMAP():
        adata = ad.AnnData(
            X=X_data.get(), var=pd.DataFrame(var_data.get()), 
            obsm=obsm_data.get(), obs=obs_data.get()
        )
        if adata:
            method = input.plottype()
            if input.umap_rb() == "Feature":
                return spac.visualization.dimensionality_reduction_plot(
                    adata, method=method, feature=input.umap_rb_feat(), point_size=50
                )
            if input.umap_rb() == "Annotation":
                return spac.visualization.dimensionality_reduction_plot(
                    adata, method=method, annotation=input.umap_rb_anno(), 
                    point_size=100
                )
        return None

    @reactive.Effect
    def umap_reactivity():
        if data_loaded.get():
            btn = input.umap_rb()
            if btn == "Annotation":
                dropdown = ui.input_select(
                    "umap_rb_anno", "Select an Annotation", choices=obs_names.get()
                )
                ui.insert_ui(
                    ui.div({"id": "inserted-rbdropdown_anno"}, dropdown), 
                    selector="#main-ump_rb_dropdown_anno", where="beforeEnd"
                )
                ui.remove_ui("#inserted-rbdropdown_feat")
            elif btn == "Feature":
                dropdown = ui.input_select(
                    "umap_rb_feat", "Select a Feature", choices=var_names.get()
                )
                ui.insert_ui(
                    ui.div({"id": "inserted-rbdropdown_feat"}, dropdown), 
                    selector="#main-ump_rb_dropdown_feat", where="beforeEnd"
                )
                ui.remove_ui("#inserted-rbdropdown_anno")
            elif btn == "None":
                ui.remove_ui("#inserted-rbdropdown_anno")
                ui.remove_ui("#inserted-rbdropdown_feat")

    @output
    @render.plot
    @reactive.event(input.go_umap2, ignore_none=True)
    def spac_UMAP2():
        adata = ad.AnnData(
            X=X_data.get(), var=pd.DataFrame(var_data.get()), 
            obsm=obsm_data.get(), obs=obs_data.get()
        )
        if adata:
            method = input.plottype2()
            if input.umap_rb2() == "Feature":
                return spac.visualization.dimensionality_reduction_plot(
                    adata, method=method, feature=input.umap_rb_feat2(), point_size=50
                )
            if input.umap_rb2() == "Annotation":
                return spac.visualization.dimensionality_reduction_plot(
                    adata, method=method, annotation=input.umap_rb_anno2(), 
                    point_size=100
                )
        return None

    @reactive.Effect
    def umap_reactivity2():
        if data_loaded.get():
            btn = input.umap_rb2()
            if btn == "Annotation":
                dropdown = ui.input_select(
                    "umap_rb_anno2", "Select an Annotation", choices=obs_names.get()
                )
                ui.insert_ui(
                    ui.div({"id": "inserted-rbdropdown_anno2"}, dropdown), 
                    selector="#main-ump_rb_dropdown_anno2", where="beforeEnd"
                )
                ui.remove_ui("#inserted-rbdropdown_feat2")
            elif btn == "Feature":
                dropdown = ui.input_select(
                    "umap_rb_feat2", "Select a Feature", choices=var_names.get()
                )
                ui.insert_ui(
                    ui.div({"id": "inserted-rbdropdown_feat2"}, dropdown), 
                    selector="#main-ump_rb_dropdown_feat2", where="beforeEnd"
                )
                ui.remove_ui("#inserted-rbdropdown_anno2")
            elif btn == "None":
                ui.remove_ui("#inserted-rbdropdown_anno2")
                ui.remove_ui("#inserted-rbdropdown_feat2")

    @output
    @render_widget
    @reactive.event(input.go_sp1, ignore_none=True)
    def spac_Spatial():
        adata = ad.AnnData(
            X=X_data.get(), obs=pd.DataFrame(obs_data.get()), obsm=obsm_data.get()
        )
        if adata:
            out = spac.visualization.interative_spatial_plot(
                adata, annotations=input.spatial_anno(), figure_width=4, 
                figure_height=4, dot_size=input.spatial_slider()
            )
            return out
        return None

    @reactive.Calc
    def get_scatterplot_names():
        if obsm_names.get() and var_names.get():
            obsm_dict = {item: item for item in obsm_names.get()}
            features_dict = {item: item for item in var_names.get()}
            return {"Annotated Tables": obsm_dict, "Features": features_dict}
        return []

    @reactive.Calc
    def get_scatterplot_coordinates_x():
        adata = ad.AnnData(
            X=X_data.get(), var=pd.DataFrame(var_data.get()), 
            obsm=obsm_data.get(), layers=layers_data.get()
        )
        obsm = obsm_names.get()
        features = var_names.get()
        layer_selection = input.scatter_layer()
        selection = input.scatter_x()

        if selection in obsm:
            coords = adata.obsm[selection]
            return coords[:, 0]  # Extract the first column for x-coordinates
        if selection in features:
            column_index = adata.var_names.get_loc(selection)
            data = adata.X if layer_selection == "Original" else adata.layers[
                layer_selection]
            return data[:, column_index]  # Extract the column for the feature
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
            return coords[:, 1]  # Extract the second column for y-coordinates
        if selection in features:
            column_index = adata.var_names.get_loc(selection)
            data = adata.X if layer_selection == "Original" else adata.layers[
                layer_selection]
            return data[:, column_index]  # Extract the column for the feature
        return None

    @reactive.Effect
    def _():
        if input.scatter_color_check():
            dropdown = ui.input_select(
                "scatter_color", "Select Feature", choices=var_names.get()
            )
            ui.insert_ui(
                ui.div({"id": "inserted-scatter_dropdown"}, dropdown), 
                selector="#main-scatter_dropdown", where="beforeEnd"
            )
        else:
            ui.remove_ui("#inserted-scatter_dropdown")

    @reactive.Calc
    def get_color_values():
        adata = ad.AnnData(X=X_data.get(), var=pd.DataFrame(var_data.get()))
        column_index = adata.var_names.get_loc(input.scatter_color())
        return adata.X[:, column_index]

    @output
    @render.plot
    @reactive.event(input.go_scatter, ignore_none=True)
    def spac_Scatter():
        x_points = get_scatterplot_coordinates_x()
        y_points = get_scatterplot_coordinates_y()
        if not input.scatter_color_check():
            fig, ax = spac.visualization.visualize_2D_scatter(x_points, y_points)
            return ax
        fig1, ax1 = spac.visualization.visualize_2D_scatter(
            x_points, y_points, labels=get_color_values()
        )
        return ax1

# Run the application
app = App(app_ui, server)
