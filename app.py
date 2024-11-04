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
from sag_py_execution_time_decorator import log_execution_time
import logging

app_ui = ui.page_fluid(


    ui.navset_tab_card(
        

    
        ui.nav_panel("Data Input",
                ui.div(
                {"style": "font-weight: bold; font-size: 30px;"},
                ui.p("SPAC Interactive Dashboard")),
                ui.input_file("input_file", "Choose a file to upload:", multiple=False),
                ui.output_text("print_rows"),
                ui.output_text("print_columns"),
                ui.output_text("print_obs_names"),
                ui.output_text("print_obsm_names"),
                ui.output_text("print_layers_names"),
                ui.output_text("print_uns_names")
            
            
        ), 
        ui.nav("Features",
    
            ui.card(
                ui.row(
                    ui.column(2,
                        ui.input_select("h1_feat", "Select a Feature", choices=[]),
                        ui.input_select("h1_layer", "Select a Table", choices=[], selected=["Original"]),
                        ui.input_checkbox("h1_group_by_check", "Group By", value=False),
                        ui.input_checkbox("h1_log_x", "Log X-axis", value=False),
                        ui.input_checkbox("h1_log_y", "Log Y-axis", value=False),
                        ui.div(id="main-h1_dropdown"),
                        ui.div(id="main-h1_check"),
                        ui.input_action_button("go_h1", "Render Plot", class_="btn-success")
                    ),
                    ui.column(10,
                        ui.output_plot("spac_Histogram_1")
                        
                    )
                ),
            )),
            ui.nav("Boxplots",
            ui.row(
                ui.column(6,
                    ui.card(
                        ui.column(12,
                            ui.input_select("bp1_anno", "Select an Annotation", choices=[]),
                            ui.input_select("bp1_layer", "Select a Table", choices=[], selected="Original"),
                            ui.input_selectize("bp1_features", "Select Features", multiple=True, choices=[], selected=[]),
                            ui.input_action_button("go_bp1", "Render Plot", class_="btn-success"),
                            ui.output_plot("spac_Boxplot_1")
                        )
                    ),
                ),
                ui.column(6,
                    ui.card(
                        ui.column(12,
                            ui.input_select("bp2_anno", "Select an Annotation", choices=[]),
                            ui.input_select("bp2_layer", "Select a Table", choices=[], selected="Original"),
                            ui.input_selectize("bp2_features", "Select Features", multiple=True, choices=[], selected=[]),
                            ui.input_action_button("go_bp2", "Render Plot", class_="btn-success"),
                            ui.output_plot("spac_Boxplot_2")
                        )
                    ),
                ),
            )
        ),
        ui.nav("Annotations",
            ui.card(
                ui.row(
                    ui.column(2,
                        ui.input_select("h2_anno", "Select an Annotation", choices=[]),
                        ui.input_checkbox("h2_group_by_check", "Group By", value=False),
                        ui.div(id="main-h2_dropdown"),
                        ui.div(id="main-h2_check"),
                        ui.input_action_button("go_h2", "Render Plot", class_="btn-success"),
                    ),
                    ui.column(10,
                        ui.output_plot("spac_Histogram_2")
                    )
                )
            )
        ),
        ui.nav("Feat. Vs Anno.",
            ui.card(
                ui.row(
                    ui.column(2,
                        ui.input_select("hm1_anno", "Select an Annotation", choices=[]),
                        ui.input_select("hm1_layer", "Select a Table", choices=[]),
                        ui.input_checkbox("dendogram", "Include Dendrogram", False),
                        ui.div(id="main-hm1_check"),
                        ui.div(id="main-hm2_check"),
                        ui.input_action_button("go_hm1", "Render Plot", class_="btn-success")
                    ),
                    ui.column(10,
                        ui.output_plot("spac_Heatmap")
                    )
                )
            )
        ),
        ui.nav("Anno. Vs Anno.",
            ui.card(
                ui.row(
                    ui.column(2,
                        ui.input_select("sk1_anno1", "Select Source Annotation", choices=[]),
                        ui.input_select("sk1_anno2", "Select Target Annotation", choices=[]),
                        ui.input_action_button("go_sk1", "Render Plot", class_="btn-success")
                    ),
                    ui.column(10,
                        output_widget("spac_Sankey")
                    )
                )
            ),
            ui.card(
                ui.row(
                    ui.column(2,
                        ui.input_select("rhm_anno1", "Select Source Annotation", choices=[], selected=[]),
                        ui.input_select("rhm_anno2", "Select Target Annotation", choices=[], selected=[]),
                        ui.input_action_button("go_rhm1", "Render Plot", class_="btn-success")
                    ),
                    ui.column(10,
                        output_widget("spac_Relational")
                    )
                )
            )
        ),
        ui.nav("Spatial",
            ui.card(
                ui.row(
                    ui.column(2,
                        ui.input_select("spatial_anno", "Select an Object", choices=[]),
                        ui.input_slider("spatial_slider", "Point Size", min=2, max=10, value=3),
                        ui.input_action_button("go_sp1", "Render Plot", class_="btn-success")
                    ),
                    ui.column(10,
                        output_widget("spac_Spatial")
                    )
                )
            )),
        ui.nav("UMAP",    
            ui.card(
                ui.row(
                    ui.column(6,
                        ui.input_radio_buttons("umap_rb", "Choose one:", ["Annotation", "Feature"]),
                        ui.input_select("plottype", "Select a plot type", choices=["umap", "pca", "tsne"]),
                        ui.div(id="main-ump_rb_dropdown_anno"),
                        ui.div(id="main-ump_rb_dropdown_feat"),
                        ui.div(id="main-ump_table_dropdown_feat"),
                        ui.input_slider("umap_slider_1", "Point Size", min=.5, max=10, value=3),
                        ui.input_action_button("go_umap1", "Render Plot", class_="btn-success"),
                        ui.output_plot("spac_UMAP")
                    ),
                    ui.column(6,
                        ui.input_radio_buttons("umap_rb2", "Choose one:", ["Annotation", "Feature"]),
                        ui.input_select("plottype2", "Select a plot type", choices=["umap", "pca", "tsne"]),
                        ui.div(id="main-ump_rb_dropdown_anno2"),
                        ui.div(id="main-ump_rb_dropdown_feat2"),
                        ui.div(id="main-ump_table_dropdown_feat2"),
                        ui.input_slider("umap_slider_2", "Point Size", min=.5, max=10, value=3),
                        ui.input_action_button("go_umap2", "Render Plot", class_="btn-success"),
                        ui.output_plot("spac_UMAP2")



                )
                
            )

            )
        ),
        ui.nav("Scatterplot",
            ui.card(
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
                        ui.output_plot("spac_Scatter")
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
            return "Obs: " + obs_str
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
            return "Obsm: " + obsm_str
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
            return "Layers: " + layers_str
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
            return "Uns: " + uns_str
        return
    
    @reactive.Calc
    @render.text 
    def print_rows():
        shape = shape_data.get()
        if shape is not None:
            return "# of Rows: " + str(shape[0])
        return 

    @reactive.Calc
    @render.text
    def print_columns():
        shape = shape_data.get()
        if shape is not None:
            return "# of Columns: " + str(shape[1])
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
                    return fig1
                else:
                    fig1 = spac.visualization.histogram(adata, feature=input.h1_feat(), log_scale=(btn_log_x, btn_log_y))
                    return fig1

            if input.h1_group_by_check() is not False:
                if input.h1_layer() != "Original":
                    fig1 = spac.visualization.histogram(adata, feature=input.h1_feat(), layer=input.h1_layer(), group_by=input.h1_anno(), together=input.h1_together_check(), log_scale=(btn_log_x, btn_log_y))
                    return fig1
                else:
                    fig1 = spac.visualization.histogram(adata, feature=input.h1_feat(), group_by=input.h1_anno(), together=input.h1_together_check(), log_scale=(btn_log_x, btn_log_y))
                    return fig1
        return None

    @reactive.effect
    def histogram_reactivity():
        btn = input.h1_group_by_check()
        if btn is True:
            dropdown = ui.input_select("h1_anno", "Select an Annotation", choices=obs_names.get())
            together_check = ui.input_checkbox("h1_together_check", "Plot Together", value=False)
            ui.insert_ui(
                ui.div({"id": "inserted-dropdown"}, dropdown),
                selector="#main-h1_dropdown",
                where="beforeEnd",
            )
            ui.insert_ui(
                ui.div({"id": "inserted-check"}, together_check),
                selector="#main-h1_check",
                where="beforeEnd",
            )
        elif btn is False:
            ui.remove_ui("#inserted-dropdown")
            ui.remove_ui("#inserted-check")

    @output
    @render.plot
    @reactive.event(input.go_bp1, ignore_none=True)
    def spac_Boxplot_1():
        adata = ad.AnnData(X=X_data.get(), obs=pd.DataFrame(obs_data.get()), var=pd.DataFrame(var_data.get()), layers=layers_data.get(), dtype=X_data.get().dtype)
        if adata is not None and adata.var is not None:
            if input.bp1_layer() != "Original" and input.bp1_anno() != "No Annotation":
                fig,ax = spac.visualization.boxplot(adata, annotation=input.bp1_anno(), layer=input.bp1_layer(), features=list(input.bp1_features()))
                return ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
            if input.bp1_layer() == "Original" and input.bp1_anno() != "No Annotation":
                fig,ax = spac.visualization.boxplot(adata, annotation=input.bp1_anno(), features=list(input.bp1_features()))
                return ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
            if input.bp1_layer() != "Original" and input.bp1_anno() == "No Annotation":
                fig,ax = spac.visualization.boxplot(adata, layer=input.bp1_layer(), features=list(input.bp1_features()))
                return ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
            if input.bp1_layer() == "Original" and input.bp1_anno() == "No Annotation":
                fig,ax = spac.visualization.boxplot(adata, features=list(input.bp1_features()))
                return ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
        return None

    @output
    @render.plot
    @reactive.event(input.go_bp2, ignore_none=True)
    def spac_Boxplot_2():
        adata = ad.AnnData(X=X_data.get(), obs=pd.DataFrame(obs_data.get()), var=pd.DataFrame(var_data.get()), layers=layers_data.get(), dtype=X_data.get().dtype)
        if adata is not None and adata.var is not None:
            if input.bp2_layer() != "Original" and input.bp2_anno() != "No Annotation":
                fig,ax = spac.visualization.boxplot(adata, annotation=input.bp2_anno(), layer=input.bp2_layer(), features=list(input.bp2_features()))
                return ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
            if input.bp2_layer() == "Original" and input.bp2_anno() != "No Annotation":
                fig,ax = spac.visualization.boxplot(adata, annotation=input.bp2_anno(), features=list(input.bp2_features()))
                return ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
            if input.bp2_layer() != "Original" and input.bp2_anno() == "No Annotation":
                fig,ax = spac.visualization.boxplot(adata, layer=input.bp2_layer(), features=list(input.bp2_features()))
                return ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
            if input.bp2_layer() == "Original" and input.bp2_anno() == "No Annotation":
                fig,ax = spac.visualization.boxplot(adata, features=list(input.bp2_features()))
                return ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
        return None

    

    @output
    @render.plot
    @reactive.event(input.go_h2, ignore_none=True)
    def spac_Histogram_2():
        adata = ad.AnnData(X=X_data.get(), obs=pd.DataFrame(obs_data.get()), var=pd.DataFrame(var_data.get()), layers=layers_data.get(), dtype=X_data.get().dtype)
        if adata is not None:
            if input.h2_group_by_check() is not False:
                fig1 = spac.visualization.histogram(adata, annotation=input.h2_anno(), group_by=input.h2_anno_1(), together=input.h2_together_check())
                return fig1
            else:
                fig = spac.visualization.histogram(adata, annotation=input.h2_anno())
                return fig
        return None    
    
    @reactive.effect
    def histogram_reactivity_2():
        btn = input.h2_group_by_check()
        if btn is not False:
            dropdown = ui.input_select("h2_anno_1", "Select an Annotation", choices=obs_names.get())
            together_check = ui.input_checkbox("h2_together_check", "Plot Together", value=False)
            ui.insert_ui(
                ui.div({"id": "inserted-dropdown-1"}, dropdown),
                selector="#main-h2_dropdown",
                where="beforeEnd",
            )
            ui.insert_ui(
                ui.div({"id": "inserted-check-1"}, together_check),
                selector="#main-h2_check",
                where="beforeEnd",
            )
        elif btn is not True:
            ui.remove_ui("#inserted-dropdown-1")
            ui.remove_ui("#inserted-check-1")

    @output
    @render.plot
    @reactive.event(input.go_hm1, ignore_none=True)
    def spac_Heatmap():
        adata = ad.AnnData(X=X_data.get(), obs=pd.DataFrame(obs_data.get()), var=pd.DataFrame(var_data.get()), layers=layers_data.get(), dtype=X_data.get().dtype)
        if adata is not None:
            if input.dendogram() is not True:
                if input.hm1_layer() != "Original":
                    df, fig, ax = spac.visualization.hierarchical_heatmap(adata, annotation=input.hm1_anno(), layer=input.hm1_layer(), z_score=None)
                    return fig
                else:
                    df, fig, ax = spac.visualization.hierarchical_heatmap(adata, annotation=input.hm1_anno(), layer=None, z_score=None)
                    return fig
            elif input.dendogram() is not False:
                cluster_annotations = input.h2_anno_dendro()  
                cluster_features = input.h2_feat_dendro()  
                if input.hm1_layer() != "Original":
                    df, fig, ax = spac.visualization.hierarchical_heatmap(adata, annotation=input.hm1_anno(), layer=input.hm1_layer(), z_score=None, cluster_annotations=cluster_annotations, cluster_feature=cluster_features)
                    return fig
                else:
                    df, fig, ax = spac.visualization.hierarchical_heatmap(adata, annotation=input.hm1_anno(), layer=None, z_score=None, cluster_annotations=cluster_annotations, cluster_feature=cluster_features)
                    return fig

        return None

    @reactive.effect
    def heatmap_reactivity():
        btn = input.dendogram()
        if btn is not False:
            annotation_check = ui.input_checkbox("h2_anno_dendro", "Annotation Cluster", value=False)
            feat_check = ui.input_checkbox("h2_feat_dendro", "Feature Cluster", value=False)
            ui.insert_ui(
                ui.div({"id": "inserted-check"}, annotation_check),
                selector="#main-hm1_check",
                where="beforeEnd",
            )
            ui.insert_ui(
                ui.div({"id": "inserted-check1"}, feat_check),
                selector="#main-hm2_check",
                where="beforeEnd",
            )
        elif btn is not True:
            ui.remove_ui("#inserted-check")
            ui.remove_ui("#inserted-check1")

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
            fig = spac.visualization.relational_heatmap(adata, source_annotation=input.rhm_anno1(), target_annotation=input.rhm_anno2())
            return fig
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
    
    @reactive.effect
    def umap_reactivity():
        flipper=data_loaded.get()
        if flipper is not False:
            btn = input.umap_rb()
            if btn == "Annotation":
                dropdown = ui.input_select("umap_rb_anno", "Select an Annotation", choices=obs_names.get())
                ui.insert_ui(
                    ui.div({"id": "inserted-rbdropdown_anno"}, dropdown),
                    selector="#main-ump_rb_dropdown_anno",
                    where="beforeEnd",
                )

                ui.remove_ui("#inserted-rbdropdown_feat")

            elif btn == "Feature":
                dropdown1 = ui.input_select("umap_rb_feat", "Select a Feature", choices=var_names.get())
                ui.insert_ui(
                    ui.div({"id": "inserted-rbdropdown_feat"}, dropdown1),
                    selector="#main-ump_rb_dropdown_feat",
                    where="beforeEnd",
                )
                new_choices = layers_names.get() + ["Original"]
                table_umap = ui.input_select("umap_layer", "Select a Table", choices=new_choices, selected=["Original"])
                ui.insert_ui(
                    ui.div({"id": "inserted-umap_table"}, table_umap),
                    selector="#main-ump_table_dropdown_feat",
                    where="beforeEnd",
                )
                ui.remove_ui("#inserted-rbdropdown_anno")
                
            elif btn == "None":
                ui.remove_ui("#inserted-rbdropdown_anno")
                ui.remove_ui("#inserted-rbdropdown_feat")
                ui.remove_ui("#inserted-umap_table")



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
    
    @reactive.effect
    def umap_reactivity2():
        flipper=data_loaded.get()
        if flipper is not False:
            btn = input.umap_rb2()
            if btn == "Annotation":
                dropdown = ui.input_select("umap_rb_anno2", "Select an Annotation", choices=obs_names.get())
                ui.insert_ui(
                    ui.div({"id": "inserted-rbdropdown_anno2"}, dropdown),
                    selector="#main-ump_rb_dropdown_anno2",
                    where="beforeEnd",
                )

                ui.remove_ui("#inserted-rbdropdown_feat2")
                ui.remove_ui("#inserted-umap_table2")

            elif btn == "Feature":
                dropdown1 = ui.input_select("umap_rb_feat2", "Select a Feature", choices=var_names.get())
                ui.insert_ui(
                    ui.div({"id": "inserted-rbdropdown_feat2"}, dropdown1),
                    selector="#main-ump_rb_dropdown_feat2",
                    where="beforeEnd",
                )
                new_choices = layers_names.get() + ["Original"]
                table_umap_1 = ui.input_select("umap_layer2", "Select a Table", choices=new_choices, selected=["Original"])
                ui.insert_ui(
                    ui.div({"id": "inserted-umap_table2"}, table_umap_1),
                    selector="#main-ump_table_dropdown_feat2",
                    where="beforeEnd",
                )
                ui.remove_ui("#inserted-rbdropdown_anno2")
                
            elif btn == "None":
                ui.remove_ui("#inserted-rbdropdown_anno2")
                ui.remove_ui("#inserted-rbdropdown_feat2")
                ui.remove_ui("#inserted-umap_table2")

                
    

    @output
    @render_widget
    @reactive.event(input.go_sp1, ignore_none=True)
    def spac_Spatial():
        adata = ad.AnnData(X=X_data.get(), obs=pd.DataFrame(obs_data.get()), obsm=obsm_data.get(), dtype=X_data.get().dtype)
        if adata is not None:
            out = spac.visualization.interative_spatial_plot(adata, annotations=input.spatial_anno(), figure_width=4, figure_height=4, dot_size=input.spatial_slider())
            out.update_xaxes(showticklabels=True, ticks="outside", tickwidth=2, ticklen=10)
            out.update_yaxes(showticklabels=True, ticks="outside", tickwidth=2, ticklen=10)
            return out
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
    

    @reactive.effect
    def _():
        btn = input.scatter_color_check()
        if btn is True:
            dropdown = ui.input_select("scatter_color", "Select Feature", choices=var_names.get())
            ui.insert_ui(
                ui.div({"id": "inserted-scatter_dropdown"}, dropdown),
                selector="#main-scatter_dropdown",
                where="beforeEnd",
            )

        elif btn is False:
            ui.remove_ui("#inserted-scatter_dropdown")

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


