from shiny import App, ui, reactive


# Screen imports
from ui import (
    data_input_ui,
    annotations_ui,
    features_ui,
    boxplot_ui,
    feat_vs_anno_ui,
    anno_vs_anno_ui,
    spatial_ui,
    umap_ui,
    scatterplot_ui,
    nearest_neighbor_ui,
    ripleyL_ui
)   

# Server imports
from server import (
    data_input_server,
    effect_update_server,
    annotations_server,
    features_server,
    boxplot_server,
    feat_vs_anno_server,
    anno_vs_anno_server,
    spatial_server,
    umap_server,
    scatterplot_server,
    nearest_neighbor_server,
    ripleyL_server
)

# Util imports
from utils.data_processing import load_data


file_path = "dev_example.pickle"  # Path to your preloaded .pickle file
preloaded_data = load_data(file_path)  # Initialize as None

app_ui = ui.page_fluid(
    ui.navset_card_tab(
        data_input_ui(),
        annotations_ui(),
        features_ui(),
        boxplot_ui(),
        feat_vs_anno_ui(),
        anno_vs_anno_ui(),
        spatial_ui(),
        umap_ui(),
        scatterplot_ui(),
        nearest_neighbor_ui(),
        ripleyL_ui(),
    )
)

def server(input, output, session):

    # Define a reactive variable to track if data is loaded
    data_loaded = reactive.Value(False)

    # Create a reactive variable for the main data
    adata_main = reactive.Value(preloaded_data)  # Initialize with preloaded data

    data_keys = [
        "X_data", 
        "obs_data", #AKA Annotations
        "obsm_data",
        "layers_data",
        "var_data", #AKA Features
        "uns_data",
        "shape_data",
        "obs_names",
        "obsm_names",
        "layers_names",
        "var_names",
        "uns_names",
        "df_heatmap",
        "df_relational",
        "df_boxplot",
        "df_histogram2",
        "df_histogram1",
        "df_ripley",
        "df_nn"
    ]

    shared = {
        "preloaded_data": preloaded_data,  # Preloaded data for initial load
        "data_loaded": data_loaded, # Reactive to track if data is loaded
        "adata_main": adata_main, # Main anndata object
    }

    # Dynamically create the reactive values for parts of the anndata object
    # and add them to the shared dictionary
    for key in data_keys:
        shared[key] = reactive.Value(None)

    # Individual server components
    data_input_server(input, output, session, shared)

    effect_update_server(input, output, session, shared)

    annotations_server(input, output, session, shared)

    features_server(input, output, session, shared)

    boxplot_server(input, output, session, shared)

    feat_vs_anno_server(input, output, session, shared)

    anno_vs_anno_server(input, output, session, shared)

    spatial_server(input, output, session, shared)

    umap_server(input, output, session, shared)

    scatterplot_server(input, output, session, shared)

    nearest_neighbor_server(input, output, session, shared)

    ripleyL_server(input, output, session, shared)


app = App(app_ui, server)
