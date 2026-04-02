from shiny import App, ui, reactive
import os

#Screen Imports
from ui import (
    getting_started_ui,
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
#Server Imports
from server import (
    getting_started_server,
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
#Utils Imports
from utils.data_processing import load_data, read_html_file
from utils.accessibility import accessible_navigation, apply_slider_accessibility_global
from utils.security import apply_security_enhancements

header_html = read_html_file("header.html")
footer_html = read_html_file("footer.html")
file_path = "dev_example.pickle"
preloaded_data = load_data(file_path)

app_ui = ui.page_fluid(
    apply_security_enhancements(),
    ui.HTML(header_html),
    accessible_navigation(),
    apply_slider_accessibility_global(),
    ui.navset_card_tab(
        getting_started_ui(),
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
    ),
    ui.HTML(footer_html)
)


def cancel_plot(plot_id, shared):
    """
    Global cancel function. Terminates process or clears thread state
    for the given plot_id, then cleans up all associated reactive values.
    """
    registry = shared['plot_registry']

    if plot_id not in registry:
        return

    entry = registry[plot_id]

    # Terminate if it's a process
    if entry['type'] == 'process':
        p = entry['proc']
        if p is not None and p.is_alive():
            p.terminate()
            p.join(timeout=1)
            if p.is_alive():
                p.kill()
    if entry.get('current_proc_reactive'):
        entry['current_proc_reactive'].set(None)

    # Clear reactive state
    entry['is_calculating'].set(False)
    entry['result'].set(None)

    # Clear associated shared data key if present
    if entry.get('data_key') and shared.get(entry['data_key']):
        shared[entry['data_key']].set(None)

    del registry[plot_id]

    ui.notification_show("Render Cancelled", type="warning")


def server(input, output, session):
    data_loaded = reactive.Value(False)
    adata_main = reactive.Value(preloaded_data)

    data_keys = [
        "X_data",
        "obs_data",
        "obsm_data",
        "layers_data",
        "var_data",
        "uns_data",
        "shape_data",
        "obs_names",
        "obsm_names",
        "layers_names",
        "var_names",
        "uns_names",
        "spatial_distance_columns",
        "df_heatmap",
        "df_relational",
        "df_boxplot",
        "df_histogram2",
        "df_histogram1",
        "df_nn",
        "df_ripley"
    ]

    shared = {
        "preloaded_data": preloaded_data,
        "data_loaded": data_loaded,
        "adata_main": adata_main,
        "plot_registry": {},
    }

    for key in data_keys:
        shared[key] = reactive.Value(None)

    # Individual server components
    getting_started_server(input, output, session, shared)
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

    # Global cancel handlers
    STOP_BUTTON_MAP = {
        'stop_h2':'histogram',
        'stop_scatterplot':'scatterplot',
        'stop_umap1':'umap1',
        'stop_umap2':'umap2',
        'stop_spatial':'spatial',
        'stop_features':'features',
        'stop_boxplot':'boxplot',
        'stop_sankey':'sankey',
        'stop_relational':'relational',
        'stop_nn':'nn',
        'stop_rl':'rl',
        'stop_hm1':'hm1',
    }

    for stop_btn, plot_id in STOP_BUTTON_MAP.items():
        def make_cancel_effect(btn=stop_btn, pid=plot_id):
            @reactive.Effect
            @reactive.event(input[btn])
            def _cancel_effect():
                cancel_plot(pid, shared)
        make_cancel_effect()


static_path = os.path.join(os.path.dirname(__file__), "www")
app = App(app_ui, server, static_assets=static_path)