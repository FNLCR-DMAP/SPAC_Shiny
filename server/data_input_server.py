from shiny import render, reactive
import pickle
import anndata as ad


def data_input_server(input, output, session, shared):
    @reactive.Effect
    def adata_filter():
        print("Calling Data")
        file_info = input.input_file()
        if not file_info:
            # Only set preloaded data if it exists
            if shared['preloaded_data'] is None:
                shared['data_loaded'].set(False)
            else:
                shared['adata_main'].set(shared['preloaded_data'])
                shared['data_loaded'].set(True)
        else:
            file_path = file_info[0]['datapath']
            with open(file_path, 'rb') as file:
                if file_path.endswith('.pickle'):
                    shared['adata_main'].set(pickle.load(file))
                elif file_path.endswith('.h5ad'):
                    shared['adata_main'].set(ad.read_h5ad(file_path))
                else:
                    shared['adata_main'].set(ad.read(file_path))
            # Set to True if a file is successfully uploaded
            shared['data_loaded'].set(True)

    @reactive.Effect
    def update_parts():
        print("Updating Parts")
        adata = shared['adata_main'].get()
        if adata is None:
            shared['obs_data'].set(None)
            shared['obsm_data'].set(None)
            shared['layers_data'].set(None)
            shared['var_data'].set(None)
            shared['uns_data'].set(None)
            shared['shape_data'].set(None)
            shared['obs_names'].set(None)
            shared['obsm_names'].set(None)
            shared['layers_names'].set(None)
            shared['var_names'].set(None)
            shared['uns_names'].set(None)
            shared['spatial_distance_columns'].set(None)

        else:
           if hasattr(adata, 'X'):
               shared['X_data'].set(adata.X)
           else:
               shared['X_data'].set(None)

           if hasattr(adata, 'obs'):
               shared['obs_data'].set(adata.obs)
           else:
               shared['obs_data'].set(None)

           if hasattr(adata, 'obsm'):
               shared['obsm_data'].set(adata.obsm)
           else:
               shared['obsm_data'].set(None)

           if hasattr(adata, 'layers'):
               shared['layers_data'].set(adata.layers)
           else:
               shared['layers_data'].set(None)

           if hasattr(adata, 'var'):
               shared['var_data'].set(adata.var)
           else:
               shared['var_data'].set(None)

           if hasattr(adata, 'uns'):
               shared['uns_data'].set(adata.uns)
           else:
               shared['uns_data'].set(None)

           shared['shape_data'].set(adata.shape)

           if hasattr(adata, 'obs'):
               shared['obs_names'].set(list(adata.obs.keys()))
           else:
               shared['obs_names'].set(None)

           if hasattr(adata, 'obsm'):
               shared['obsm_names'].set(list(adata.obsm.keys()))
           else:
               shared['obsm_names'].set(None)

           if hasattr(adata, 'layers'):
               shared['layers_names'].set(list(adata.layers.keys()))
           else:
               shared['layers_names'].set(None)

           if hasattr(adata, 'var'):
               shared['var_names'].set(list(adata.var.index.tolist()))
           else:
               shared['var_names'].set(None)

           if hasattr(adata, 'uns'):
               shared['uns_names'].set(list(adata.uns.keys()))
           else:
               shared['uns_names'].set(None)

           # Extract spatial_distance column names if available via helper
           from utils.data_processing import get_spatial_distance_columns

           spatial_cols = get_spatial_distance_columns(adata)
           shared['spatial_distance_columns'].set(spatial_cols)


    @reactive.Calc
    @render.text
    def print_obs_names():
        obs = shared['obs_names'].get()
        if not obs:
            return "Annotations: None"
        if obs is None:
            return "Empty"
        else:
             if len(obs) > 1:
                 obs_str = ", ".join(obs)
             else:
                 obs_str = obs[0] if obs else ""
             return "Annotations: " + obs_str
        return

    @reactive.Calc
    @render.text
    def print_obsm_names():
        obsm = shared['obsm_names'].get()
        if not obsm:
            return "Associated Tables: None"
        if obsm is None:
            return "Empty"
        else:

            if len(obsm) > 1:
                obsm_str = ", ".join(obsm)
            else:
                obsm_str = obsm[0] if obsm else ""
            return "Associated Tables: " + obsm_str
        return

    @reactive.Calc
    @render.text
    def print_layers_names():
        layers = shared['layers_names'].get()
        # If there are no layers at all, just say "None"
        if not layers:
            return "Tables: None"
        # If there's more than one layer
        if len(layers) > 1:
            layers_str = ", ".join(layers)
        # If there's exactly one layer
        else:
            layers_str = layers[0]
        return "Tables: " + layers_str

    @reactive.Calc
    @render.text
    def print_uns_names():
        uns = shared['uns_names'].get()
        if not uns:
            return "Unstructured Data: None"
        if uns is None:
            return
        else:
            if len(uns) > 1:
                uns_str = ", ".join(uns)
            else:
                uns_str = uns[0] if uns else ""
            return "Unstructured Data: " + uns_str
        return

    @reactive.Calc
    @render.text
    def print_rows():
        shape = shared['shape_data'].get()
        if not shape:
            return "None"
        if shape is None:
            return "Empty"
        else:
            return str(shape[0])
        return

    @reactive.Calc
    @render.text
    def print_columns():
        shape = shared['shape_data'].get()
        if not shape:
            return "None"
        if shape is None:
            return "Empty"
        else:
            return  str(shape[1])
        return

    # Formatted UI outputs for better display
    @reactive.Calc
    @render.ui
    def formatted_obs_names():
        from shiny import ui
        obs = shared['obs_names'].get()
        if not obs:
            return ui.div("No annotations available", class_="text-muted")
        if obs is None:
           return ui.div("No data loaded", class_="text-muted")
        elif len(obs) > 0:
            items = [
                ui.div(
                    ui.span("•", class_="metric-bullet"),
                    ui.span(name, class_="metric-text"),
                    class_="metric-item"
                ) for name in obs
            ]
            return ui.div(*items)

    @reactive.Calc
    @render.ui
    def formatted_obsm_names():
        from shiny import ui
        obsm = shared['obsm_names'].get()
        if not obsm:
            return ui.div("No associated tables available", class_="text-muted")
        if obsm is None:
            return ui.div("No data loaded", class_="text-muted")
        elif len(obsm) > 0:
            items = [
                ui.div(
                    ui.span("•", class_="metric-bullet"),
                    ui.span(name, class_="metric-text"),
                    class_="metric-item"
                ) for name in obsm
            ]
            return ui.div(*items)

    @reactive.Calc
    @render.ui
    def formatted_layers_names():
        from shiny import ui
        layers = shared['layers_names'].get()
        if not layers:
            return ui.div("No tables available", class_="text-muted")
        if layers is None:
            return ui.div("No data loaded", class_="text-muted")
            items = [
                ui.div(
                    ui.span("•", class_="metric-bullet"),
                    ui.span(name, class_="metric-text"),
                    class_="metric-item"
                ) for name in layers
            ]
            return ui.div(*items)
        elif len(layers) > 0:
            items = [
                ui.div(
                    ui.span("•", class_="metric-bullet"),
                    ui.span(name, class_="metric-text"),
                    class_="metric-item"
                ) for name in layers
            ]
            return ui.div(*items)

    @reactive.Calc
    @render.ui
    def formatted_uns_names():
        from shiny import ui
        uns = shared['uns_names'].get()
        if not uns:
            return ui.div("No unstructured data available", class_="text-muted")
        if uns is None:
            return ui.div("No data loaded", class_="text-muted")
        elif len(uns) > 0:
           items = [
               ui.div(
                   ui.span("•", class_="metric-bullet"),
                   ui.span(name, class_="metric-text"),
                   class_="metric-item"
               ) for name in uns
           ]
           return ui.div(*items)