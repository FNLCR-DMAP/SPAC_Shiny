from shiny import App, Inputs, Outputs, Session, render, ui, reactive
import pickle
import anndata as ad


app_ui = ui.page_fluid(


    ui.navset_tab_card(
        

    
        ui.nav_panel("Data Input",
                ui.input_file("input_file", "Choose a file to upload:", multiple=False),
                ui.output_text("print_rows"),
                ui.output_text("print_columns"),
                ui.output_text("print_obs_names"),
                ui.output_text("print_obsm_names"),
                ui.output_text("print_layers_names"),
                ui.output_text("print_uns_names")
            
            
        )
    )
)


def server(input, output, session):



    # Update the reactive variable when a new file is uploaded
    @reactive.Effect
    def adata_filter():
        with ui.Progress(min=1, max=10) as p:
            p.set(message="Calculation in progress", detail="This may take a while...")
            file_info = input.input_file()
            if file_info:
                file_path = file_info[0]['datapath']
                with open(file_path, 'rb') as file:
                    if file_path.endswith('.pickle'):
                        adata_main.set(pickle.load(file))
                    elif file_path.endswith('.h5ad'):
                        adata_main.set(ad.read_h5ad(file_path))
                    else:
                        adata_main.set(ad.read(file_path))


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
    

    


app = App(app_ui, server)


