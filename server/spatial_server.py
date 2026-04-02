from shiny import ui, reactive
from shinywidgets import output_widget, render_widget
import anndata as ad
import pandas as pd
import spac.visualization


def spatial_server(input, output, session, shared):
    slide_ui_initialized = reactive.Value(False)

    @reactive.effect
    def slide_reactivity():
        btn = input.slide_select_check()
        ui_initialized = slide_ui_initialized.get()

        if btn and not ui_initialized:
            dropdown_slide = ui.input_select(
                "slide_select_drop",
                "Select the Slide Annotation",
                choices=shared['obs_names'].get())
            ui.insert_ui(
                ui.div({"id": "inserted-slide_dropdown"}, dropdown_slide),
                selector="#main-slide_dropdown",
                where="beforeEnd",
            )

            dropdown_label = ui.input_select(
                "slide_select_label",
                "Select a Slide",
                choices=[]
            )
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
        adata = ad.AnnData(obs=shared['obs_data'].get())
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
            dropdown_region = ui.input_select(
                "region_select_drop",
                "Select the Region Annotation",
                choices=shared['obs_names'].get())
            ui.insert_ui(
                ui.div({"id": "inserted-region_dropdown"}, dropdown_region),
                selector="#main-region_dropdown",
                where="beforeEnd",
            )

            dropdown_label = ui.input_select(
                "region_label_select",
                "Select a Region",
                choices=[]
            )
            ui.insert_ui(
                ui.div(
                    {"id": "inserted-region_label_select_dropdown"},
                    dropdown_label
                ),
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
        adata = ad.AnnData(obs=shared['obs_data'].get())
        if input.region_select_drop():
            selected_anno = input.region_select_drop()
            labels = adata.obs[selected_anno].unique().tolist()
            ui.update_select("region_label_select", choices=labels)

    @output
    @render_widget
    @reactive.event(input.go_sp1, ignore_none=True)
    def spac_Spatial():
        x_data = shared['X_data'].get()
        if x_data is None:
            return None

        cache = shared['cache']
        version = shared['dataset_version'].get()

        slide_check = input.slide_select_check()
        region_check = input.region_select_check()
        mode = input.spatial_rb()

        # Safely read dynamic inputs that may not be rendered yet
        try:
            slide_drop = input.slide_select_drop() if slide_check else None
            slide_label = input.slide_select_label() if slide_check else None
        except Exception:
            slide_drop, slide_label = None, None

        try:
            region_drop = input.region_select_drop() if region_check else None
            region_label = input.region_label_select() if region_check else None
        except Exception:
            region_drop, region_label = None, None

        try:
            annotation = input.spatial_anno() if mode == "Annotation" else None
        except Exception:
            annotation = None

        try:
            feature = input.spatial_feat() if mode == "Feature" else None
            sp_layer = (
                None if input.spatial_layer() == "Original"
                else input.spatial_layer()
            ) if mode == "Feature" else None
        except Exception:
            feature, sp_layer = None, None

        dot_size = input.spatial_slider()

        params = {
            'slide_check': slide_check,
            'slide_drop': slide_drop,
            'slide_label': slide_label,
            'region_check': region_check,
            'region_drop': region_drop,
            'region_label': region_label,
            'mode': mode,
            'annotation': annotation,
            'feature': feature,
            'layer': sp_layer,
            'dot_size': dot_size,
        }

        def compute():
            adata = ad.AnnData(
                X=x_data,
                var=pd.DataFrame(shared['var_data'].get()),
                obsm=shared['obsm_data'].get(),
                obs=shared['obs_data'].get(),
                dtype=x_data.dtype,
                layers=shared['layers_data'].get()
            )

            # Apply slide / region subsetting
            if not slide_check and not region_check:
                adata_subset = adata
            elif slide_check and not region_check:
                adata_subset = adata[
                    adata.obs[slide_drop] == slide_label
                ].copy()
            elif slide_check and region_check:
                adata_subset = adata[
                    (adata.obs[slide_drop] == slide_label) &
                    (adata.obs[region_drop] == region_label)
                ].copy()
            elif not slide_check and region_check:
                adata_subset = adata[
                    adata.obs[region_drop] == region_label
                ].copy()
            else:
                return None, None

            if mode == "Feature":
                if not feature:
                    return None, None
                out = spac.visualization.interactive_spatial_plot(
                    adata_subset,
                    feature=feature,
                    layer=sp_layer,
                    figure_width=5.5,
                    figure_height=5,
                    dot_size=dot_size
                )
            elif mode == "Annotation":
                if not annotation:
                    return None, None
                out = spac.visualization.interactive_spatial_plot(
                    adata_subset,
                    annotations=annotation,
                    figure_width=5.5,
                    figure_height=5,
                    dot_size=dot_size
                )
            else:
                return None, None

            fig = out[0]['image_object']
            fig.update_xaxes(
                showticklabels=True,
                ticks="outside",
                tickwidth=2,
                ticklen=10
            )
            fig.update_yaxes(
                showticklabels=True,
                ticks="outside",
                tickwidth=2,
                ticklen=10
            )
            return fig, None

        fig, _ = cache.get_or_compute('spatial', version, params, compute)
        return fig

    #Track UI State
    spatial_annotation_initialized = reactive.Value(False)
    spatial_feature_initialized = reactive.Value(False)

    @reactive.effect
    def spatial_reactivity():
        flipper = shared['data_loaded'].get()
        if flipper is not False:
            btn = input.spatial_rb()

            if btn == "Annotation":
                if not spatial_annotation_initialized.get():
                    dropdown = ui.input_select(
                        "spatial_anno", "Select an Annotation",
                        choices=shared['obs_names'].get()
                    )
                    ui.insert_ui(
                        ui.div(
                            {"id": "inserted-spatial_dropdown_anno"},
                            dropdown
                        ),
                        selector="#main-spatial_dropdown_anno",
                        where="beforeEnd"
                    )
                    spatial_annotation_initialized.set(True)

                if spatial_feature_initialized.get():
                    ui.remove_ui("#inserted-spatial_dropdown_feat")
                    ui.remove_ui("#inserted-spatial_table")
                    spatial_feature_initialized.set(False)

            elif btn == "Feature":
                if not spatial_feature_initialized.get():
                    dropdown = ui.input_select(
                        "spatial_feat",
                        "Select a Feature",
                        choices=shared['var_names'].get()
                    )
                    ui.insert_ui(
                        ui.div(
                            {"id": "inserted-spatial_dropdown_feat"}, dropdown
                            ),
                        selector="#main-spatial_dropdown_feat",
                        where="beforeEnd"
                    )
                    table_select = ui.input_select(
                        "spatial_layer",
                        "Select a Table",
                        choices=shared['layers_names'].get() + ["Original"],
                        selected="Original"
                    )
                    ui.insert_ui(
                        ui.div({"id": "inserted-spatial_table"}, table_select),
                        selector="#main-spatial_table_dropdown_feat",
                        where="beforeEnd"
                    )
                    spatial_feature_initialized.set(True)

                if spatial_annotation_initialized.get():
                    ui.remove_ui("#inserted-spatial_dropdown_anno")
                    spatial_annotation_initialized.set(False)
