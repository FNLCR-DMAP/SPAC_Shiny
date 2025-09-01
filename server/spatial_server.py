from shiny import ui, reactive
from shinywidgets import output_widget, render_widget
import anndata as ad
import pandas as pd
import spac.visualization
# Added...
import matplotlib.pyplot as plt


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
        adata = ad.AnnData(
            X=shared['X_data'].get(),
            var=pd.DataFrame(shared['var_data'].get()),
            obsm=shared['obsm_data'].get(),
            obs=shared['obs_data'].get(),
            dtype=shared['X_data'].get().dtype,
            layers=shared['layers_data'].get()
        )
        slide_check = input.slide_select_check()
        region_check = input.region_select_check()
        # Added...
        font_size = input.spatial_font_size()

        # Added...
        plt.rcParams.update({'font.size': font_size})

        if adata is not None:
            if slide_check is False and region_check is False:
                adata_subset = adata
            elif slide_check is True and region_check is False:
                adata_subset = adata[
                    adata.obs[
                        input.slide_select_drop()
                    ] == input.slide_select_label()
                ].copy()
            elif slide_check is True and region_check is True:
                adata_subset = adata[
                    (adata.obs[
                        input.slide_select_drop()
                    ] == input.slide_select_label()) &
                    (adata.obs[
                        input.region_select_drop()
                    ] == input.region_label_select())
                ].copy()
            elif slide_check is False and region_check is True:
                adata_subset = adata[
                    adata.obs[
                        input.region_select_drop()
                    ] == input.region_label_select()
                ].copy()
            else:
                return None
            if input.spatial_rb() == "Feature":
                if "spatial_feat" not in input or input.spatial_feat() is None:
                    return None
                layer = (
                    None if input.spatial_layer() == "Original"
                    else input.spatial_layer()
                )
                out = spac.visualization.interactive_spatial_plot(
                    adata_subset,
                    feature=input.spatial_feat(),
                    layer=layer,
                    figure_width=5.5,
                    figure_height=5,
                    dot_size=input.spatial_slider()
                )
            elif input.spatial_rb() == "Annotation":
                if "spatial_anno" not in input or input.spatial_anno() is None:
                    return None
                out = spac.visualization.interactive_spatial_plot(
                    adata_subset,
                    annotations=input.spatial_anno(),
                    figure_width=5.5,
                    figure_height=5,
                    dot_size=input.spatial_slider()
                )
            else:
                return None
            out[0]['image_object'].update_xaxes(
                showticklabels=True,
                ticks="outside",
                tickwidth=2,
                ticklen=10
            )
            out[0]['image_object'].update_yaxes(
                showticklabels=True,
                ticks="outside",
                tickwidth=2,
                ticklen=10
            )
            # Added...
            out[0]['image_object'].update_layout(font=dict(size=font_size))
            return out[0]['image_object']

        return None

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
                        "spatial_anno", "Select an Annotation", choices=shared['obs_names'].get()
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