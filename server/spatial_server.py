from shiny import ui, render, reactive
from shinywidgets import render_widget
import anndata as ad
import pandas as pd
import spac.visualization
from utils.plot_manager import PlotManager


def spatial_server(input, output, session, shared):
    pm = PlotManager('spatial', shared, plot_type='thread', data_key=None)

    # ── Slide UI ──────────────────────────────────────────────────────────
    slide_ui_initialized = reactive.Value(False)

    @reactive.effect
    def slide_reactivity():
        btn = input.slide_select_check()
        ui_initialized = slide_ui_initialized.get()

        if btn and not ui_initialized:
            ui.insert_ui(
                ui.div(
                    {"id": "inserted-slide_dropdown"},
                    ui.input_select(
                        "slide_select_drop",
                        "Select the Slide Annotation",
                        choices=shared['obs_names'].get()
                    )
                ),
                selector="#main-slide_dropdown",
                where="beforeEnd",
            )
            ui.insert_ui(
                ui.div(
                    {"id": "inserted-label_dropdown"},
                    ui.input_select(
                        "slide_select_label",
                        "Select a Slide",
                        choices=[]
                    )
                ),
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
            labels = adata.obs[input.slide_select_drop()].unique().tolist()
            ui.update_select("slide_select_label", choices=labels)

    # ── Region UI ─────────────────────────────────────────────────────────
    region_ui_initialized = reactive.Value(False)

    @reactive.effect
    def region_reactivity():
        btn = input.region_select_check()
        ui_initialized = region_ui_initialized.get()

        if btn and not ui_initialized:
            ui.insert_ui(
                ui.div(
                    {"id": "inserted-region_dropdown"},
                    ui.input_select(
                        "region_select_drop",
                        "Select the Region Annotation",
                        choices=shared['obs_names'].get()
                    )
                ),
                selector="#main-region_dropdown",
                where="beforeEnd",
            )
            ui.insert_ui(
                ui.div(
                    {"id": "inserted-region_label_select_dropdown"},
                    ui.input_select(
                        "region_label_select",
                        "Select a Region",
                        choices=[]
                    )
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
            labels = adata.obs[input.region_select_drop()].unique().tolist()
            ui.update_select("region_label_select", choices=labels)

    # ── Spatial annotation/feature UI ────────────────────────────────────
    spatial_annotation_initialized = reactive.Value(False)
    spatial_feature_initialized = reactive.Value(False)

    @reactive.effect
    def spatial_reactivity():
        flipper = shared['data_loaded'].get()
        if flipper is not False:
            btn = input.spatial_rb()

            if btn == "Annotation":
                if not spatial_annotation_initialized.get():
                    ui.insert_ui(
                        ui.div(
                            {"id": "inserted-spatial_dropdown_anno"},
                            ui.input_select(
                                "spatial_anno",
                                "Select an Annotation",
                                choices=shared['obs_names'].get()
                            )
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
                    ui.insert_ui(
                        ui.div(
                            {"id": "inserted-spatial_dropdown_feat"},
                            ui.input_select(
                                "spatial_feat",
                                "Select a Feature",
                                choices=shared['var_names'].get()
                            )
                        ),
                        selector="#main-spatial_dropdown_feat",
                        where="beforeEnd"
                    )
                    ui.insert_ui(
                        ui.div(
                            {"id": "inserted-spatial_table"},
                            ui.input_select(
                                "spatial_layer",
                                "Select a Table",
                                choices=shared['layers_names'].get() + ["Original"],
                                selected="Original"
                            )
                        ),
                        selector="#main-spatial_table_dropdown_feat",
                        where="beforeEnd"
                    )
                    spatial_feature_initialized.set(True)
                if spatial_annotation_initialized.get():
                    ui.remove_ui("#inserted-spatial_dropdown_anno")
                    spatial_annotation_initialized.set(False)

    # ── Plot ──────────────────────────────────────────────────────────────
    @reactive.Effect
    @reactive.event(input.go_sp1, ignore_none=True)
    def start_spatial_task():
        if pm.is_calculating.get():
            return

        # CAPTURE ALL INPUTS BEFORE THREAD STARTS
        adata = ad.AnnData(
            X=shared['X_data'].get(),
            var=pd.DataFrame(shared['var_data'].get()),
            obsm=shared['obsm_data'].get(),
            obs=shared['obs_data'].get(),
            dtype=shared['X_data'].get().dtype,
            layers=shared['layers_data'].get()
        )
        if adata is None:
            return

        slide_check = input.slide_select_check()
        region_check = input.region_select_check()
        spatial_mode = input.spatial_rb()
        dot_size = input.spatial_slider()

        slide_drop = input.slide_select_drop() if slide_check else None
        slide_label = input.slide_select_label() if slide_check else None
        region_drop = input.region_select_drop() if region_check else None
        region_label = input.region_label_select() if region_check else None

        spatial_feat = input.spatial_feat() if spatial_mode == "Feature" else None
        spatial_anno = input.spatial_anno() if spatial_mode == "Annotation" else None
        spatial_layer = input.spatial_layer() if spatial_mode == "Feature" else None

        def worker():
            try:
                # Subset adata
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
                    pm.result_queue.put("Error: Invalid subset combination")
                    return

                if spatial_mode == "Feature":
                    if spatial_feat is None:
                        pm.result_queue.put("Error: No feature selected")
                        return
                    layer = None if spatial_layer == "Original" else spatial_layer
                    out = spac.visualization.interactive_spatial_plot(
                        adata_subset,
                        feature=spatial_feat,
                        layer=layer,
                        figure_width=5.5,
                        figure_height=5,
                        dot_size=dot_size
                    )
                elif spatial_mode == "Annotation":
                    if spatial_anno is None:
                        pm.result_queue.put("Error: No annotation selected")
                        return
                    out = spac.visualization.interactive_spatial_plot(
                        adata_subset,
                        annotations=spatial_anno,
                        figure_width=5.5,
                        figure_height=5,
                        dot_size=dot_size
                    )
                else:
                    pm.result_queue.put("Error: Invalid spatial mode")
                    return

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
                pm.result_queue.put(fig)

            except Exception as e:
                pm.result_queue.put(f"Error: {str(e)}")

        pm.start_thread(worker)

    @reactive.Effect
    def check_status():
        def on_result(res):
            pm.result.set(res)
        pm.check_thread(on_result)

    @render_widget
    def spac_Spatial():
        return pm.result.get()

    @render.ui
    def spatial_stop_button_ui():
        return pm.stop_button_ui('stop_spatial')