from shiny import ui, render, reactive
import anndata as ad
import pandas as pd
import io
import tempfile
import multiprocessing
import matplotlib.pyplot as plt
import spac.visualization
from utils.plot_manager import PlotManager


def run_umap_worker(queue, adata, method, mode, point_size,
                    feature, layer, annotation):
    try:
        plt.clf()
        plt.close('all')

        if mode == "Feature":
            fig, ax = spac.visualization.dimensionality_reduction_plot(
                adata, method=method, feature=feature,
                layer=layer, point_size=point_size
            )
            ax.set_title(f"{method.upper()}: {feature}", fontsize=14)
            ax.set_xlabel(f"{method.upper()} 1")
            ax.set_ylabel(f"{method.upper()} 2")
            for extra_ax in fig.axes:
                if hasattr(extra_ax, "get_ylabel") and extra_ax != ax:
                    extra_ax.set_ylabel(
                        f"Colored by: {feature.upper()}", fontsize=12
                    )

        elif mode == "Annotation":
            fig, ax = spac.visualization.dimensionality_reduction_plot(
                adata, method=method, annotation=annotation,
                point_size=point_size
            )
            ax.set_title(f"{method.upper()}: {annotation}", fontsize=14)
            ax.set_xlabel(f"{method.upper()} 1")
            ax.set_ylabel(f"{method.upper()} 2")
        else:
            queue.put("Error: Invalid mode")
            return

        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=140,
                    bbox_inches='tight', pad_inches=0.2)
        buf.seek(0)
        img_bytes = buf.read()
        plt.close(fig)

        queue.put(img_bytes)

    except Exception as e:
        queue.put(f"Error: {str(e)}")


def umap_server(input, output, session, shared):
    pm1 = PlotManager('umap1', shared, plot_type='process', data_key=None)
    pm2 = PlotManager('umap2', shared, plot_type='process', data_key=None)

    # ── UMAP 1 ────────────────────────────────────────────────────────────
    umap_annotation_initialized = reactive.Value(False)
    umap_feature_initialized = reactive.Value(False)

    @reactive.effect
    def umap_reactivity():
        flipper = shared['data_loaded'].get()
        if flipper is not False:
            btn = input.umap_rb()

            if btn == "Annotation":
                if not umap_annotation_initialized.get():
                    ui.insert_ui(
                        ui.div(
                            {"id": "inserted-rbdropdown_anno"},
                            ui.input_select(
                                "umap_rb_anno",
                                "Select an Annotation",
                                choices=shared['obs_names'].get()
                            )
                        ),
                        selector="#main-ump_rb_dropdown_anno",
                        where="beforeEnd",
                    )
                    umap_annotation_initialized.set(True)
                if umap_feature_initialized.get():
                    ui.remove_ui("#inserted-rbdropdown_feat")
                    ui.remove_ui("#inserted-umap_table")
                    umap_feature_initialized.set(False)

            elif btn == "Feature":
                if not umap_feature_initialized.get():
                    ui.insert_ui(
                        ui.div(
                            {"id": "inserted-rbdropdown_feat"},
                            ui.input_select(
                                "umap_rb_feat",
                                "Select a Feature",
                                choices=shared['var_names'].get()
                            )
                        ),
                        selector="#main-ump_rb_dropdown_feat",
                        where="beforeEnd",
                    )
                    ui.insert_ui(
                        ui.div(
                            {"id": "inserted-umap_table"},
                            ui.input_select(
                                "umap_layer",
                                "Select a Table",
                                choices=shared['layers_names'].get() + ["Original"],
                                selected=["Original"]
                            )
                        ),
                        selector="#main-ump_table_dropdown_feat",
                        where="beforeEnd",
                    )
                    umap_feature_initialized.set(True)
                if umap_annotation_initialized.get():
                    ui.remove_ui("#inserted-rbdropdown_anno")
                    umap_annotation_initialized.set(False)

            elif btn == "None":
                if umap_annotation_initialized.get():
                    ui.remove_ui("#inserted-rbdropdown_anno")
                    umap_annotation_initialized.set(False)
                if umap_feature_initialized.get():
                    ui.remove_ui("#inserted-rbdropdown_feat")
                    ui.remove_ui("#inserted-umap_table")
                    umap_feature_initialized.set(False)

    @reactive.Effect
    @reactive.event(input.go_umap1, ignore_none=True)
    def start_umap1_task():
        if pm1.is_calculating.get():
            return

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

        method = input.plottype()
        point_size = input.umap_slider_1()
        mode = input.umap_rb()
        feature = input.umap_rb_feat() if mode == "Feature" else None
        layer = (None if input.umap_layer() == "Original"
                 else input.umap_layer()) if mode == "Feature" else None
        annotation = input.umap_rb_anno() if mode == "Annotation" else None

        pm1.start_process(
            run_umap_worker,
            args=(adata, method, mode, point_size, feature, layer, annotation)
        )

    @reactive.Effect
    def check_umap1_status():
        def on_result(res):
            pm1.result.set(res)
        pm1.check_process(on_result)

    @output
    @render.image
    def spac_UMAP():
        img_bytes = pm1.result.get()
        if img_bytes is None:
            return None
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            tmp.write(img_bytes)
            tmp_path = tmp.name
        return {
            "src": tmp_path,
            "contentType": "image/png",
            "style": "max-width: 100%; height: auto;"
        }

    @render.ui
    def umap1_stop_button_ui():
        return pm1.stop_button_ui('stop_umap1')

    # ── UMAP 2 ────────────────────────────────────────────────────────────
    umap2_annotation_initialized = reactive.Value(False)
    umap2_feature_initialized = reactive.Value(False)

    @reactive.effect
    def umap_reactivity2():
        flipper = shared['data_loaded'].get()
        if flipper is not False:
            btn = input.umap_rb2()

            if btn == "Annotation":
                if not umap2_annotation_initialized.get():
                    ui.insert_ui(
                        ui.div(
                            {"id": "inserted-rbdropdown_anno2"},
                            ui.input_select(
                                "umap_rb_anno2",
                                "Select an Annotation",
                                choices=shared['obs_names'].get()
                            )
                        ),
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
                    ui.insert_ui(
                        ui.div(
                            {"id": "inserted-rbdropdown_feat2"},
                            ui.input_select(
                                "umap_rb_feat2",
                                "Select a Feature",
                                choices=shared['var_names'].get()
                            )
                        ),
                        selector="#main-ump_rb_dropdown_feat2",
                        where="beforeEnd",
                    )
                    ui.insert_ui(
                        ui.div(
                            {"id": "inserted-umap_table2"},
                            ui.input_select(
                                "umap_layer2",
                                "Select a Table",
                                choices=shared['layers_names'].get() + ["Original"],
                                selected=["Original"]
                            )
                        ),
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

    @reactive.Effect
    @reactive.event(input.go_umap2, ignore_none=True)
    def start_umap2_task():
        if pm2.is_calculating.get():
            return

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

        method = input.plottype2()
        point_size = input.umap_slider_2()
        mode = input.umap_rb2()
        feature = input.umap_rb_feat2() if mode == "Feature" else None
        layer = (None if input.umap_layer2() == "Original"
                 else input.umap_layer2()) if mode == "Feature" else None
        annotation = input.umap_rb_anno2() if mode == "Annotation" else None

        pm2.start_process(
            run_umap_worker,
            args=(adata, method, mode, point_size, feature, layer, annotation)
        )

    @reactive.Effect
    def check_umap2_status():
        def on_result(res):
            pm2.result.set(res)
        pm2.check_process(on_result)

    @output
    @render.image
    def spac_UMAP2():
        img_bytes = pm2.result.get()
        if img_bytes is None:
            return None
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            tmp.write(img_bytes)
            tmp_path = tmp.name
        return {
            "src": tmp_path,
            "contentType": "image/png",
            "style": "max-width: 100%; height: auto;"
        }

    @render.ui
    def umap2_stop_button_ui():
        return pm2.stop_button_ui('stop_umap2')

    @render.ui
    def download_umap1_plot_button_ui():
        return pm1.plot_download_button_ui('download_umap1_plot')

    @render.download(filename="umap1_plot.png")
    def download_umap1_plot():
        return pm1.create_plot_download_handler()()

    @render.ui
    def download_umap2_plot_button_ui():
        return pm2.plot_download_button_ui('download_umap2_plot')

    @render.download(filename="umap2_plot.png")
    def download_umap2_plot():
        return pm2.create_plot_download_handler()()