from shiny import ui, render, reactive
import anndata as ad
import pandas as pd
import spac.visualization
from utils.plot_utils import fig_to_png_bytes, png_bytes_to_figure


def umap_server(input, output, session, shared):
    @output
    @render.plot
    @reactive.event(input.go_umap1, ignore_none=True)
    def spac_UMAP():
        x_data = shared['X_data'].get()
        if x_data is None:
            return None

        cache = shared['cache']
        version = shared['dataset_version'].get()

        mode = input.umap_rb()
        method = input.plottype()
        point_size = input.umap_slider_1()

        try:
            feature = input.umap_rb_feat() if mode == "Feature" else None
        except Exception:
            feature = None

        try:
            layer = (
                None if input.umap_layer() == "Original"
                else input.umap_layer()
            ) if mode == "Feature" else None
        except Exception:
            layer = None

        try:
            annotation = input.umap_rb_anno() if mode == "Annotation" else None
        except Exception:
            annotation = None

        params = {
            'mode': mode,
            'method': method,
            'point_size': point_size,
            'feature': feature,
            'layer': layer,
            'annotation': annotation,
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

            if mode == "Feature" and feature:
                fig, ax = spac.visualization.dimensionality_reduction_plot(
                    adata,
                    method=method,
                    feature=feature,
                    layer=layer,
                    point_size=point_size
                )
                ax.set_title(f"{method.upper()}: {feature}", fontsize=14)
                ax.set_xlabel(f"{method.upper()} 1")
                ax.set_ylabel(f"{method.upper()} 2")
                for extra_ax in fig.axes:
                    if hasattr(extra_ax, "get_ylabel") and extra_ax != ax:
                        extra_ax.set_ylabel(
                            f"Colored by: {feature.upper()}", fontsize=12
                        )
                return fig_to_png_bytes(fig), None

            elif mode == "Annotation" and annotation:
                fig, ax = spac.visualization.dimensionality_reduction_plot(
                    adata,
                    method=method,
                    annotation=annotation,
                    point_size=point_size
                )
                ax.set_title(f"{method.upper()}: {annotation}", fontsize=14)
                ax.set_xlabel(f"{method.upper()} 1")
                ax.set_ylabel(f"{method.upper()} 2")
                return fig_to_png_bytes(fig), None

            return None, None

        img_bytes, _ = cache.get_or_compute('umap1', version, params, compute)
        if img_bytes is None:
            return None
        return png_bytes_to_figure(img_bytes)

    # Track the UI state
    umap_annotation_initialized = reactive.Value(False)
    umap_feature_initialized = reactive.Value(False)

    @reactive.effect
    def umap_reactivity():
        flipper = shared['data_loaded'].get()
        if flipper is not False:
            btn = input.umap_rb()

            if btn == "Annotation":
                if not umap_annotation_initialized.get():
                    # Create the Annotation dropdown
                    dropdown = ui.input_select(
                        "umap_rb_anno",
                        "Select an Annotation",
                        choices=shared['obs_names'].get(),
                    )
                    ui.insert_ui(
                        ui.div({"id": "inserted-rbdropdown_anno"}, dropdown),
                        selector="#main-ump_rb_dropdown_anno",
                        where="beforeEnd",
                    )
                    # Update the state
                    umap_annotation_initialized.set(True)
                # Remove the Feature dropdown and table
                if umap_feature_initialized.get():
                    ui.remove_ui("#inserted-rbdropdown_feat")
                    ui.remove_ui("#inserted-umap_table")
                    umap_feature_initialized.set(False)

            elif btn == "Feature":
                if not umap_feature_initialized.get():
                    # Create the Feature dropdown
                    dropdown1 = ui.input_select(
                        "umap_rb_feat",
                        "Select a Feature",
                        choices=shared['var_names'].get()
                    )
                    ui.insert_ui(
                        ui.div({"id": "inserted-rbdropdown_feat"}, dropdown1),
                        selector="#main-ump_rb_dropdown_feat",
                        where="beforeEnd",
                    )

                    # Create the Table dropdown
                    new_choices = shared['layers_names'].get() + ["Original"]
                    table_umap = ui.input_select(
                        "umap_layer",
                        "Select a Table",
                        choices=new_choices,
                        selected=["Original"]
                    )
                    ui.insert_ui(
                        ui.div({"id": "inserted-umap_table"}, table_umap),
                        selector="#main-ump_table_dropdown_feat",
                        where="beforeEnd",
                    )
                    # Update the state
                    umap_feature_initialized.set(True)
                # Remove the Annotation dropdown
                if umap_annotation_initialized.get():
                    ui.remove_ui("#inserted-rbdropdown_anno")
                    umap_annotation_initialized.set(False)

            elif btn == "None":
                # Remove all dropdowns and reset states
                if umap_annotation_initialized.get():
                    ui.remove_ui("#inserted-rbdropdown_anno")
                    umap_annotation_initialized.set(False)
                if umap_feature_initialized.get():
                    ui.remove_ui("#inserted-rbdropdown_feat")
                    ui.remove_ui("#inserted-umap_table")
                    umap_feature_initialized.set(False)


    @output
    @render.plot
    @reactive.event(input.go_umap2, ignore_none=True)
    def spac_UMAP2():
        x_data = shared['X_data'].get()
        if x_data is None:
            return None

        cache = shared['cache']
        version = shared['dataset_version'].get()

        mode = input.umap_rb2()
        method = input.plottype2()
        point_size = input.umap_slider_2()

        try:
            feature = input.umap_rb_feat2() if mode == "Feature" else None
        except Exception:
            feature = None

        try:
            layer = (
                None if input.umap_layer2() == "Original"
                else input.umap_layer2()
            ) if mode == "Feature" else None
        except Exception:
            layer = None

        try:
            annotation = input.umap_rb_anno2() if mode == "Annotation" else None
        except Exception:
            annotation = None

        params = {
            'mode': mode,
            'method': method,
            'point_size': point_size,
            'feature': feature,
            'layer': layer,
            'annotation': annotation,
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

            if mode == "Feature" and feature:
                fig, ax = spac.visualization.dimensionality_reduction_plot(
                    adata,
                    method=method,
                    feature=feature,
                    layer=layer,
                    point_size=point_size
                )
                ax.set_title(f"{method.upper()}: {feature}", fontsize=14)
                ax.set_xlabel(f"{method.upper()} 1")
                ax.set_ylabel(f"{method.upper()} 2")
                for extra_ax in fig.axes:
                    if hasattr(extra_ax, "get_ylabel") and extra_ax != ax:
                        extra_ax.set_ylabel(
                            f"Colored by: {feature}", fontsize=12
                        )
                return fig_to_png_bytes(fig), None

            elif mode == "Annotation" and annotation:
                fig, ax = spac.visualization.dimensionality_reduction_plot(
                    adata,
                    method=method,
                    annotation=annotation,
                    point_size=point_size
                )
                ax.set_title(f"{method.upper()}: {annotation}", fontsize=14)
                ax.set_xlabel(f"{method.upper()} 1")
                ax.set_ylabel(f"{method.upper()} 2")
                return fig_to_png_bytes(fig), None

            return None, None

        img_bytes, _ = cache.get_or_compute('umap2', version, params, compute)
        if img_bytes is None:
            return None
        return png_bytes_to_figure(img_bytes)

    # Track the UI state
    umap2_annotation_initialized = reactive.Value(False)
    umap2_feature_initialized = reactive.Value(False)

    @reactive.effect
    def umap_reactivity2():
        flipper = shared['data_loaded'].get()
        if flipper is not False:
            btn = input.umap_rb2()

            if btn == "Annotation":
                if not umap2_annotation_initialized.get():
                    dropdown = ui.input_select(
                        "umap_rb_anno2",
                        "Select an Annotation",
                        choices=shared['obs_names'].get()
                    )
                    ui.insert_ui(
                        ui.div({"id": "inserted-rbdropdown_anno2"}, dropdown),
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
                    dropdown1 = ui.input_select(
                        "umap_rb_feat2",
                        "Select a Feature",
                        choices=shared['var_names'].get()
                    )
                    ui.insert_ui(
                        ui.div({"id": "inserted-rbdropdown_feat2"}, dropdown1),
                        selector="#main-ump_rb_dropdown_feat2",
                        where="beforeEnd",
                    )

                    new_choices = shared['layers_names'].get() + ["Original"]
                    table_umap_1 = ui.input_select(
                        "umap_layer2",
                        "Select a Table",
                        choices=new_choices,
                        selected=["Original"]
                    )
                    ui.insert_ui(
                        ui.div({"id": "inserted-umap_table2"}, table_umap_1),
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
