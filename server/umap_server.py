from shiny import ui, render, reactive, req

from utils.template_wrapper import (
    register_memory_object,
    unregister_memory_object,
)
from spac.templates.umap_tsne_pca_template import run_from_json


def _prevent_label_clipping(fig):
    """Add safe margins so axis labels are fully visible."""
    if fig is None:
        return
    try:
        fig.tight_layout(pad=1.4)
    except Exception:
        pass
    try:
        # Keep a small outer frame in case tight_layout is insufficient.
        fig.subplots_adjust(left=0.12, right=0.98, bottom=0.12, top=0.94)
    except Exception:
        pass


def _set_feature_colorbar_fontsize(fig, fontsize):
    """Apply feature-mode legend font size to matplotlib colorbar ticks."""
    if fig is None:
        return
    try:
        size = float(fontsize)
    except Exception:
        return
    for ax in getattr(fig, "axes", []):
        # Colorbar axes are typically very narrow (vertical bar) or short (horizontal bar).
        bbox = ax.get_position()
        is_colorbar_axis = (bbox.width < 0.08) or (bbox.height < 0.08)
        if not is_colorbar_axis:
            continue
        try:
            ax.tick_params(axis="both", labelsize=size)
        except Exception:
            continue


def umap_server(input, output, session, shared):
    @reactive.calc
    def get_adata():
        return shared["adata_main"].get()

    def _figure_legend_params(panel: int):
        """panel 1 or 2. Use pre-template method-like display defaults when collapsed."""
        if panel == 1:
            if not input.umap_show_figure_config():
                return {
                    "Figure_Width": 10,
                    "Figure_Height": 8,
                    "Figure_DPI": 110,
                    "Font_Size": 8,
                    "Legend_Location": "upper right",
                    "Legend_Font_Size": 8,
                    "Legend_Marker_Size": 1.5,
                }
            return {
                "Figure_Width": input.umap_fig_width(),
                "Figure_Height": input.umap_fig_height(),
                "Figure_DPI": input.umap_fig_dpi(),
                "Font_Size": input.umap_font_size(),
                "Legend_Location": "upper right",
                "Legend_Font_Size": input.umap_legend_font_size(),
                "Legend_Marker_Size": input.umap_legend_marker_scale(),
            }
        if not input.umap_show_figure_config2():
            return {
                "Figure_Width": 10,
                "Figure_Height": 8,
                "Figure_DPI": 110,
                "Font_Size": 8,
                "Legend_Location": "upper right",
                "Legend_Font_Size": 8,
                "Legend_Marker_Size": 1.5,
            }
        return {
            "Figure_Width": input.umap_fig_width2(),
            "Figure_Height": input.umap_fig_height2(),
            "Figure_DPI": input.umap_fig_dpi2(),
            "Font_Size": input.umap_font_size2(),
            "Legend_Location": "upper right",
            "Legend_Font_Size": input.umap_legend_font_size2(),
            "Legend_Marker_Size": input.umap_legend_marker_scale2(),
        }

    def _value_range_params(panel: int, mode: str):
        if mode != "Feature":
            return "None", "None"
        if panel == 1:
            if not input.umap_use_val_range():
                return "None", "None"
            return str(input.umap_val_min()), str(input.umap_val_max())
        if not input.umap_use_val_range2():
            return "None", "None"
        return str(input.umap_val_min2()), str(input.umap_val_max2())

    def _run_umap_plot(*, plottype_id: str, rb_id: str, slider_id: str, panel: int):
        adata = get_adata()
        if adata is None:
            return None

        method = getattr(input, plottype_id)()
        point_size = getattr(input, slider_id)()
        mode = getattr(input, rb_id)()

        if mode == "Feature":
            if panel == 1:
                req(input.umap_rb_feat())
                feature = input.umap_rb_feat()
                layer_raw = input.umap_layer()
            else:
                req(input.umap_rb_feat2())
                feature = input.umap_rb_feat2()
                layer_raw = input.umap_layer2()
            layer_str = "Original" if layer_raw == "Original" else layer_raw
            color_by = "Feature"
            annotation_hl = "None"
            feature_hl = feature
        elif mode == "Annotation":
            if panel == 1:
                req(input.umap_rb_anno())
                annotation = input.umap_rb_anno()
            else:
                req(input.umap_rb_anno2())
                annotation = input.umap_rb_anno2()
            color_by = "Annotation"
            annotation_hl = annotation
            feature_hl = "None"
            layer_str = "Original"
        else:
            return None

        fl = _figure_legend_params(panel)
        v_min, v_max = _value_range_params(panel, mode)

        params = {
            "Upstream_Analysis": None,
            "Color_By": color_by,
            "Annotation_to_Highlight": annotation_hl,
            "Feature_to_Highlight": feature_hl,
            "Table": layer_str,
            "Dimension_Reduction_Method": method,
            "Dot_Size": point_size,
            "Value_Min": v_min,
            "Value_Max": v_max,
            **fl,
        }

        virtual_path = None
        try:
            virtual_path = register_memory_object(adata)
            params["Upstream_Analysis"] = virtual_path
            fig, _df = run_from_json(
                json_path=params,
                save_results=False,
                show_plot=False,
            )
            if fig is None:
                return None
            if mode == "Feature":
                _set_feature_colorbar_fontsize(fig, fl["Legend_Font_Size"])
            _prevent_label_clipping(fig)
            return fig
        except Exception:
            import traceback

            traceback.print_exc()
            return None
        finally:
            if virtual_path:
                unregister_memory_object(virtual_path)

    @output
    @render.plot
    @reactive.event(input.go_umap1, ignore_none=True)
    def spac_UMAP():
        return _run_umap_plot(
            plottype_id="plottype",
            rb_id="umap_rb",
            slider_id="umap_slider_1",
            panel=1,
        )

    umap_annotation_initialized = reactive.Value(False)
    umap_feature_initialized = reactive.Value(False)

    @reactive.effect
    def umap_reactivity():
        flipper = shared["data_loaded"].get()
        if flipper is not False:
            btn = input.umap_rb()

            if btn == "Annotation":
                if not umap_annotation_initialized.get():
                    dropdown = ui.input_select(
                        "umap_rb_anno",
                        "Select an Annotation",
                        choices=shared["obs_names"].get(),
                    )
                    ui.insert_ui(
                        ui.div({"id": "inserted-rbdropdown_anno"}, dropdown),
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
                    dropdown1 = ui.input_select(
                        "umap_rb_feat",
                        "Select a Feature",
                        choices=shared["var_names"].get(),
                    )
                    ui.insert_ui(
                        ui.div({"id": "inserted-rbdropdown_feat"}, dropdown1),
                        selector="#main-ump_rb_dropdown_feat",
                        where="beforeEnd",
                    )

                    new_choices = shared["layers_names"].get() + ["Original"]
                    table_umap = ui.input_select(
                        "umap_layer",
                        "Select a Table",
                        choices=new_choices,
                        selected=["Original"],
                    )
                    ui.insert_ui(
                        ui.div({"id": "inserted-umap_table"}, table_umap),
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

    @output
    @render.plot
    @reactive.event(input.go_umap2, ignore_none=True)
    def spac_UMAP2():
        return _run_umap_plot(
            plottype_id="plottype2",
            rb_id="umap_rb2",
            slider_id="umap_slider_2",
            panel=2,
        )

    umap2_annotation_initialized = reactive.Value(False)
    umap2_feature_initialized = reactive.Value(False)

    @reactive.effect
    def umap_reactivity2():
        flipper = shared["data_loaded"].get()
        if flipper is not False:
            btn = input.umap_rb2()

            if btn == "Annotation":
                if not umap2_annotation_initialized.get():
                    dropdown = ui.input_select(
                        "umap_rb_anno2",
                        "Select an Annotation",
                        choices=shared["obs_names"].get(),
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
                        choices=shared["var_names"].get(),
                    )
                    ui.insert_ui(
                        ui.div({"id": "inserted-rbdropdown_feat2"}, dropdown1),
                        selector="#main-ump_rb_dropdown_feat2",
                        where="beforeEnd",
                    )

                    new_choices = shared["layers_names"].get() + ["Original"]
                    table_umap_1 = ui.input_select(
                        "umap_layer2",
                        "Select a Table",
                        choices=new_choices,
                        selected=["Original"],
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
