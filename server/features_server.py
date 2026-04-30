from shiny import ui, render, reactive
import anndata as ad
import numpy as np
import pandas as pd
import spac.visualization


def features_server(input, output, session, shared):
    # def on_layer_check():
    #     return input.h1_layer() if input.h1_layer() != "Original" else None

    def get_bins_value():
        bins_type = input.h1_bins_type()

        if bins_type == "auto":
            return None

        elif bins_type == "number":
            val = input.h1_bins_number()
            return int(val) if val is not None else None

        elif bins_type == "list":
            raw = input.h1_bins_list()
            if not raw or not raw.strip():
                return None
            try:
                return [float(x.strip()) for x in raw.split(",") if x.strip()]
            except ValueError:
                return None
    @reactive.calc
    def get_layer():
        """
        Return None for 'Original', otherwise return selected layer.

        Returns
        -------
        str or None
            Selected layer name or None for original data
        """
        layer = input.h1_layer()
        return None if layer == "Original" else layer
    
    @reactive.calc
    def get_group_by():
        """
        Return group_by annotation if Group By is checked, else None.

        Returns
        -------
        str or None
            Annotation column name or None
        """
        if input.h1_group_by_check():
            return input.h1_anno()
        return None
    
    @reactive.calc
    def get_element():
        """
        Return element value, defaulting to 'bars' if hidden/empty.

        Returns
        -------
        str
            Element type: 'bars', 'step', or 'poly'
        """
        if not input.h1_show_element():
            return "bars"  # default when section is hidden
        val = input.h1_element()
        return val if val else "bars"
    
    @reactive.calc
    def get_together():
        if not input.h1_group_by_check():
            return False
        if input.h1_facet():        # ← Facet takes priority
            return False
        return input.h1_together_check()

    @reactive.calc
    def get_multiple():
        """
        Return stack type if Plot Together is checked, else None.

        Returns
        -------
        str or None
            Stack type string or None
        """
        if input.h1_group_by_check() and input.h1_together_check():
            return input.h1_together_drop()
        return None
    
        # ▼▼▼ NEW: Facet ↔ Plot Together mutual exclusion
    @reactive.effect
    @reactive.event(input.h1_facet)
    def on_facet_changed():
        if input.h1_facet():
            try:
                ui.update_checkbox("h1_together_check", value=False)
            except Exception:
                pass

    @reactive.effect
    @reactive.event(input.h1_together_check)
    def on_together_changed():
        try:
            if input.h1_together_check():
                ui.update_checkbox("h1_facet", value=False)
        except Exception:
            pass
    # ▲▲▲ END NEW

    @reactive.effect
    @reactive.event(input.feat_slider)
    def sync_slider_to_num():
        ui.update_numeric("feat_slider_num", value=input.feat_slider())

    @reactive.effect
    @reactive.event(input.feat_slider_num)
    def sync_num_to_slider():
        val = input.feat_slider_num()
        if val is not None and 0 <= val <= 90:
            ui.update_slider("feat_slider", value=val)

    @output
    @render.plot
    @reactive.event(input.go_h1, ignore_none=True)
    def spac_Histogram_1():
        adata = ad.AnnData(
            X=shared['X_data'].get(),
            obs=pd.DataFrame(shared['obs_data'].get()),
            var=pd.DataFrame(shared['var_data'].get()),
            layers=shared['layers_data'].get(),
            dtype=shared['X_data'].get().dtype
        )

        if adata is None:
            return None

        # feature = input.h1_feat()
        # rotation = input.feat_slider()
        # btn_log_x = input.h1_log_x()
        # btn_log_y = input.h1_log_y()
        # layer = on_layer_check()
        # stat = input.h1_stat()
        # element = input.h1_element()

        # call get_bins_value() and add to kwargs
        bins_val = get_bins_value()

        # kwargs = {
        #     "adata": adata,
        #     "feature": feature,
        #     "layer": layer,
        #     "x_log_scale": btn_log_x,
        #     "y_log_scale": btn_log_y,
        #     "element": element,
        #     "stat": stat,
        # }
        # Register adata in memory and get virtual path
        
        from utils.template_wrapper import (
                register_memory_object,
                unregister_memory_object
            )
        from spac.templates.histogram_template import run_from_json

        virtual_path = register_memory_object(adata)
        params = {
                        "Upstream_Analysis": virtual_path,
                        "Feature": input.h1_feat(),
                        "Layer": get_layer() or "None",
                        "X_Log_Scale": input.h1_log_x(),
                        "Y_Log_Scale": input.h1_log_y(),
                        "Element": input.h1_element(),
                        "Stat": input.h1_stat(),
                        "Bins": get_bins_value() if get_bins_value() is not None else "auto",
                        "Group_by": get_group_by() or "None",
                        "Together": (
                            input.h1_together_check()
                            if input.h1_group_by_check()
                            else False
                        ),
                        "Multiple": get_multiple() or "dodge",
                        "X_Axis_Label_Rotation": input.feat_slider(),
                        "Figure_Width": input.h1_figure_width(),
                        "Figure_Height": input.h1_figure_height(),
                        "Figure_DPI": input.h1_figure_dpi(),
                        "Font_Size": input.h1_font_size(),
                        "Plot_By": "Feature",
                        "Facet": input.h1_facet(),
                    }
        
        # only pass bins if not None, avoids overriding auto behaviour ▼▼▼
        if bins_val is not None:
            params["Bins"] = bins_val

        if input.h1_group_by_check():
            params["Group_by"] = input.h1_anno()
            params["Together"] = input.h1_together_check()
            if input.h1_together_check():
                params["Multiple"] = input.h1_together_drop()
        try:
                # Call run_from_json with virtual path
            figs, df_data = run_from_json(
                json_path=params,
                save_to_disk=False,
                show_plot=False
            )
        finally:
                    # Always clean up memory registry
            unregister_memory_object(virtual_path)
    
        # fig1, ax, df = spac.visualization.histogram(**params).values()

        # axes = ax if isinstance(ax, (list, np.ndarray)) else [ax]
        # for a in axes:
        #     a.tick_params(axis='x', rotation=params["X_Axis_Label_Rotation"], labelsize=10)

        shared['df_histogram1'].set(df_data)
        return figs

    histogram_ui_initialized = reactive.Value(False)

    @render.download(filename="features_histogram_data.csv")
    def download_histogram1_df():
        df = shared['df_histogram1'].get()
        if df is not None:
            csv_string = df.to_csv(index=False)
            csv_bytes = csv_string.encode("utf-8")
            return csv_bytes, "text/csv"
        return None

    @render.ui
    @reactive.event(input.go_h1, ignore_none=True)
    def download_histogram1_button_ui():
        if shared['df_histogram1'].get() is not None:
            return ui.download_button(
                "download_histogram1_df",
                "Download Data",
                class_="btn-warning"
            )
        return None

    @reactive.effect
    def histogram_reactivity():
        btn = input.h1_group_by_check()
        ui_initialized = histogram_ui_initialized.get()

        if btn and not ui_initialized:
            dropdown = ui.input_select(
                "h1_anno",
                "Select an Annotation",
                choices=shared['obs_names'].get()
            )
            ui.insert_ui(
                ui.div({"id": "inserted-dropdown"}, dropdown),
                selector="#main-h1_dropdown",
                where="beforeEnd",
            )

            together_check = ui.input_checkbox(
                "h1_together_check",
                "Plot Together",
                value=True
            )
            ui.insert_ui(
                ui.div({"id": "inserted-check"}, together_check),
                selector="#main-h1_check",
                where="beforeEnd",
            )
            histogram_ui_initialized.set(True)

        elif not btn and ui_initialized:
            ui.remove_ui("#inserted-dropdown")
            ui.remove_ui("#inserted-check")
            ui.remove_ui("#inserted-dropdown_together")
            histogram_ui_initialized.set(False)

    @reactive.effect
    @reactive.event(input.h1_together_check)
    def update_stack_type_dropdown():
        if input.h1_together_check():
            dropdown_together = ui.input_select(
                "h1_together_drop",
                "Select Stack Type",
                choices=['stack', 'layer', 'dodge', 'fill'],
                selected='stack'
            )
            ui.insert_ui(
                ui.div(
                    {"id": "inserted-dropdown_together"},
                    dropdown_together
                ),
                selector="#main-h1_together_drop",
                where="beforeEnd",
            )
        else:
            ui.remove_ui("#inserted-dropdown_together")