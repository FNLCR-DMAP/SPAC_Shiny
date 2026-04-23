from shiny import App, ui, reactive, render
import os
import json
import re
import httpx
from rag import retrieve_context


# Screen imports
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

# Server imports
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

# Util imports
from utils.data_processing import load_data, read_html_file
from utils.accessibility import accessible_navigation, apply_slider_accessibility_global
from utils.security import apply_security_enhancements




DATA_PATH       = os.getenv("DATA_PATH", "dev_example.pickle")
OLLAMA_URL      = os.getenv("OLLAMA_URL", "http://host.docker.internal:11434")
OLLAMA_MODEL    = os.getenv("OLLAMA_MODEL", "gemma3:latest")
RAG_N_RESULTS   = int(os.getenv("RAG_N_RESULTS", "3"))
OLLAMA_TIMEOUT  = int(os.getenv("OLLAMA_TIMEOUT", "120"))


# ---------------------------------------------------------------------------
# Plot JSON schema
# ---------------------------------------------------------------------------
# Maps each plot type to the exact Shiny input IDs the server expects.
# When the LLM generates a plot spec it must use these field names so that
# the "Apply to App" button can call ui.update_* on the correct inputs.

PLOT_SCHEMA = {
    "boxplot": {
        "description": "Compare continuous feature distributions across annotation groups.",
        "tab_id": "Boxplot",
        # Field names here are the LLM-facing semantic names.
        # They are mapped to actual Shiny input IDs in _apply_plot_params.
        "fields": {
            "features":     "One or more features/proteins to plot (list of strings from adata.var_names). Can be a single item list.",
            "annotation":   "Categorical column from adata.obs to group cells by (string). Use 'No Annotation' to skip grouping.",
            "layer":        "Expression table/layer to use (string). Use 'Original' for the default adata.X matrix.",
            "outliers":     "How to show outliers: 'all' to show all, 'downsample' for a subset, 'none' to hide (default 'none').",
            "log_scale":    "Whether to log-transform the axis (boolean, default false).",
            "horizontal":   "Whether to flip to a horizontal orientation (boolean, default false).",
            "interactive":  "Whether to render an interactive Plotly plot vs a static image (boolean, default true)."
        }
    },
    "feat_vs_anno": {
        "description": "Violin / bar plot of a feature split by an annotation category.",
        "tab_id": "Feature vs. Annotation",
        "fields": {
            "feature":    "Continuous feature to visualise (string).",
            "annotation": "Categorical annotation column (string)."
        }
    },
    "anno_vs_anno": {
        "description": "Visualise relationships between two categorical annotations via a Sankey diagram or a Relational Heatmap.",
        "tab_id": "Anno. Vs Anno.",
        "fields": {
            "source_annotation": "Source (left / row) categorical annotation column (string).",
            "target_annotation": "Target (right / column) categorical annotation column (string).",
            "plot_subtype":      "Which visualisation to generate: 'sankey' (default) or 'relational_heatmap'."
        }
    },
    "spatial": {
        "description": "Cells plotted in physical tissue coordinates coloured by a feature or annotation.",
        "tab_id": "Spatial",
        "fields": {
            "color_by":    "Feature name or annotation column to colour cells (string).",
            "point_size":  "Marker size in pixels (number, default 3).",
            "alpha":       "Marker opacity 0-1 (number, default 0.8).",
            "spatial_key": "Key in adata.obsm containing 2-D spatial coordinates (string, default 'spatial')."
        }
    },
    "umap": {
        "description": "UMAP dimensionality-reduction embedding coloured by a feature or annotation.",
        "tab_id": "UMAP",
        "fields": {
            "color_by":   "Feature name or annotation column to colour cells (string).",
            "umap_key":   "Key in adata.obsm for UMAP coordinates (string, default 'X_umap').",
            "point_size": "Marker size in pixels (number, default 3)."
        }
    },
    "scatterplot": {
        "description": "Scatter plot of two continuous features against each other.",
        "tab_id": "Scatterplot",
        "fields": {
            "feature_x":  "Feature for the x-axis (string).",
            "feature_y":  "Feature for the y-axis (string).",
            "color_by":   "Optional annotation or feature to colour points (string or null).",
            "log_x":      "Log-transform x-axis (boolean, default false).",
            "log_y":      "Log-transform y-axis (boolean, default false)."
        }
    },
    "nearest_neighbor": {
        "description": "Heatmap of how often cell types appear as spatial nearest neighbours.",
        "tab_id": "Nearest Neighbor",
        "fields": {
            "annotation": "Categorical annotation column used to label cell types (string).",
            "k":          "Number of nearest neighbours to consider (integer, default 5)."
        }
    },
    "ripleyL": {
        "description": "Ripley's L curve testing spatial clustering or dispersion of a cell type.",
        "tab_id": "Ripley's L",
        "fields": {
            "annotation":  "Categorical annotation column (string).",
            "cell_type":   "Specific category value to test (string).",
            "max_distance":"Maximum distance for the L-curve (number).",
            "n_simulations":"Number of Monte-Carlo envelope simulations (integer, default 99)."
        }
    }
}


PLOT_JSON_SYSTEM_PROMPT = f"""
## Generating Plot Specifications

When a user asks you to generate, create, or show a plot, you MUST emit a fenced code block with the language tag `plot_json` containing a valid JSON object.

The JSON must follow this structure:
{{
  "plot_type": "<one of the keys below>",
  "parameters": {{
    "<field>": <value>,
    ...
  }},
  "title": "<short human-readable title for this plot>"
}}

Available plot types and their required/optional fields:

{json.dumps(PLOT_SCHEMA, indent=2)}

Rules:
- Only include fields that are relevant; omit optional fields if the user did not mention them.
- Field values must match the dataset (use exact feature names and annotation column names from the app context).
- Place the ```plot_json ... ``` block at the END of your response, after any explanation.
- If the user's request is ambiguous, ask one clarifying question instead of guessing.
- If no plot is being requested, do NOT emit a plot_json block.

Example:
User: "Show me a boxplot of CD3 expression by cell type"
Response:
Sure! Here is a boxplot comparing CD3 expression across the cell types in your dataset.

```plot_json
{{
  "plot_type": "boxplot",
  "parameters": {{
    "feature": "CD3",
    "annotation": "cell_type",
    "show_points": true,
    "log_scale": false
  }},
  "title": "CD3 expression by cell type"
}}
```
"""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def safe_read_html(filename: str, fallback: str = "") -> str:
    try:
        return read_html_file(filename)
    except Exception:
        print(f"[WARNING] Could not load '{filename}' — skipping.")
        return fallback

header_html = safe_read_html("header.html")
footer_html = safe_read_html("footer.html")


def safe_load_data(path: str):
    if not os.path.exists(path):
        print(f"[WARNING] Data file '{path}' not found. "
              f"Set the DATA_PATH environment variable to point to your file. "
              f"You can still load data at runtime via the Data Input tab.")
        return None
    try:
        data = load_data(path)
        print(f"[INFO] Preloaded dataset from '{path}'.")
        return data
    except Exception as e:
        print(f"[ERROR] Failed to load data from '{path}': {e}")
        return None

preloaded_data = safe_load_data(DATA_PATH)


def get_app_context(input, shared) -> str:
    context_parts = []

    try:
        active_tab = input.active_tab()
        context_parts.append(f"User is currently on the '{active_tab}' tab.")
    except Exception:
        pass

    try:
        adata = shared["adata_main"].get()
        if adata is not None:
            context_parts.append(f"Loaded dataset: {adata.shape[0]} cells x {adata.shape[1]} features.")

            if len(adata.obs.columns) > 0:
                context_parts.append(f"Annotation columns: {', '.join(adata.obs.columns.tolist())}.")
                for col in adata.obs.columns:
                    try:
                        if adata.obs[col].dtype.name in ['category', 'object']:
                            val_counts = adata.obs[col].value_counts()
                            top = val_counts.head(10)
                            summary = ", ".join([f"{k} (n={v})" for k, v in top.items()])
                            context_parts.append(
                                f"  - '{col}' has {adata.obs[col].nunique()} unique values: "
                                f"{summary}{'...' if len(val_counts) > 10 else ''}."
                            )
                        else:
                            context_parts.append(
                                f"  - '{col}' is numeric: min={adata.obs[col].min():.3f}, "
                                f"max={adata.obs[col].max():.3f}, mean={adata.obs[col].mean():.3f}."
                            )
                    except Exception:
                        pass

            context_parts.append(f"All features ({adata.shape[1]} total): {', '.join(adata.var_names.tolist())}.")

            if hasattr(adata, 'layers') and len(adata.layers) > 0:
                context_parts.append(f"Available layers: {', '.join(adata.layers.keys())}.")

            if hasattr(adata, 'obsm') and len(adata.obsm) > 0:
                context_parts.append(f"Available embeddings (obsm): {', '.join(adata.obsm.keys())}.")

            if hasattr(adata, 'uns') and len(adata.uns) > 0:
                context_parts.append(f"Unstructured metadata keys: {', '.join(adata.uns.keys())}.")

            spatial_keys = (
                [k for k in adata.obsm.keys() if 'spatial' in k.lower()]
                if hasattr(adata, 'obsm') else []
            )
            if spatial_keys:
                context_parts.append(f"Spatial coordinate keys: {', '.join(spatial_keys)}.")

    except Exception as e:
        context_parts.append(f"Error reading dataset context: {str(e)}")

    return " ".join(context_parts) if context_parts else "No dataset currently loaded."


# ---------------------------------------------------------------------------
# Parse a plot_json block out of an LLM response
# Returns (clean_text, plot_spec_dict | None)
# ---------------------------------------------------------------------------

PLOT_JSON_RE = re.compile(r"```plot_json\s*([\s\S]*?)```", re.IGNORECASE)

def extract_plot_json(raw_text: str):
    """
    Scan `raw_text` for a fenced ```plot_json ... ``` block.
    Returns (text_without_block, parsed_dict) or (original_text, None).
    """
    match = PLOT_JSON_RE.search(raw_text)
    if not match:
        return raw_text, None

    json_str = match.group(1).strip()
    clean_text = PLOT_JSON_RE.sub("", raw_text).strip()

    try:
        spec = json.loads(json_str)
        # Basic validation
        if "plot_type" not in spec or spec["plot_type"] not in PLOT_SCHEMA:
            return clean_text, None
        if "parameters" not in spec or not isinstance(spec["parameters"], dict):
            return clean_text, None
        return clean_text, spec
    except json.JSONDecodeError as e:
        print(f"[WARNING] Could not parse plot_json block: {e}")
        return clean_text, None


# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

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
        id="main_tabs"
    ),

    # ── Floating chat button ──────────────────────────────────────────────
    ui.input_action_button(
        "my_fixed_btn",
        "💬",
        class_="fixed-button btn btn-primary",
        onclick="toggleChatPanel()"
    ),

    # ── Chat panel ───────────────────────────────────────────────────────
    ui.div(
        # Header row
        ui.div(
            ui.h3("Chat with SPAC!", style="margin-top: 0;"),
            ui.tags.button(
                "×",
                onclick="document.getElementById('chat_panel').style.display='none'",
                class_="close-chat-btn",
                type="button"
            ),
            style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px;"
        ),

        # Message history
        ui.div(
            ui.output_ui("chat_display"),
            id="chat_messages",
            style=(
                "max-height: 300px; overflow-y: auto; margin-bottom: 15px; "
                "padding: 10px; background: #f8f9fa; border-radius: 5px;"
            )
        ),

        # ── Plot-JSON action bar (hidden until a spec is available) ───────
        # apply_plot_params — navigates to the tab & pre-fills inputs
        ui.div(
            ui.output_ui("plot_action_bar"),
            id="plot_action_bar_wrapper"
        ),

        # Text input + submit
        ui.input_text_area(
            "user_input",
            "Type your message here:",
            placeholder="Enter text...",
            rows=3,
            width="100%"
        ),
        ui.input_action_button(
            "submit_input", "Submit", class_="btn-primary",
            style="width: 100%; margin-top: 10px;"
        ),

        id="chat_panel",
        class_="chat-panel",
        style="display: none;"
    ),

    # ── JS helpers ────────────────────────────────────────────────────────
    ui.tags.script("""
        function toggleChatPanel() {
            var panel = document.getElementById('chat_panel');
            panel.style.display = (panel.style.display === 'none' || panel.style.display === '') ? 'block' : 'none';
        }
    """),

    ui.tags.script("""
        // Handle programmatic button clicks sent from the server via
        // session.send_custom_message("click_button", {"id": "..."}).
        // A short delay lets Shiny flush updated input values first.
        Shiny.addCustomMessageHandler("click_button", function(msg) {
            setTimeout(function() {
                var btn = document.getElementById(msg.id);
                if (btn) btn.click();
            }, 150);
        });
    """),
    ui.tags.script("""
        Shiny.addCustomMessageHandler("switch_tab", function(msg) {
            setTimeout(function() {
                // Find the tab with matching text and click it
                var tabs = document.querySelectorAll('[role="tab"]');
                tabs.forEach(function(tab) {
                    if (tab.innerText.trim() === msg.tab) tab.click();
                });
            }, 300);
        });
    """),

    ui.tags.script("""
        // Auto-scroll chat messages
        const observer = new MutationObserver(() => {
            const chatMessages = document.getElementById('chat_messages');
            if (chatMessages) chatMessages.scrollTop = chatMessages.scrollHeight;
        });
        setTimeout(() => {
            const chatMessages = document.getElementById('chat_messages');
            if (chatMessages) observer.observe(chatMessages, { childList: true, subtree: true });
        }, 1000);
    """),

    ui.tags.script("""
        // Track the active tab and push value to Shiny
        function getActiveTab() {
            const activeTab = document.querySelector('[role="tab"][aria-selected="true"]');
            if (activeTab) Shiny.setInputValue('active_tab', activeTab.innerText.trim(), {priority: 'event'});
        }
        const tabObserver = new MutationObserver(function(mutations) {
            mutations.forEach(function(mutation) {
                if (mutation.attributeName === 'aria-selected') getActiveTab();
            });
        });
        setTimeout(function() {
            const tabs = document.querySelectorAll('[role="tab"]');
            tabs.forEach(function(tab) { tabObserver.observe(tab, { attributes: true }); });
            getActiveTab();
        }, 1000);
    """),

    # ── Styles ────────────────────────────────────────────────────────────
    ui.tags.style("""
        .fixed-button {
            position: fixed; bottom: 20px; right: 15px; z-index: 1000;
            border-radius: 50%; width: 60px; height: 60px;
            min-width: 60px; min-height: 60px; max-width: 60px; max-height: 60px;
            font-size: 24px; text-align: center; line-height: 60px; padding: 0;
            background: linear-gradient(45deg, #17a2b8, #20c997);
            transition: all 0.3s ease; border: none; outline: none; box-sizing: border-box;
        }
        .fixed-button:hover { transform: scale(1.2); }
        .fixed-button:focus { outline: none !important; box-shadow: none !important; }

        .chat-panel {
            position: fixed; bottom: 90px; right: 15px;
            width: 350px; max-width: 90vw;
            background: white; border-radius: 10px; padding: 20px;
            z-index: 999; box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }
        .close-chat-btn {
            background: none; border: none; font-size: 24px;
            cursor: pointer; padding: 0; width: 30px; height: 30px; color: #666;
        }
        .close-chat-btn:hover { color: #000; transform: scale(1.2); }

        /* Plot action bar */
        .plot-action-bar {
            display: flex; gap: 8px; margin-bottom: 12px;
            padding: 10px; background: #e8f5e9;
            border: 1px solid #a5d6a7; border-radius: 6px;
            flex-wrap: wrap; align-items: center;
        }
        .plot-action-bar .plot-title {
            flex: 1; font-size: 12px; font-weight: 600;
            color: #2e7d32; min-width: 100%;
            margin-bottom: 4px;
        }
        .plot-action-bar .btn { font-size: 12px; padding: 4px 10px; }
    """),

    ui.HTML(footer_html)
)


# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

def server(input, output, session):

    # ── Reactive state ────────────────────────────────────────────────────

    chat_history = reactive.Value([
        {
            "role": "system",
            "content": """You are SPAC (Spatial Proteomics Analysis Companion), a scientific assistant embedded in an interactive analysis tool for spatial omics data.

            ## Your Role
            Help researchers understand, navigate, and interpret their spatial transcriptomics or proteomics data. You are knowledgeable, concise, and always ground your answers in the data or analysis currently available in the app.
            Simplify concepts to make them easily understable for all users including those who have little to no knowledge about the subject matter.

            ## The App's Capabilities
            The app has the following analysis tabs the user can interact with:
            - **Data Input**: Load AnnData (.h5ad) or pickle files containing spatial omics datasets.
            - **Annotations**: View and explore cell-level metadata (e.g., cell type, tissue region, sample ID).
            - **Features**: Explore molecular features (e.g., protein or gene expression levels per cell).
            - **Boxplot**: Compare feature distributions across annotation groups.
            - **Feature vs. Annotation**: Visualise how a continuous feature varies across categorical annotations.
            - **Annotation vs. Annotation**: Explore relationships between two categorical annotation columns.
            - **Spatial**: View cells plotted in their physical tissue coordinates, colored by annotation or feature.
            - **UMAP**: Explore dimensionality-reduced embeddings to identify clusters and structure.
            - **Scatterplot**: Plot any two features against each other to identify correlations.
            - **Nearest Neighbor**: Analyze cellular neighborhoods — which cell types tend to be spatially adjacent.
            - **Ripley's L**: A spatial statistics tool to test whether cell types are clustered, dispersed, or randomly distributed in tissue.

            ## Data Format
            The app works with AnnData objects:
            - `adata.X` — expression/intensity matrix (cells × features)
            - `adata.obs` — per-cell metadata
            - `adata.var` — per-feature metadata
            - `adata.obsm` — multi-dimensional embeddings
            - `adata.layers` — alternative expression matrices
            - `adata.uns` — unstructured metadata

            ## How to Respond
            - Explain plots in plain scientific language.
            - Guide users to the correct tab and inputs.
            - Answer general bioinformatics / spatial biology questions clearly.
            - Keep responses concise unless the user asks for detail.
            - Do not make up results or data values.
            - Use the current app state provided in each message for context-aware answers.

            ## Tone
            Professional but approachable. You are a knowledgeable lab colleague, not a textbook.
            """
                       + PLOT_JSON_SYSTEM_PROMPT
        }
    ])

    # Stores the most recently parsed plot spec dict (or None)
    current_plot_spec = reactive.Value(None)

    # ── Chat display ──────────────────────────────────────────────────────

    @output(suspend_when_hidden=False)
    @render.ui
    def chat_display():
        history = chat_history.get()
        messages = []
        for msg in history[1:]:
            if msg["role"] == "user":
                messages.append(
                    ui.div(
                        ui.strong("You: "),
                        msg["content"],
                        style="margin:8px 0;padding:10px;background:#e3f2fd;border-radius:8px;border-left:3px solid #2196F3;"
                    )
                )
            elif msg["role"] == "assistant":
                messages.append(
                    ui.div(
                        ui.strong("SPAC: "),
                        msg["content"],
                        style="margin:8px 0;padding:10px;background:#f5f5f5;border-radius:8px;border-left:3px solid #4CAF50;"
                    )
                )
        return ui.div(*messages) if messages else ui.p("No messages yet. Start chatting!", style="color:#999;")

    # ── Plot action bar (shown only when a spec is available) ─────────────

    @output
    @render.ui
    def plot_action_bar():
        spec = current_plot_spec.get()
        if spec is None:
            return ui.div()  # render nothing

        plot_type  = spec.get("plot_type", "")
        title      = spec.get("title", plot_type)
        tab_label  = PLOT_SCHEMA.get(plot_type, {}).get("tab_id", "")

        return ui.div(
            ui.span(f"📊 Plot ready: {title}", class_="plot-title"),
            # Apply button – navigates to the target tab and fills inputs
            ui.input_action_button(
                "apply_plot_params",
                f"▶ Apply to '{tab_label}' tab",
                class_="btn btn-sm btn-success"
            ),
            # Clear button
            ui.input_action_button(
                "clear_plot_spec",
                "✕ Clear",
                class_="btn btn-sm btn-outline-secondary"
            ),
            class_="plot-action-bar"
        )

    # ── Apply plot params → navigate to tab + pre-fill inputs ─────────────

    @reactive.effect
    @reactive.event(input.apply_plot_params)
    def _apply_plot_params():
        spec = current_plot_spec.get()
        if spec is None:
            return

        plot_type  = spec["plot_type"]
        params     = spec.get("parameters", {})
        tab_label  = PLOT_SCHEMA.get(plot_type, {}).get("tab_id")

        # Navigate to the correct tab
        if tab_label:
            session.send_custom_message("switch_tab", {"tab": tab_label})

        # ------------------------------------------------------------------
        # Pre-fill inputs for each plot type.
        # The input IDs below must match the actual Shiny input IDs used in
        # your server/* modules.  Adjust them to your codebase as needed.
        # ------------------------------------------------------------------
        try:
            if plot_type == "boxplot":
                adata = shared["adata_main"].get()

                # Rebuild choices from live data before setting selected
                if adata is not None:
                    anno_choices = ["No Annotation"] + list(adata.obs.columns)
                    layer_choices = ["Original"] + list(adata.layers.keys())
                    feature_choices = list(adata.var_names)

                    anno_val = params.get("annotation", "No Annotation")
                    ui.update_select("bp_anno",
                                     choices=anno_choices,
                                     selected=anno_val
                                     )

                    features = params.get("features", [])
                    if isinstance(features, str):
                        features = [features]
                    ui.update_selectize("bp_features",
                                        choices=feature_choices,
                                        selected=features
                                        )

                    layer_val = params.get("layer", "Original")
                    ui.update_select("bp_layer",
                                     choices=layer_choices,
                                     selected=layer_val
                                     )
            elif plot_type == "feat_vs_anno":
                if "feature"    in params: ui.update_select("fva_feature",    selected=params["feature"])
                if "annotation" in params: ui.update_select("fva_annotation", selected=params["annotation"])

            elif plot_type == "anno_vs_anno":
                subtype = params.get("plot_subtype", "sankey")
                src = params.get("source_annotation")
                tgt = params.get("target_annotation")

                if subtype == "relational_heatmap":
                    if src: ui.update_select("rhm_anno1", selected=src)
                    if tgt: ui.update_select("rhm_anno2", selected=tgt)
                    session.send_custom_message("click_button", {"id": "go_rhm1"})
                else:  # default: sankey
                    if src: ui.update_select("sk1_anno1", selected=src)
                    if tgt: ui.update_select("sk1_anno2", selected=tgt)
                    session.send_custom_message("click_button", {"id": "go_sk1"})

            elif plot_type == "spatial":
                if "color_by"    in params: ui.update_select("spatial_color_by",    selected=params["color_by"])
                if "point_size"  in params: ui.update_slider("spatial_point_size",  value=params["point_size"])
                if "alpha"       in params: ui.update_slider("spatial_alpha",        value=params["alpha"])
                if "spatial_key" in params: ui.update_select("spatial_key",          selected=params["spatial_key"])

            elif plot_type == "umap":
                if "color_by"   in params: ui.update_select("umap_color_by",   selected=params["color_by"])
                if "umap_key"   in params: ui.update_select("umap_key",         selected=params["umap_key"])
                if "point_size" in params: ui.update_slider("umap_point_size", value=params["point_size"])

            elif plot_type == "scatterplot":
                if "feature_x" in params: ui.update_select("scatter_feature_x", selected=params["feature_x"])
                if "feature_y" in params: ui.update_select("scatter_feature_y", selected=params["feature_y"])
                if "color_by"  in params: ui.update_select("scatter_color_by",  selected=params["color_by"])
                if "log_x"     in params: ui.update_checkbox("scatter_log_x",   value=params["log_x"])
                if "log_y"     in params: ui.update_checkbox("scatter_log_y",   value=params["log_y"])

            elif plot_type == "nearest_neighbor":
                if "annotation" in params: ui.update_select("nn_annotation", selected=params["annotation"])
                if "k"          in params: ui.update_slider("nn_k",           value=params["k"])

            elif plot_type == "ripleyL":
                if "annotation"    in params: ui.update_select("ripley_annotation",    selected=params["annotation"])
                if "cell_type"     in params: ui.update_select("ripley_cell_type",     selected=params["cell_type"])
                if "max_distance"  in params: ui.update_slider("ripley_max_distance",  value=params["max_distance"])
                if "n_simulations" in params: ui.update_slider("ripley_n_simulations", value=params["n_simulations"])

        except Exception as e:
            ui.notification_show(
                f"Could not pre-fill all inputs: {e}. Please set any remaining fields manually.",
                type="warning",
                duration=6
            )

    # ── Clear plot spec ───────────────────────────────────────────────────

    @reactive.effect
    @reactive.event(input.clear_plot_spec)
    def _clear_plot_spec():
        current_plot_spec.set(None)

    # ── Main chat submission handler ──────────────────────────────────────

    @reactive.effect
    @reactive.event(input.submit_input)
    async def handle_submission():
        user_text = input.user_input()
        if not user_text:
            return

        ui.update_text_area("user_input", value="")
        history = list(chat_history.get())
        history.append({"role": "user", "content": user_text})
        chat_history.set(history)

        relevant_context = retrieve_context(user_text, n_results=RAG_N_RESULTS)
        app_context      = get_app_context(input, shared)

        rag_message = {
            "role": "system",
            "content": (
                f"Current app state: {app_context}\n\n"
                f"Relevant sections from the SPAC paper for this query:\n\n{relevant_context}\n\n"
                "Use the above as your primary reference when answering."
                "IMPORTANT: Only emit a plot_json block if the user is explicitly "
                "asking for a NEW plot. Do not repeat or reference previous plot specs."
            )
        }
        messages_to_send = [history[0], rag_message] + history[1:]

        try:
            async with httpx.AsyncClient(timeout=OLLAMA_TIMEOUT) as client:
                response = await client.post(
                    f"{OLLAMA_URL}/api/chat",
                    json={"model": OLLAMA_MODEL, "messages": messages_to_send, "stream": False}
                )
                response.raise_for_status()
                data = response.json()
                raw_reply = data["message"]["content"]

            # ── Parse out any plot_json block ──────────────────────────
            clean_reply, plot_spec = extract_plot_json(raw_reply)

            if plot_spec is not None:
                current_plot_spec.set(plot_spec)
                # Append a note so the user knows a spec was generated
                plot_note = (
                    f"\n\n*(Plot spec generated for **{plot_spec.get('plot_type')}**. "
                    "Use the action bar below to apply it directly to the app.)*"
                )
                clean_reply = clean_reply + plot_note

            history = list(chat_history.get())
            history.append({"role": "assistant", "content": clean_reply})
            chat_history.set(history)

        except httpx.ConnectError:
            ui.notification_show(
                f"Could not connect to Ollama at {OLLAMA_URL}. "
                "Is Ollama running? If running outside Docker, set OLLAMA_URL=http://localhost:11434.",
                type="error", duration=8
            )
        except httpx.TimeoutException:
            ui.notification_show(
                f"Request to Ollama timed out after {OLLAMA_TIMEOUT}s. "
                "Try increasing OLLAMA_TIMEOUT or using a smaller model.",
                type="error", duration=8
            )
        except Exception as e:
            ui.notification_show(f"Unexpected error: {str(e)}", type="error", duration=5)

    # ── Shared data state (unchanged from original) ───────────────────────

    data_loaded = reactive.Value(False)
    adata_main  = reactive.Value(preloaded_data)

    data_keys = [
        "X_data", "obs_data", "obsm_data", "layers_data", "var_data", "uns_data",
        "shape_data", "obs_names", "obsm_names", "layers_names", "var_names",
        "uns_names", "spatial_distance_columns", "df_heatmap", "df_relational",
        "df_boxplot", "df_histogram2", "df_histogram1", "df_nn", "df_ripley"
    ]

    shared = {
        "preloaded_data": preloaded_data,
        "data_loaded":    data_loaded,
        "adata_main":     adata_main,
    }
    for key in data_keys:
        shared[key] = reactive.Value(None)

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


static_path = os.path.join(os.path.dirname(__file__), "www")
app = App(app_ui, server, static_assets=static_path)