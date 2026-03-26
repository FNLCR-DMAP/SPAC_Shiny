from shiny import App, ui, reactive, render
import os
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


header_html = read_html_file("header.html")
footer_html = read_html_file("footer.html")

file_path = "dev_example.pickle"
preloaded_data = load_data(file_path)


def get_app_context(input, shared) -> str:
    context_parts = []

    try:
        active_tab = input.active_tab()
        context_parts.append(f"User is currently on the '{active_tab}' tab.")
    except:
        pass

    try:
        adata = shared["adata_main"].get()
        if adata is not None:
            # Shape
            context_parts.append(f"Loaded dataset: {adata.shape[0]} cells x {adata.shape[1]} features.")

            # All annotation columns + their unique values and counts
            if len(adata.obs.columns) > 0:
                context_parts.append(f"Annotation columns: {', '.join(adata.obs.columns.tolist())}.")
                for col in adata.obs.columns:
                    try:
                        if adata.obs[col].dtype.name in ['category', 'object']:
                            val_counts = adata.obs[col].value_counts()
                            top = val_counts.head(10)
                            summary = ", ".join([f"{k} (n={v})" for k, v in top.items()])
                            context_parts.append(f"  - '{col}' has {adata.obs[col].nunique()} unique values: {summary}{'...' if len(val_counts) > 10 else ''}.")
                        else:
                            # Numeric annotation
                            context_parts.append(
                                f"  - '{col}' is numeric: min={adata.obs[col].min():.3f}, "
                                f"max={adata.obs[col].max():.3f}, mean={adata.obs[col].mean():.3f}."
                            )
                    except:
                        pass

            # All features
            context_parts.append(f"All features ({adata.shape[1]} total): {', '.join(adata.var_names.tolist())}.")

            # Layers
            if hasattr(adata, 'layers') and len(adata.layers) > 0:
                context_parts.append(f"Available layers: {', '.join(adata.layers.keys())}.")

            # Embeddings (obsm)
            if hasattr(adata, 'obsm') and len(adata.obsm) > 0:
                context_parts.append(f"Available embeddings (obsm): {', '.join(adata.obsm.keys())}.")

            # Unstructured metadata keys
            if hasattr(adata, 'uns') and len(adata.uns) > 0:
                context_parts.append(f"Unstructured metadata keys: {', '.join(adata.uns.keys())}.")

            # Spatial coordinates availability
            spatial_keys = [k for k in adata.obsm.keys() if 'spatial' in k.lower()] if hasattr(adata, 'obsm') else []
            if spatial_keys:
                context_parts.append(f"Spatial coordinate keys: {', '.join(spatial_keys)}.")

    except Exception as e:
        context_parts.append(f"Error reading dataset context: {str(e)}")

    return " ".join(context_parts) if context_parts else "No dataset currently loaded."


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

    ui.input_action_button(
        "my_fixed_btn",
        "💬",
        class_="fixed-button btn btn-primary",
        onclick="toggleChatPanel()"
    ),

    ui.div(
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
        ui.div(
            ui.output_ui("chat_display"),
            id="chat_messages",
            style="max-height: 300px; overflow-y: auto; margin-bottom: 15px; padding: 10px; background: #f8f9fa; border-radius: 5px;"
        ),
        ui.input_text_area(
            "user_input",
            "Type your message here:",
            placeholder="Enter text...",
            rows=3,
            width="100%"
        ),
        ui.input_action_button("submit_input", "Submit", class_="btn-primary", style="width: 100%; margin-top: 10px;"),
        id="chat_panel",
        class_="chat-panel",
        style="display: none;"
    ),

    ui.tags.script("""
        function toggleChatPanel() {
            var panel = document.getElementById('chat_panel');
            if (panel.style.display === 'none' || panel.style.display === '') {
                panel.style.display = 'block';
            } else {
                panel.style.display = 'none';
            }
        }
    """),

    ui.tags.script("""
        const observer = new MutationObserver(() => {
            const chatMessages = document.getElementById('chat_messages');
            if (chatMessages) {
                chatMessages.scrollTop = chatMessages.scrollHeight;
            }
        });
        
        setTimeout(() => {
            const chatMessages = document.getElementById('chat_messages');
            if (chatMessages) {
                observer.observe(chatMessages, { childList: true, subtree: true });
            }
        }, 1000);
    """),

    ui.tags.style("""
        .fixed-button {
            position: fixed;
            bottom: 20px;
            right: 15px;
            z-index: 1000;
            border-radius: 50%;
            width: 60px;
            height: 60px;
            min-width: 60px;
            min-height: 60px;
            max-width: 60px;
            max-height: 60px;
            font-size: 24px;
            text-align: center;
            line-height: 60px;
            padding: 0;
            background: linear-gradient(45deg, #17a2b8, #20c997);
            transition: all 0.3s ease;
            border: none;
            outline: none;
            box-sizing: border-box;
        }
        .fixed-button:hover {
            transform: scale(1.2);
            outline: none !important;
            box-shadow: none !important;
        }
        .fixed-button:focus {
            outline: none !important;
            box-shadow: none !important;
        }
        .chat-panel {
            position: fixed;
            bottom: 90px;
            right: 15px;
            width: 350px;
            max-width: 90vw;
            background: white;
            border-radius: 10px;
            padding: 20px;
            z-index: 999;
            box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        }
        .close-chat-btn {
            background: none;
            border: none;
            font-size: 24px;
            cursor: pointer;
            padding: 0;
            width: 30px;
            height: 30px;
            color: #666666;
        }
        .close-chat-btn:hover {
            color: #000000;
            transform: scale(1.2);
        }
    """),
    ui.tags.script("""
    function getActiveTab() {
        const activeTab = document.querySelector('[role="tab"][aria-selected="true"]');
        if (activeTab) {
            Shiny.setInputValue('active_tab', activeTab.innerText.trim(), {priority: 'event'});
        }
    }

    const tabObserver = new MutationObserver(function(mutations) {
        mutations.forEach(function(mutation) {
            if (mutation.attributeName === 'aria-selected') {
                getActiveTab();
            }
        });
    });

    setTimeout(function() {
        const tabs = document.querySelectorAll('[role="tab"]');
        tabs.forEach(function(tab) {
            tabObserver.observe(tab, { attributes: true });
        });
        getActiveTab();
    }, 1000);
"""),
    ui.HTML(footer_html)
)

def server(input, output, session):
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
            - **Feature vs. Annotation**: Visualize how a continuous feature varies across categorical annotations.
            - **Annotation vs. Annotation**: Explore relationships between two categorical annotation columns.
            - **Spatial**: View cells plotted in their physical tissue coordinates, colored by annotation or feature.
            - **UMAP**: Explore dimensionality-reduced embeddings to identify clusters and structure.
            - **Scatterplot**: Plot any two features against each other to identify correlations.
            - **Nearest Neighbor**: Analyze cellular neighborhoods — which cell types tend to be spatially adjacent.
            - **Ripley's L**: A spatial statistics tool to test whether cell types are clustered, dispersed, or randomly distributed in tissue.
            
            ## Data Format
            The app works with AnnData objects, a standard format in single-cell and spatial omics analysis:
            - `adata.X` — the main expression/intensity matrix (cells × features)
            - `adata.obs` — per-cell metadata (annotations like cell type, region, sample)
            - `adata.var` — per-feature metadata (e.g., gene/protein names)
            - `adata.obsm` — multi-dimensional embeddings (e.g., spatial coordinates, UMAP)
            - `adata.layers` — alternative expression matrices (e.g., normalized, raw counts)
            - `adata.uns` — unstructured metadata (e.g., color palettes, analysis results)
            
            ## How to Respond
            - If a user asks what a plot or analysis means, explain it in plain scientific language.
            - If a user asks how to do something in the app, guide them to the correct tab and inputs.
            - If a user asks a general bioinformatics or spatial biology question, answer clearly and accurately.
            - If a user asks something outside your scope (e.g., unrelated coding questions), politely redirect them.
            - Keep responses concise unless the user asks for detail. Avoid unnecessary jargon.
            - Do not make up results or data values — if you don't know the current state of the data, say so and suggest they check the relevant tab.
            - Use the current app state provided in each message to give context-aware answers.
            
            ## Tone
            Professional but approachable. You are a knowledgeable lab colleague, not a textbook.
            """
        }
    ])

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
                        style="margin: 8px 0; padding: 10px; background: #e3f2fd; border-radius: 8px; border-left: 3px solid #2196F3;"
                    )
                )
            elif msg["role"] == "assistant":
                messages.append(
                    ui.div(
                        ui.strong("SPAC: "),
                        msg["content"],
                        style="margin: 8px 0; padding: 10px; background: #f5f5f5; border-radius: 8px; border-left: 3px solid #4CAF50;"
                    )
                )

        return ui.div(*messages) if messages else ui.p(
            "No messages yet. Start chatting!", style="color: #999;"
        )

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

        relevant_context = retrieve_context(user_text, n_results=3)
        app_context = get_app_context(input, shared)

        rag_message = {
            "role": "system",
            "content": f"""Current app state: {app_context}

Relevant sections from the SPAC paper for this query:

{relevant_context}

Use the above as your primary reference when answering."""
        }
        messages_to_send = [history[0], rag_message] + history[1:]

        try:
            async with httpx.AsyncClient(timeout=120) as client:
                response = await client.post(
                    "http://host.docker.internal:11434/api/chat",
                    json={
                        "model": "gemma3:latest",
                        "messages": messages_to_send,
                        "stream": False
                    }
                )
                response.raise_for_status()
                data = response.json()

                bot_reply = data["message"]["content"]

                history = list(chat_history.get())
                history.append({"role": "assistant", "content": bot_reply})
                chat_history.set(history)

        except Exception as e:
            ui.notification_show(f"Error: {str(e)}", type="error", duration=5)

    data_loaded = reactive.Value(False)
    adata_main = reactive.Value(preloaded_data)

    data_keys = [
        "X_data", "obs_data", "obsm_data", "layers_data", "var_data", "uns_data",
        "shape_data", "obs_names", "obsm_names", "layers_names", "var_names",
        "uns_names", "spatial_distance_columns", "df_heatmap", "df_relational",
        "df_boxplot", "df_histogram2", "df_histogram1", "df_nn", "df_ripley"
    ]

    shared = {
        "preloaded_data": preloaded_data,
        "data_loaded": data_loaded,
        "adata_main": adata_main,
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