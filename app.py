from shiny import App, ui, reactive, render
import os
import requests


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


# Read header and footer HTML content
header_html = read_html_file("header.html")
footer_html = read_html_file("footer.html")

file_path = "dev_example.pickle"
preloaded_data = load_data(file_path)


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

    ui.HTML(footer_html)
)

def server(input, output, session):
    chat_history = reactive.Value([
        {
            "role": "system",
            "content": "You are SPAC, a helpful scientific assistant."
        }
    ])

    @output
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

        return ui.div(*messages) if messages else ui.p("No messages yet. Start chatting!", style="color: #999;")

    @reactive.effect
    @reactive.event(input.submit_input)
    def handle_submission():
        user_text = input.user_input()

        if not user_text:
            return

        # Add user message
        history = chat_history.get().copy()
        history.append({"role": "user", "content": user_text})
        chat_history.set(history)

        try:
            response = requests.post(
                "http://host.docker.internal:11434/api/chat",
                json={
                    "model": "gemma3:latest",
                    "messages": history,
                    "stream": False
                },
                timeout=120
            )

            print("STATUS:", response.status_code)
            print("RAW:", response.text[:300])

            response.raise_for_status()

            data = response.json()
            bot_reply = data.get("message", {}).get("content", "No response")

            history = chat_history.get().copy()
            history.append({"role": "assistant", "content": bot_reply})
            chat_history.set(history)

        except Exception as e:
            print("ERROR:", e)
            ui.notification_show(f"Error: {str(e)}", type="error", duration=5)

        ui.update_text_area("user_input", value="")
    # Data initialization - THIS WAS IN THE WRONG PLACE
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

    # Individual server components
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