# SPAC Shiny - Reactive Patterns Quick Reference

A condensed guide for quick lookup while developing.

## Table of Contents
- [Reactive Decorators Cheat Sheet](#reactive-decorators-cheat-sheet)
- [Common Patterns](#common-patterns)
- [Shared State API](#shared-state-api)
- [Quick Troubleshooting](#quick-troubleshooting)

---

## Reactive Decorators Cheat Sheet

| Decorator | Purpose | When to Use | Example |
|-----------|---------|-------------|---------|
| `@reactive.Value()` | Mutable state container | Store app state | `data = reactive.Value(None)` |
| `@reactive.calc` | Cached computation | Expensive operations | Compute derived values |
| `@reactive.effect` | Side effect (lowercase) | Update UI, logging | Synchronize UI elements |
| `@reactive.Effect` | Side effect (uppercase) | Same as effect | Synchronize UI elements |
| `@reactive.event()` | Event-triggered reaction | Button clicks only | Wait for user action |
| `@render.plot` | Render plot output | Display matplotlib/plotly | Create visualizations |
| `@render.text` | Render text output | Display strings | Show computed text |
| `@render.ui` | Render dynamic UI | Conditional UI elements | Create dynamic inputs |
| `@render_widget` | Render plotly widget | Interactive plotly plots | Display interactive viz |
| `@render.download()` | File download | CSV/data export | Provide data downloads |

---

## Common Patterns

### Pattern 1: Read/Write Shared State

```python
# ✅ Reading
value = shared['adata_main'].get()

# ✅ Writing
shared['adata_main'].set(new_value)

# ❌ Don't do this
value = shared['adata_main']  # Wrong!
shared['adata_main'].value = x  # Wrong!
```

### Pattern 2: Update UI Input Choices

```python
@reactive.Effect
def update_dropdown():
    choices = shared['obs_names'].get()
    ui.update_select("my_select", choices=choices)
    
    # Set default selection
    if choices and len(choices) > 0:
        ui.update_select("my_select", selected=choices[0])
```

### Pattern 3: Button-Triggered Plot

```python
@output
@render.plot
@reactive.event(input.go_button, ignore_none=True)
def my_plot():
    """Only runs when button clicked."""
    adata = reconstruct_adata(shared)
    fig = create_plot(adata)
    
    # Store data for download
    shared['df_output'].set(extract_data(fig))
    
    return fig
```

### Pattern 4: Cached Computation

```python
@reactive.calc
def expensive_computation():
    """Result cached until dependencies change."""
    data = shared['X_data'].get()
    return costly_transform(data)

# Reuse across multiple renders
@output
@render.plot
def plot1():
    result = expensive_computation()  # Uses cache
    return plot(result)

@output
@render.plot
def plot2():
    result = expensive_computation()  # Same cache!
    return different_plot(result)
```

### Pattern 5: Conditional UI Insertion

```python
ui_initialized = reactive.Value(False)

@reactive.effect
def dynamic_ui():
    show = input.show_checkbox()
    initialized = ui_initialized.get()
    
    if show and not initialized:
        # Insert UI
        new_element = ui.input_select("dynamic_input", "Label", choices=[])
        ui.insert_ui(
            ui.div({"id": "inserted-element"}, new_element),
            selector="#target-container",
            where="beforeEnd"
        )
        ui_initialized.set(True)
        
    elif not show and initialized:
        # Remove UI
        ui.remove_ui("#inserted-element")
        ui_initialized.set(False)
```

### Pattern 6: Reconstruct AnnData

```python
@reactive.calc
def get_adata():
    """Reconstruct AnnData from shared components."""
    return ad.AnnData(
        X=shared['X_data'].get(),
        obs=pd.DataFrame(shared['obs_data'].get()),
        var=pd.DataFrame(shared['var_data'].get()),
        layers=shared['layers_data'].get(),
        dtype=shared['X_data'].get().dtype
    )
```

### Pattern 7: CSV Download

```python
@render.download(filename="data.csv")
def download_data():
    """Provide CSV download of stored dataframe."""
    df = shared['df_boxplot'].get()
    if df is not None:
        csv_string = df.to_csv(index=False)
        csv_bytes = csv_string.encode("utf-8")
        return csv_bytes, "text/csv"
    return None
```

### Pattern 8: Input Validation

```python
from shiny import req

@output
@render.plot
@reactive.event(input.go_btn, ignore_none=True)
def my_plot():
    # Ensure inputs are available before proceeding
    req(input.feature_select())
    req(shared['adata_main'].get())
    
    # Safe to use inputs now
    return create_plot()
```

---

## Shared State API

### Data Components

| Key | Type | Description |
|-----|------|-------------|
| `adata_main` | AnnData | Main data object |
| `X_data` | np.ndarray | Expression matrix |
| `obs_data` | pd.DataFrame | Cell annotations |
| `var_data` | pd.DataFrame | Gene/feature info |
| `layers_data` | dict | Data layers |
| `uns_data` | dict | Unstructured data |
| `obsm_data` | dict | Multi-dimensional annotations |

### Metadata

| Key | Type | Description |
|-----|------|-------------|
| `obs_names` | list | Annotation column names |
| `var_names` | list | Feature names |
| `layers_names` | list | Layer names |
| `obsm_names` | list | Obsm keys |
| `uns_names` | list | Uns keys |
| `shape_data` | tuple | Data dimensions |
| `spatial_distance_columns` | list | Spatial analysis columns |

### Output DataFrames

| Key | Type | Description |
|-----|------|-------------|
| `df_boxplot` | pd.DataFrame | Boxplot data |
| `df_histogram1` | pd.DataFrame | Histogram data (features) |
| `df_histogram2` | pd.DataFrame | Histogram data (annotations) |
| `df_heatmap` | pd.DataFrame | Heatmap data |
| `df_relational` | pd.DataFrame | Relational plot data |
| `df_nn` | pd.DataFrame | Nearest neighbor data |
| `df_ripley` | pd.DataFrame | Ripley's L data |

### State Flags

| Key | Type | Description |
|-----|------|-------------|
| `data_loaded` | bool | Whether data is loaded |
| `preloaded_data` | AnnData | Initial demo data |

---

## Quick Troubleshooting

### Problem: Dropdown not updating after data load

**Cause**: Effect not watching correct reactive value

**Solution**:
```python
# In effect_update_server.py
@reactive.Effect
def update_my_dropdown():
    choices = shared['obs_names'].get()  # Watch this value
    ui.update_select("my_dropdown", choices=choices)
```

---

### Problem: Plot not refreshing when data changes

**Cause**: Using `@reactive.event` blocks auto-refresh OR reading value outside render function

**Solution**:
```python
# Option 1: Remove @reactive.event to make auto-reactive
@output
@render.plot
def my_plot():
    data = shared['adata_main'].get()  # Read inside
    return plot(data)

# Option 2: Keep event but read fresh data
@output
@render.plot
@reactive.event(input.go_btn, ignore_none=True)
def my_plot():
    data = shared['adata_main'].get()  # Read inside
    return plot(data)
```

---

### Problem: "Infinite reactive loop" or app freezes

**Cause**: Effect reads and writes same reactive value

**Solution**:
```python
# ❌ BAD: Creates loop
@reactive.Effect
def bad():
    val = shared['counter'].get()
    shared['counter'].set(val + 1)  # Triggers itself!

# ✅ GOOD: Use event to break loop
@reactive.effect
@reactive.event(input.increment_btn)
def good():
    val = shared['counter'].get()
    shared['counter'].set(val + 1)  # Only on button
```

---

### Problem: Expensive computation running too often

**Cause**: Not using `@reactive.calc` for caching

**Solution**:
```python
# ❌ BAD: Runs every time
@output
@render.plot
def plot1():
    data = shared['X_data'].get()
    result = expensive_transform(data)  # Recomputes!
    return plot(result)

# ✅ GOOD: Cached
@reactive.calc
def transformed_data():
    data = shared['X_data'].get()
    return expensive_transform(data)  # Cached!

@output
@render.plot
def plot1():
    return plot(transformed_data())  # Uses cache
```

---

### Problem: Error "NoneType has no attribute..."

**Cause**: Reactive value is None when accessed

**Solution**:
```python
# Use req() to ensure value exists
from shiny import req

@output
@render.plot
def my_plot():
    data = shared['adata_main'].get()
    req(data)  # Stops execution if None
    
    # Safe to use data now
    return plot(data)
```

---

### Problem: Download button downloads wrong/stale data

**Cause**: Generating new data instead of using stored data OR data not stored

**Solution**:
```python
# ✅ Store data when plot is generated
@output
@render.plot
@reactive.event(input.go_btn, ignore_none=True)
def my_plot():
    fig, df = create_plot_with_data()
    shared['df_output'].set(df)  # Store it!
    return fig

# ✅ Retrieve stored data for download
@render.download(filename="data.csv")
def download():
    df = shared['df_output'].get()  # Get stored data
    if df is not None:
        return df.to_csv(index=False).encode("utf-8"), "text/csv"
    return None
```

---

### Problem: UI element inserted multiple times

**Cause**: Not tracking insertion state

**Solution**:
```python
# Track whether UI is inserted
ui_inserted = reactive.Value(False)

@reactive.effect
def insert_ui():
    show = input.show_checkbox()
    inserted = ui_inserted.get()
    
    if show and not inserted:
        # Insert only if not already inserted
        ui.insert_ui(...)
        ui_inserted.set(True)
        
    elif not show and inserted:
        # Remove only if inserted
        ui.remove_ui(...)
        ui_inserted.set(False)
```

---

## Quick Module Template

### Minimal New Feature Module

**UI Module** (`ui/my_feature_ui.py`):
```python
from shiny import ui

def my_feature_ui():
    return ui.nav_panel(
        "My Feature",
        ui.layout_sidebar(
            ui.sidebar(
                ui.input_select("mf_anno", "Annotation:", choices=[]),
                ui.input_action_button("go_mf", "Generate", class_="btn-primary")
            ),
            ui.output_plot("mf_plot")
        )
    )
```

**Server Module** (`server/my_feature_server.py`):
```python
from shiny import render, reactive
import anndata as ad
import pandas as pd

def my_feature_server(input, output, session, shared):
    @output
    @render.plot
    @reactive.event(input.go_mf, ignore_none=True)
    def mf_plot():
        adata = ad.AnnData(
            X=shared['X_data'].get(),
            obs=pd.DataFrame(shared['obs_data'].get()),
            var=pd.DataFrame(shared['var_data'].get())
        )
        return create_my_plot(adata)
```

**Effect Update** (`server/effect_update_server.py`):
```python
@reactive.Effect
def update_my_feature():
    choices = shared['obs_names'].get()
    ui.update_select("mf_anno", choices=choices)
```

**Register in app.py**:
```python
# Import
from ui import my_feature_ui
from server import my_feature_server

# Add to UI
ui.navset_card_tab(
    # ... existing
    my_feature_ui(),
)

# Add to server
my_feature_server(input, output, session, shared)
```

---

## Key Takeaways

1. **Always use `.get()` and `.set()`** for reactive values
2. **Use `@reactive.calc`** for expensive computations
3. **Use `@reactive.effect`** for UI updates
4. **Use `@reactive.event`** for button-triggered actions
5. **Store output data** for downloads in shared state
6. **Track UI insertion state** to prevent duplicates
7. **Read reactive values inside** render functions
8. **Validate inputs** with `req()` before using

---

**For full documentation, see**: [REACTIVE_ARCHITECTURE.md](REACTIVE_ARCHITECTURE.md)
