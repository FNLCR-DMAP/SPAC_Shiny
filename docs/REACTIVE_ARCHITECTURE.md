# SPAC Shiny App - Reactive Architecture Guide

## Table of Contents
1. [Overview](#overview)
2. [Architecture Diagram](#architecture-diagram)
3. [Core Reactive Components](#core-reactive-components)
4. [Reactive Patterns](#reactive-patterns)
5. [Data Flow](#data-flow)
6. [Module Communication](#module-communication)
7. [Developer Guide](#developer-guide)
8. [Best Practices](#best-practices)

---

## Overview

The SPAC Shiny application uses **Shiny for Python's reactive programming framework** to create an interactive single-cell analysis platform. The architecture is built on a **shared state pattern** where reactive values propagate changes automatically throughout the application.

### Key Principles

- **Modular Architecture**: UI and server logic are separated into dedicated modules
- **Shared State**: A central `shared` dictionary contains reactive values accessible to all modules
- **Automatic Updates**: Changes to reactive values trigger dependent computations and UI updates
- **Separation of Concerns**: Effect updates are isolated in a dedicated module

---

## Architecture Diagram

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              SPAC Shiny Application                          │
│                                                                               │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │                           app.py (Main)                              │   │
│  │  • Initializes shared reactive state                                 │   │
│  │  • Coordinates all UI and server modules                             │   │
│  │  • Manages app lifecycle                                             │   │
│  └────────────────────────┬────────────────────────────────────────────┘   │
│                           │                                                  │
│         ┌─────────────────┴─────────────────┐                               │
│         │                                     │                               │
│         ▼                                     ▼                               │
│  ┌─────────────┐                      ┌─────────────┐                       │
│  │   UI Layer  │                      │Server Layer │                       │
│  │  (ui/*.py)  │                      │(server/*.py)│                       │
│  └─────────────┘                      └─────────────┘                       │
│         │                                     │                               │
│         │                                     │                               │
└─────────┼─────────────────────────────────────┼───────────────────────────────┘
          │                                     │
          │     ┌───────────────────────────────┴────────────┐
          │     │                                              │
          │     ▼                                              ▼
    ┌─────┴─────────────────┐              ┌──────────────────────────────┐
    │   User Interactions    │              │    Shared Reactive State     │
    │   • File uploads       │──────────────▶│  (Central State Store)       │
    │   • Button clicks      │              │                              │
    │   • Input selections   │              │  ┌────────────────────────┐ │
    │   • Parameter changes  │              │  │ Data Reactive Values:  │ │
    └────────────────────────┘              │  │ • adata_main           │ │
                                            │  │ • X_data               │ │
                                            │  │ • obs_data             │ │
                                            │  │ • var_data             │ │
                                            │  │ • layers_data          │ │
                                            │  │ • uns_data             │ │
           Triggers                         │  │ • obs_names            │ │
              │                             │  │ • var_names            │ │
              ▼                             │  │ • layers_names         │ │
    ┌──────────────────────┐               │  │ • spatial_distance_    │ │
    │  data_input_server   │───────────────▶│ │   columns              │ │
    │  • Loads data files  │               │  │                        │ │
    │  • Updates adata_main│               │  └────────────────────────┘ │
    │  • Extracts parts    │               │                              │
    └──────────────────────┘               │  ┌────────────────────────┐ │
              │                             │  │ Output Data Frames:    │ │
              │                             │  │ • df_boxplot           │ │
              │ Propagates Changes          │  │ • df_histogram1        │ │
              ▼                             │  │ • df_histogram2        │ │
    ┌──────────────────────┐               │  │ • df_heatmap           │ │
    │ effect_update_server │◀──────────────┤  │ • df_relational        │ │
    │  • Watches reactive  │               │  │ • df_nn                │ │
    │    value changes     │               │  │ • df_ripley            │ │
    │  • Updates UI inputs │               │  │                        │ │
    │  • Syncs dropdowns   │               │  └────────────────────────┘ │
    └──────────────────────┘               │                              │
              │                             │  ┌────────────────────────┐ │
              │ Updates available to        │  │ State Flags:           │ │
              ▼                             │  │ • data_loaded          │ │
    ┌──────────────────────┐               │  │ • preloaded_data       │ │
    │ Feature Server       │◀──────────────┤  └────────────────────────┘ │
    │ Modules:             │               │                              │
    │ • annotations_server │               └──────────────────────────────┘
    │ • features_server    │                              │
    │ • boxplot_server     │                              │
    │ • spatial_server     │              Reads/Writes    │
    │ • umap_server        │◀─────────────────────────────┘
    │ • scatterplot_server │
    │ • nearest_neighbor_  │
    │   server             │
    │ • ripleyL_server     │
    │ • anno_vs_anno_      │
    │   server             │
    │ • feat_vs_anno_      │
    │   server             │
    └──────────────────────┘
              │
              │ Generates
              ▼
    ┌──────────────────────┐
    │   Outputs/Renders    │
    │   • Plots            │
    │   • Tables           │
    │   • Downloads        │
    │   • Dynamic UI       │
    └──────────────────────┘
```

---

## Core Reactive Components

### 1. Shared State Dictionary

The `shared` dictionary in `app.py` is the **central nervous system** of the application. It contains all reactive values that modules use to communicate.

```python
# In app.py
shared = {
    # Preloaded data
    "preloaded_data": preloaded_data,
    
    # State tracking
    "data_loaded": reactive.Value(False),
    
    # Main data object
    "adata_main": reactive.Value(adata_main),
    
    # AnnData components (dynamically created)
    "X_data": reactive.Value(None),
    "obs_data": reactive.Value(None),
    "var_data": reactive.Value(None),
    "layers_data": reactive.Value(None),
    "uns_data": reactive.Value(None),
    "obsm_data": reactive.Value(None),
    
    # Metadata
    "shape_data": reactive.Value(None),
    "obs_names": reactive.Value(None),
    "var_names": reactive.Value(None),
    "layers_names": reactive.Value(None),
    "obsm_names": reactive.Value(None),
    "uns_names": reactive.Value(None),
    "spatial_distance_columns": reactive.Value(None),
    
    # Output dataframes for downloads
    "df_boxplot": reactive.Value(None),
    "df_histogram1": reactive.Value(None),
    "df_histogram2": reactive.Value(None),
    "df_heatmap": reactive.Value(None),
    "df_relational": reactive.Value(None),
    "df_nn": reactive.Value(None),
    "df_ripley": reactive.Value(None),
}
```

### 2. Reactive Value Types

| Type | Purpose | Example Use |
|------|---------|-------------|
| `reactive.Value()` | Stores mutable state | `adata_main = reactive.Value(None)` |
| `@reactive.calc` | Computed values (cached) | Process data when inputs change |
| `@reactive.effect` | Side effects (lowercase) | Update UI elements dynamically |
| `@reactive.Effect` | Side effects (uppercase) | Update UI inputs, same as `@reactive.effect` |
| `@reactive.event()` | Event-triggered reactions | Run only when button clicked |
| `@render.*` | Output rendering | Display plots, tables, text |

---

## Reactive Patterns

### Pattern 1: Data Loading and Decomposition

**Location**: `data_input_server.py`

This pattern handles file uploads and breaks the AnnData object into reactive components.

```python
# Pattern: File Upload → Decompose → Update Shared State

@reactive.Effect
def adata_filter():
    """Load data from file upload or use preloaded data."""
    file_info = input.input_file()
    
    if not file_info:
        # Use preloaded data
        if shared['preloaded_data'] is not None:
            shared['adata_main'].set(shared['preloaded_data'])
            shared['data_loaded'].set(True)
    else:
        # Load uploaded file
        file_path = file_info[0]['datapath']
        with open(file_path, 'rb') as file:
            if file_path.endswith('.pickle'):
                shared['adata_main'].set(pickle.load(file))
            elif file_path.endswith('.h5ad'):
                shared['adata_main'].set(ad.read_h5ad(file_path))
        shared['data_loaded'].set(True)

@reactive.Effect
def update_parts():
    """Decompose adata_main into reactive components."""
    adata = shared['adata_main'].get()
    
    if adata is not None:
        # Extract and store each component
        shared['X_data'].set(adata.X)
        shared['obs_data'].set(adata.obs)
        shared['var_data'].set(adata.var)
        shared['layers_data'].set(adata.layers)
        shared['uns_data'].set(adata.uns)
        shared['obsm_data'].set(adata.obsm)
        
        # Extract metadata
        shared['obs_names'].set(list(adata.obs.keys()))
        shared['var_names'].set(list(adata.var.index.tolist()))
        shared['layers_names'].set(list(adata.layers.keys()))
        # ... etc
```

**Flow**:
1. User uploads file → `adata_filter()` triggers
2. `adata_main` reactive value updated
3. `update_parts()` automatically runs (depends on `adata_main`)
4. All component reactive values updated
5. Downstream effects cascade automatically

---

### Pattern 2: Effect-Based UI Updates

**Location**: `effect_update_server.py`

This pattern watches reactive values and updates UI inputs accordingly.

```python
# Pattern: Watch Reactive Value → Update UI Input

@reactive.Effect
def update_select_input_feat():
    """Update feature selection dropdowns when var_names changes."""
    choices = shared['var_names'].get()
    
    # Update all feature-related dropdowns
    ui.update_select("h1_feat", choices=choices)
    ui.update_select("umap_feat", choices=choices)
    
    if choices is not None:
        ui.update_select("bp_features", choices=choices)
        ui.update_selectize("bp_features", selected=choices[:2])

@reactive.Effect
def update_select_input_anno():
    """Update annotation selection dropdowns when obs_names changes."""
    choices = shared['obs_names'].get()
    
    # Update all annotation-related dropdowns
    ui.update_select("bp_anno", choices=choices)
    ui.update_select("h2_anno", choices=choices)
    ui.update_select("spatial_anno", choices=choices)
    
    if choices is not None and len(choices) > 1:
        ui.update_select("sk1_anno1", choices=choices)
        ui.update_selectize("sk1_anno1", selected=choices[0])
        ui.update_select("sk1_anno2", choices=choices)
        ui.update_selectize("sk1_anno2", selected=choices[1])
```

**Flow**:
1. Data loaded → `obs_names` or `var_names` updated
2. Effects automatically trigger
3. All related UI inputs synchronized
4. User sees updated dropdown options

---

### Pattern 3: Event-Driven Computation

**Location**: Various server modules (e.g., `boxplot_server.py`)

This pattern performs expensive computations only when explicitly triggered.

```python
# Pattern: Button Click → Compute → Render → Store Result

@output
@render_widget
@reactive.event(input.go_bp, ignore_none=True)
def spac_Boxplot():
    """Generate boxplot only when 'Generate' button is clicked."""
    
    if not input.bp_output_type():
        return None
    
    # Reconstruct AnnData from reactive components
    adata = ad.AnnData(
        X=shared['X_data'].get(),
        obs=pd.DataFrame(shared['obs_data'].get()),
        var=pd.DataFrame(shared['var_data'].get()),
        layers=shared['layers_data'].get(),
        dtype=shared['X_data'].get().dtype
    )
    
    # Generate plot
    fig, df = spac.visualization.boxplot_interactive(
        adata,
        annotation=input.bp_anno(),
        features=list(input.bp_features()),
        # ... other parameters
    ).values()
    
    # Store dataframe for download
    shared['df_boxplot'].set(df)
    
    return fig
```

**Flow**:
1. User configures parameters
2. User clicks "Generate" button
3. `spac_Boxplot()` function executes
4. Plot generated and displayed
5. Result dataframe stored for download

---

### Pattern 4: Computed Reactive Values

**Location**: `nearest_neighbor_server.py`, `scatterplot_server.py`

This pattern creates derived values that update automatically but are cached.

```python
# Pattern: Multiple Inputs → Computed Value (Cached)

@reactive.calc
def get_adata():
    """Get the main AnnData object (cached)."""
    return shared['adata_main'].get()

@reactive.calc
def process_target_labels():
    """
    Process target label selection.
    Result is cached until inputs change.
    """
    target_labels = input.nn_target_label()
    if target_labels and len(target_labels) > 0:
        return target_labels
    return None

@reactive.calc
def get_plot_type():
    """Get plot type based on method selection."""
    method = input.nn_plot_method()
    if method == "numeric":
        return input.nn_plot_type_numeric()
    else:
        return input.nn_plot_type_distribution()

# Used in render function
@output
@render.plot
def nn_plot():
    adata = get_adata()  # Reuses cached value
    labels = process_target_labels()  # Reuses cached value
    plot_type = get_plot_type()  # Reuses cached value
    
    # Generate plot using computed values
    # ...
```

**Benefits**:
- Expensive computations run once and are cached
- Automatic invalidation when dependencies change
- Cleaner code organization
- Better performance

---

### Pattern 5: Dynamic UI Generation

**Location**: `spatial_server.py`, `nearest_neighbor_server.py`

This pattern creates UI elements conditionally based on user actions or data.

```python
# Pattern: User Action → Insert/Remove UI Elements

slide_ui_initialized = reactive.Value(False)

@reactive.effect
def slide_reactivity():
    """Insert slide selection UI when checkbox is checked."""
    btn = input.slide_select_check()
    ui_initialized = slide_ui_initialized.get()
    
    if btn and not ui_initialized:
        # Create dropdown
        dropdown_slide = ui.input_select(
            "slide_select_drop",
            "Select the Slide Annotation",
            choices=shared['obs_names'].get()
        )
        
        # Insert into DOM
        ui.insert_ui(
            ui.div({"id": "inserted-slide_dropdown"}, dropdown_slide),
            selector="#main-slide_dropdown",
            where="beforeEnd",
        )
        
        slide_ui_initialized.set(True)
        
    elif not btn and ui_initialized:
        # Remove from DOM
        ui.remove_ui("#inserted-slide_dropdown")
        slide_ui_initialized.set(False)

@reactive.effect
def update_slide_select_drop():
    """Populate slide labels based on selected annotation."""
    adata = ad.AnnData(obs=shared['obs_data'].get())
    
    if input.slide_select_drop():
        selected_anno = input.slide_select_drop()
        labels = adata.obs[selected_anno].unique().tolist()
        ui.update_select("slide_select_label", choices=labels)
```

**Flow**:
1. User checks/unchecks checkbox
2. UI elements dynamically inserted/removed
3. Dependent dropdowns populated from data
4. State tracked to prevent duplicate insertions

---

## Data Flow

### Primary Data Flow Path

```
File Upload
    │
    ▼
data_input_server.adata_filter()
    │
    ├─► shared['adata_main'].set(adata)
    │
    ▼
data_input_server.update_parts()
    │
    ├─► shared['X_data'].set(adata.X)
    ├─► shared['obs_data'].set(adata.obs)
    ├─► shared['var_data'].set(adata.var)
    ├─► shared['obs_names'].set(list(adata.obs.keys()))
    ├─► shared['var_names'].set(list(adata.var.index))
    └─► ... (other components)
    │
    ▼
effect_update_server (multiple effects trigger)
    │
    ├─► update_select_input_feat() → Updates feature dropdowns
    ├─► update_select_input_anno() → Updates annotation dropdowns
    ├─► update_select_input_layer() → Updates layer dropdowns
    └─► update_nearest_neighbor_choices() → Updates NN phenotype choices
    │
    ▼
Feature Servers (ready to use updated data)
    │
    ├─► boxplot_server
    ├─► spatial_server
    ├─► umap_server
    └─► ... (all other feature servers)
```

### Analysis Flow Example (Boxplot)

```
User Input
    │
    ├─► Select features (input.bp_features)
    ├─► Select annotation (input.bp_anno)
    ├─► Select layer (input.bp_layer)
    ├─► Configure options (outliers, orientation, etc.)
    │
    ▼
User clicks "Generate" button (input.go_bp)
    │
    ▼
boxplot_server.spac_Boxplot() triggered (@reactive.event)
    │
    ├─► Retrieve data from shared state:
    │   ├─► X_data
    │   ├─► obs_data
    │   ├─► var_data
    │   └─► layers_data
    │
    ├─► Reconstruct AnnData object
    │
    ├─► Call spac.visualization.boxplot_interactive()
    │
    ├─► Store result: shared['df_boxplot'].set(df)
    │
    └─► Return figure to UI
```

---

## Module Communication

### Inter-Module Communication via Shared State

Modules **never communicate directly**. All communication happens through the `shared` dictionary.

```python
# ❌ BAD: Direct module communication (not possible)
# Module A calling Module B directly - DON'T DO THIS

# ✅ GOOD: Communication via shared state

# Module A writes to shared state
def module_a_server(input, output, session, shared):
    @reactive.Effect
    def compute_something():
        result = expensive_computation()
        shared['computed_result'].set(result)

# Module B reads from shared state
def module_b_server(input, output, session, shared):
    @reactive.calc
    def use_computed_result():
        result = shared['computed_result'].get()
        return transform(result)
```

### Key Communication Patterns

#### 1. **Producer-Consumer Pattern**
- **Producer**: `data_input_server` produces decomposed data
- **Consumer**: All feature servers consume the data

#### 2. **Broadcast Pattern**
- **Broadcaster**: `data_input_server` updates `obs_names`
- **Listeners**: `effect_update_server` updates all relevant UI inputs

#### 3. **Store-Retrieve Pattern**
- **Store**: Feature servers store output dataframes
- **Retrieve**: Download functions retrieve dataframes for CSV export

---

## Developer Guide

### Adding a New Feature Module

Follow these steps to add a new feature (e.g., a new plot type):

#### Step 1: Create UI Module

**File**: `ui/new_feature_ui.py`

```python
"""
New Feature UI Module

This module provides the user interface for the new feature.
"""

from shiny import ui

def new_feature_ui():
    """
    Create the new feature UI.
    
    Returns
    -------
    shiny.ui.TagChild
        UI components for the new feature
    """
    return ui.nav_panel(
        "New Feature",
        ui.layout_sidebar(
            ui.sidebar(
                ui.input_select(
                    "nf_annotation",
                    "Select Annotation:",
                    choices=[]  # Will be populated by effects
                ),
                ui.input_selectize(
                    "nf_features",
                    "Select Features:",
                    choices=[],
                    multiple=True
                ),
                ui.input_action_button(
                    "go_nf",
                    "Generate Plot",
                    class_="btn-primary"
                ),
                width=300
            ),
            ui.output_plot("nf_plot"),
            ui.download_button("download_nf", "Download Data")
        )
    )
```

#### Step 2: Create Server Module

**File**: `server/new_feature_server.py`

```python
"""
New Feature Server Module

This module handles server-side logic for the new feature.
"""

from shiny import render, reactive
import anndata as ad
import pandas as pd

def new_feature_server(input, output, session, shared):
    """
    Server logic for new feature.
    
    Parameters
    ----------
    input : shiny.session.Inputs
        Shiny input object
    output : shiny.session.Outputs
        Shiny output object
    session : shiny.session.Session
        Shiny session object
    shared : dict
        Shared reactive values across server modules
    """
    
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
    
    @output
    @render.plot
    @reactive.event(input.go_nf, ignore_none=True)
    def nf_plot():
        """
        Generate new feature plot.
        
        Returns
        -------
        matplotlib.figure.Figure
            The generated plot figure
        """
        adata = get_adata()
        
        if adata is None:
            return None
        
        # Your plotting logic here
        fig = create_new_plot(
            adata,
            annotation=input.nf_annotation(),
            features=list(input.nf_features())
        )
        
        # Store data for download
        df = extract_plot_data(fig)
        shared['df_new_feature'].set(df)
        
        return fig
    
    @render.download(filename="new_feature_data.csv")
    def download_nf():
        """Provide CSV download of plot data."""
        df = shared['df_new_feature'].get()
        if df is not None:
            return df.to_csv(index=False).encode("utf-8"), "text/csv"
        return None
```

#### Step 3: Add Effect Updates

**File**: `server/effect_update_server.py`

Add a new effect to update your feature's UI inputs:

```python
@reactive.Effect
def update_new_feature_inputs():
    """Update new feature dropdown choices when data changes."""
    obs_choices = shared['obs_names'].get()
    var_choices = shared['var_names'].get()
    
    if obs_choices is not None:
        ui.update_select("nf_annotation", choices=obs_choices)
    
    if var_choices is not None:
        ui.update_selectize("nf_features", choices=var_choices)
```

#### Step 4: Register in Main App

**File**: `app.py`

```python
# Add imports at top
from ui import new_feature_ui
from server import new_feature_server

# In app_ui, add to navset_card_tab:
app_ui = ui.page_fluid(
    # ... existing code
    ui.navset_card_tab(
        getting_started_ui(),
        data_input_ui(),
        # ... other modules
        new_feature_ui(),  # ← Add your UI
    ),
    # ...
)

# In server function, add to data_keys:
data_keys = [
    # ... existing keys
    "df_new_feature",  # ← Add your dataframe key
]

# In server function, register server module:
def server(input, output, session):
    # ... existing setup
    
    # Register your server module
    new_feature_server(input, output, session, shared)
    
    # ... rest of code
```

#### Step 5: Test Reactive Flow

1. **Upload data** → Verify dropdowns populate
2. **Select inputs** → Verify selections work
3. **Click generate** → Verify plot appears
4. **Download data** → Verify CSV downloads

---

### Common Reactive Patterns Checklist

When implementing a new feature, use this checklist:

- [ ] **UI Module Created** (`ui/feature_ui.py`)
  - [ ] Uses semantic HTML elements
  - [ ] Follows accessibility guidelines
  - [ ] Input IDs are unique and descriptive

- [ ] **Server Module Created** (`server/feature_server.py`)
  - [ ] Has NumPy-style docstrings
  - [ ] Uses `@reactive.calc` for expensive computations
  - [ ] Uses `@reactive.event` for button-triggered actions
  - [ ] Reconstructs AnnData from shared state
  - [ ] Stores output dataframes in shared state

- [ ] **Effect Updates Added** (`server/effect_update_server.py`)
  - [ ] Dropdown choices updated when data changes
  - [ ] Default selections set appropriately

- [ ] **Registered in app.py**
  - [ ] UI module imported and added to navset
  - [ ] Server module imported and registered
  - [ ] Output dataframe key added to `data_keys`

- [ ] **Testing Complete**
  - [ ] Data loading triggers updates
  - [ ] UI inputs work correctly
  - [ ] Plot/output generates successfully
  - [ ] Download functionality works

---

## Best Practices

### 1. Reactive Value Management

#### ✅ DO:
```python
# Read reactive value
value = shared['adata_main'].get()

# Write reactive value
shared['adata_main'].set(new_value)

# Use @reactive.calc for derived values
@reactive.calc
def processed_data():
    raw = shared['X_data'].get()
    return process(raw)
```

#### ❌ DON'T:
```python
# Don't access shared values without .get()
value = shared['adata_main']  # Wrong! This is the reactive object

# Don't modify reactive values directly
shared['adata_main'].value = new_value  # Wrong! Use .set()

# Don't compute expensive operations without @reactive.calc
def processed_data():
    # This runs every time, not cached!
    raw = shared['X_data'].get()
    return expensive_process(raw)
```

### 2. Effect Design

#### ✅ DO:
```python
# Keep effects focused and single-purpose
@reactive.Effect
def update_feature_dropdown():
    """Update only feature dropdown."""
    choices = shared['var_names'].get()
    ui.update_select("feature_input", choices=choices)

# Use effects for UI updates
@reactive.Effect
def update_ui_based_on_data():
    """Update UI when data changes."""
    data = shared['adata_main'].get()
    if data is not None:
        ui.update_select("option", choices=data.obs.keys())
```

#### ❌ DON'T:
```python
# Don't mix multiple concerns in one effect
@reactive.Effect
def update_everything():
    # Too much in one effect!
    update_features()
    update_annotations()
    update_layers()
    compute_something()
    update_ui()

# Don't use effects for computations (use @reactive.calc)
@reactive.Effect
def compute_something():
    # Wrong! Effects are for side effects, not computations
    result = expensive_computation()
    return result  # Returns are ignored in effects
```

### 3. Event Handling

#### ✅ DO:
```python
# Use @reactive.event for button clicks
@output
@render.plot
@reactive.event(input.generate_btn, ignore_none=True)
def my_plot():
    """Generate plot only when button clicked."""
    return create_plot()

# Validate inputs before processing
@output
@render.plot
@reactive.event(input.go_btn, ignore_none=True)
def my_plot():
    req(input.feature_select())  # Ensure input exists
    req(shared['adata_main'].get())  # Ensure data loaded
    return create_plot()
```

#### ❌ DON'T:
```python
# Don't make expensive plots reactive to every input change
@output
@render.plot
def my_plot():
    # This reruns on EVERY input change - expensive!
    feature = input.feature_select()
    annotation = input.anno_select()
    layer = input.layer_select()
    return expensive_plot()  # Runs too often!
```

### 4. State Synchronization

#### ✅ DO:
```python
# Store plot data for downloads
@output
@render_widget
@reactive.event(input.go_bp, ignore_none=True)
def plot():
    fig, df = create_plot_with_data()
    shared['df_boxplot'].set(df)  # Store for download
    return fig

# Provide download access
@render.download(filename="data.csv")
def download():
    df = shared['df_boxplot'].get()
    if df is not None:
        return df.to_csv(index=False).encode("utf-8"), "text/csv"
    return None
```

#### ❌ DON'T:
```python
# Don't regenerate data for downloads
@render.download(filename="data.csv")
def download():
    # Wrong! This regenerates the plot just for download
    fig, df = create_plot_with_data()  # Expensive!
    return df.to_csv(index=False).encode("utf-8"), "text/csv"
```

### 5. Module Organization

#### ✅ DO:
```python
# Keep server modules focused
def spatial_server(input, output, session, shared):
    """Only spatial visualization logic here."""
    # ... spatial-specific code

def boxplot_server(input, output, session, shared):
    """Only boxplot logic here."""
    # ... boxplot-specific code

# Use utility functions for shared logic
from utils.data_processing import reconstruct_adata

def my_server(input, output, session, shared):
    adata = reconstruct_adata(shared)  # Reusable utility
```

#### ❌ DON'T:
```python
# Don't mix concerns across modules
def spatial_server(input, output, session, shared):
    # Also handles boxplots - wrong module!
    @output
    @render.plot
    def boxplot():  # This belongs in boxplot_server
        return create_boxplot()
```

### 6. Performance Optimization

#### ✅ DO:
```python
# Use @reactive.calc for expensive operations
@reactive.calc
def expensive_computation():
    """This result is cached until dependencies change."""
    data = shared['X_data'].get()
    return costly_transform(data)

# Reuse computed values
@output
@render.plot
def plot_a():
    result = expensive_computation()  # Cached result
    return plot(result)

@output
@render.plot
def plot_b():
    result = expensive_computation()  # Same cached result!
    return different_plot(result)
```

#### ❌ DON'T:
```python
# Don't recompute unnecessarily
@output
@render.plot
def plot_a():
    data = shared['X_data'].get()
    result = costly_transform(data)  # Computed again!
    return plot(result)

@output
@render.plot
def plot_b():
    data = shared['X_data'].get()
    result = costly_transform(data)  # Computed again!
    return different_plot(result)
```

---

## Debugging Reactive Applications

### Common Issues and Solutions

#### Issue 1: Dropdown Not Updating

**Symptom**: After loading data, dropdown stays empty

**Causes**:
- Effect not properly watching reactive value
- Effect not registered in `effect_update_server.py`
- Reactive value not being set in `data_input_server.py`

**Solution**:
```python
# 1. Verify data is being set
@reactive.Effect
def update_parts():
    adata = shared['adata_main'].get()
    if adata is not None:
        shared['obs_names'].set(list(adata.obs.keys()))
        print(f"Set obs_names: {shared['obs_names'].get()}")  # Debug

# 2. Verify effect is triggering
@reactive.Effect
def update_dropdown():
    choices = shared['obs_names'].get()
    print(f"Updating dropdown with: {choices}")  # Debug
    ui.update_select("my_dropdown", choices=choices)
```

#### Issue 2: Infinite Reactive Loop

**Symptom**: Application freezes or becomes unresponsive

**Causes**:
- Effect reads and writes same reactive value
- Circular dependency between reactive values

**Solution**:
```python
# ❌ BAD: Creates infinite loop
@reactive.Effect
def bad_effect():
    value = shared['counter'].get()
    shared['counter'].set(value + 1)  # Triggers itself!

# ✅ GOOD: Use event to break loop
@reactive.effect
@reactive.event(input.increment_btn)
def good_effect():
    value = shared['counter'].get()
    shared['counter'].set(value + 1)  # Only on button click
```

#### Issue 3: Render Function Not Updating

**Symptom**: Plot doesn't update when data changes

**Causes**:
- Using `@reactive.event` blocks automatic reactivity
- Not reading reactive values inside render function

**Solution**:
```python
# ❌ BAD: Reactive value read outside render function
adata = shared['adata_main'].get()  # Read once at module load

@output
@render.plot
def my_plot():
    return plot(adata)  # Uses stale data!

# ✅ GOOD: Read reactive value inside render function
@output
@render.plot
def my_plot():
    adata = shared['adata_main'].get()  # Read on each render
    return plot(adata)  # Uses fresh data
```

---

## Advanced Topics

### Conditional Reactivity

Use `@reactive.event` with `ignore_none=True` to prevent execution when inputs are None:

```python
@output
@render.plot
@reactive.event(input.go_btn, ignore_none=True)
def my_plot():
    """Only runs when button clicked AND button is not None."""
    return create_plot()
```

### Multiple Event Triggers

React to multiple events:

```python
@reactive.effect
@reactive.event(input.btn1, input.btn2)
def handle_either_button():
    """Runs when either btn1 OR btn2 is clicked."""
    if input.btn1():
        do_something()
    elif input.btn2():
        do_something_else()
```

### Isolating Reactive Dependencies

Use `with reactive.isolate()` to read values without creating dependencies:

```python
@reactive.calc
def computed_value():
    # Creates dependency on value_a
    a = shared['value_a'].get()
    
    # Does NOT create dependency on value_b
    with reactive.isolate():
        b = shared['value_b'].get()
    
    return a + b  # Only re-runs when value_a changes
```

---

## Summary

The SPAC Shiny app uses a **centralized reactive state architecture** where:

1. **Shared State**: All modules communicate through a `shared` dictionary of reactive values
2. **Data Flow**: `data_input_server` → `effect_update_server` → Feature servers
3. **Modularity**: Each feature has separate UI and server modules
4. **Reactivity**: Changes propagate automatically through the dependency graph
5. **Event-Driven**: Expensive operations triggered explicitly by user actions

This architecture provides:
- ✅ Clear separation of concerns
- ✅ Automatic UI synchronization
- ✅ Efficient computation caching
- ✅ Easy addition of new features
- ✅ Maintainable codebase

For new contributors, start by reading the [Developer Guide](#developer-guide) section and follow the step-by-step instructions to add new features.

---

## Additional Resources

- [Shiny for Python Documentation](https://shiny.rstudio.com/py/)
- [Reactive Programming Guide](https://shiny.rstudio.com/py/docs/reactive-programming.html)
- [Project Copilot Instructions](.github/copilot-instructions.md)
- [Contributing Guidelines](SCSAWorkflow/CONTRIBUTING.md)

---

**Version**: 1.0  
**Last Updated**: November 21, 2025  
**Maintained By**: SPAC Development Team
