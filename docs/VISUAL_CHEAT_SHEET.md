# SPAC Shiny - Reactive Framework Visual Cheat Sheet

**Print this page for quick reference while developing!**

---

## 📦 Shared State Dictionary Structure

```
shared = {
    ┌─────────────────────────────────────────────────┐
    │ MAIN DATA                                       │
    ├─────────────────────────────────────────────────┤
    │ adata_main          → Full AnnData object       │
    │ preloaded_data      → Demo data on startup      │
    │ data_loaded         → Boolean flag              │
    └─────────────────────────────────────────────────┘
    
    ┌─────────────────────────────────────────────────┐
    │ DATA COMPONENTS                                 │
    ├─────────────────────────────────────────────────┤
    │ X_data              → Expression matrix         │
    │ obs_data            → Cell annotations          │
    │ var_data            → Feature information       │
    │ layers_data         → Additional data layers    │
    │ uns_data            → Unstructured data         │
    │ obsm_data           → Multi-dim annotations     │
    │ shape_data          → Data dimensions           │
    └─────────────────────────────────────────────────┘
    
    ┌─────────────────────────────────────────────────┐
    │ METADATA (Column Names)                         │
    ├─────────────────────────────────────────────────┤
    │ obs_names           → ["cell_type", "region"]   │
    │ var_names           → ["CD3", "CD8", "CD20"]    │
    │ layers_names        → ["raw", "normalized"]     │
    │ obsm_names          → ["X_umap", "X_pca"]       │
    │ uns_names           → ["neighbors", "colors"]   │
    │ spatial_distance_   → NN analysis columns       │
    │   columns                                       │
    └─────────────────────────────────────────────────┘
    
    ┌─────────────────────────────────────────────────┐
    │ OUTPUT DATAFRAMES (for downloads)               │
    ├─────────────────────────────────────────────────┤
    │ df_boxplot          → Boxplot data              │
    │ df_histogram1       → Feature histogram data    │
    │ df_histogram2       → Annotation histogram data │
    │ df_heatmap          → Heatmap data              │
    │ df_relational       → Relational plot data      │
    │ df_nn               → Nearest neighbor data     │
    │ df_ripley           → Ripley's L data           │
    └─────────────────────────────────────────────────┘
}
```

---

## 🎯 Reactive Decorators Quick Reference

| Decorator | Use For | Returns | Caches | Example |
|-----------|---------|---------|--------|---------|
| `@reactive.Value()` | Store mutable state | N/A | N/A | `x = reactive.Value(0)` |
| `@reactive.calc` | Compute derived values | ✅ Yes | ✅ Yes | Process data |
| `@reactive.effect` | Side effects (UI updates) | ❌ No | ❌ No | Update dropdowns |
| `@reactive.event()` | Wait for specific trigger | N/A | N/A | Button clicks only |
| `@render.plot` | Display matplotlib plot | ✅ Yes | ❌ No | Show visualization |
| `@render_widget` | Display plotly widget | ✅ Yes | ❌ No | Interactive plot |
| `@render.text` | Display text output | ✅ Yes | ❌ No | Show strings |
| `@render.ui` | Dynamic UI elements | ✅ Yes | ❌ No | Conditional inputs |
| `@render.download()` | File downloads | ✅ Yes | ❌ No | CSV export |

---

## 🔄 Data Flow Diagram (Simplified)

```
┌──────────────┐
│ User Uploads │
│   File       │
└──────┬───────┘
       │
       ▼
┌──────────────────────┐
│ data_input_server    │ ◄─── Loads data from disk
│ • adata_filter()     │
│ • update_parts()     │
└──────┬───────────────┘
       │ Sets reactive values
       ▼
┌─────────────────────────────┐
│ Shared State Dictionary     │
│ • adata_main                │
│ • X_data, obs_data, etc.    │
│ • obs_names, var_names      │
└──────┬──────────────────────┘
       │ Changes detected
       ▼
┌──────────────────────────────┐
│ effect_update_server         │ ◄─── Syncs UI with data
│ • update_select_input_feat() │
│ • update_select_input_anno() │
│ • update_select_input_layer()│
└──────┬───────────────────────┘
       │ Updates UI
       ▼
┌──────────────────┐
│ All UI Dropdowns │ ◄─── User sees options
│ Auto-Populated   │
└──────────────────┘
       │
       │ User makes selections
       ▼
┌──────────────────────────┐
│ Feature Server Modules   │ ◄─── Generate outputs
│ • boxplot_server         │
│ • spatial_server         │
│ • umap_server            │
└──────┬───────────────────┘
       │ Stores results
       ▼
┌─────────────────────┐
│ df_* in Shared      │ ◄─── For downloads
│ State               │
└─────────────────────┘
       │
       ▼
┌─────────────────────┐
│ User Downloads CSV  │
└─────────────────────┘
```

---

## ✅ DO's and ❌ DON'Ts

### Reading/Writing Shared State

```python
# ✅ DO
value = shared['adata_main'].get()
shared['adata_main'].set(new_value)

# ❌ DON'T
value = shared['adata_main']           # Missing .get()
shared['adata_main'].value = x        # Wrong method
shared['adata_main'] = x               # Replaces reactive object!
```

### Using @reactive.calc for Performance

```python
# ✅ DO - Cached computation
@reactive.calc
def expensive_data():
    raw = shared['X_data'].get()
    return transform(raw)  # Computed once!

# ❌ DON'T - Recomputes every time
def expensive_data():
    raw = shared['X_data'].get()
    return transform(raw)  # No caching!
```

### Button-Triggered Actions

```python
# ✅ DO - Explicit trigger
@output
@render.plot
@reactive.event(input.go_btn, ignore_none=True)
def my_plot():
    return expensive_plot()  # Only on button click

# ❌ DON'T - Runs on every input change
@output
@render.plot
def my_plot():
    return expensive_plot()  # Too expensive!
```

### UI Updates

```python
# ✅ DO - Use effects
@reactive.Effect
def update_dropdown():
    choices = shared['obs_names'].get()
    ui.update_select("my_select", choices=choices)

# ❌ DON'T - Compute in effect
@reactive.Effect
def bad_effect():
    result = expensive_computation()
    return result  # Effects ignore returns!
```

---

## 🔧 Common Code Patterns

### Pattern: Reconstruct AnnData

```python
@reactive.calc
def get_adata():
    return ad.AnnData(
        X=shared['X_data'].get(),
        obs=pd.DataFrame(shared['obs_data'].get()),
        var=pd.DataFrame(shared['var_data'].get()),
        layers=shared['layers_data'].get(),
        dtype=shared['X_data'].get().dtype
    )
```

### Pattern: Button-Triggered Plot

```python
@output
@render.plot
@reactive.event(input.go_btn, ignore_none=True)
def my_plot():
    adata = get_adata()
    fig = create_plot(adata, ...)
    shared['df_output'].set(extract_data(fig))
    return fig
```

### Pattern: CSV Download

```python
@render.download(filename="data.csv")
def download():
    df = shared['df_output'].get()
    if df is not None:
        return df.to_csv(index=False).encode("utf-8"), "text/csv"
    return None
```

### Pattern: Update Dropdown Choices

```python
@reactive.Effect
def update_choices():
    choices = shared['obs_names'].get()
    ui.update_select("my_input", choices=choices)
    if choices and len(choices) > 0:
        ui.update_select("my_input", selected=choices[0])
```

### Pattern: Dynamic UI Insertion

```python
ui_inserted = reactive.Value(False)

@reactive.effect
def toggle_ui():
    show = input.show_checkbox()
    inserted = ui_inserted.get()
    
    if show and not inserted:
        ui.insert_ui(
            ui.div({"id": "my-ui"}, ui.input_select(...)),
            selector="#container",
            where="beforeEnd"
        )
        ui_inserted.set(True)
    elif not show and inserted:
        ui.remove_ui("#my-ui")
        ui_inserted.set(False)
```

### Pattern: Input Validation

```python
from shiny import req

@output
@render.plot
@reactive.event(input.go_btn, ignore_none=True)
def my_plot():
    req(input.feature_select())        # Stop if None
    req(shared['adata_main'].get())    # Stop if None
    # Safe to proceed
    return create_plot()
```

---

## 🐛 Common Errors & Solutions

| Error | Cause | Solution |
|-------|-------|----------|
| `AttributeError: 'NoneType'` | Accessing None value | Use `req()` to validate |
| Dropdown not updating | Effect not watching value | Check reactive dependency |
| Infinite loop / freeze | Effect reads & writes same value | Use `@reactive.event` |
| Plot not refreshing | Reading value outside render | Read inside render function |
| Expensive recomputation | Not using `@reactive.calc` | Add `@reactive.calc` |
| UI inserted multiple times | Not tracking state | Use `reactive.Value()` flag |
| Download has stale data | Regenerating instead of using stored | Store in shared, retrieve for download |

---

## 📝 Module Template Checklist

When adding a new feature:

- [ ] **UI Module** (`ui/feature_ui.py`)
  - [ ] Uses `ui.nav_panel()`
  - [ ] Input IDs are unique and descriptive
  - [ ] Includes action button for generation
  - [ ] Has output placeholder

- [ ] **Server Module** (`server/feature_server.py`)
  - [ ] Has NumPy-style docstrings
  - [ ] Uses `@reactive.calc` for expensive ops
  - [ ] Uses `@reactive.event` for button triggers
  - [ ] Reconstructs AnnData from shared state
  - [ ] Stores output in `shared['df_*']`

- [ ] **Effect Update** (`server/effect_update_server.py`)
  - [ ] Adds `@reactive.Effect` for UI syncs
  - [ ] Updates relevant dropdown choices
  - [ ] Sets sensible defaults

- [ ] **Registration** (`app.py`)
  - [ ] UI module imported
  - [ ] Server module imported
  - [ ] UI added to `navset_card_tab`
  - [ ] Server registered with shared state
  - [ ] Output dataframe key in `data_keys`

---

## 📊 Performance Tips

| Tip | Why | Example |
|-----|-----|---------|
| Use `@reactive.calc` | Caching prevents recomputation | Expensive transforms |
| Limit effect scope | Fewer UI updates | One effect per concern |
| Use `@reactive.event` | Prevents constant re-rendering | Wait for button click |
| Store outputs | Avoid regenerating for downloads | `shared['df_*'].set(df)` |
| Batch UI updates | Fewer DOM manipulations | One effect for multiple inputs |

---

## 🔗 Quick Links

- **Full Guide**: [REACTIVE_ARCHITECTURE.md](REACTIVE_ARCHITECTURE.md)
- **Quick Ref**: [REACTIVE_PATTERNS_QUICK_REF.md](REACTIVE_PATTERNS_QUICK_REF.md)
- **Demo**: [DEMO_PRESENTATION.md](DEMO_PRESENTATION.md)
- **Diagrams**: [architecture_diagram.mmd](architecture_diagram.mmd), [data_flow_diagrams.mmd](data_flow_diagrams.mmd)

---

## 🎓 Learning Steps

1. **Understand shared state** - Central dictionary with reactive values
2. **Learn decorators** - `@reactive.calc`, `@reactive.effect`, `@reactive.event`
3. **Practice patterns** - Read/write state, update UI, trigger actions
4. **Add a feature** - Follow module template checklist
5. **Debug issues** - Use error table and validation tips

---

**Keep this handy while coding!** 📌

*Version 1.0 - Last Updated: November 21, 2025*
