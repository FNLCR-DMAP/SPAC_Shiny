# SPAC Shiny - Reactive Framework Demo

**Presentation Summary for Colleagues & New Contributors**

---

## What is Reactive Programming?

**Traditional Programming**:
```python
# You manually update everything
data = load_data()
process_data(data)
update_ui(data)
# If data changes, you must remember to update everything again
```

**Reactive Programming**:
```python
# You define relationships once
data = reactive.Value(load_data())
# Everything updates automatically when data changes!
```

---

## SPAC App Architecture Overview

### The Big Picture

```
User Uploads File
    ↓
Data Input Server (loads & decomposes data)
    ↓
Shared Reactive State (central store)
    ↓
Effect Update Server (synchronizes UI)
    ↓
Feature Servers (generate visualizations)
    ↓
User sees results
```

### Key Innovation: **Shared State Pattern**

All modules communicate through a central `shared` dictionary containing reactive values.

**Benefits**:
- ✅ No module coupling - easy to add features
- ✅ Automatic updates - change data once, UI updates everywhere
- ✅ Clear data flow - easy to debug
- ✅ Performance - smart caching prevents waste

---

## The Shared State Dictionary

Think of it as a "reactive database" that all modules can access:

```python
shared = {
    # Main data
    "adata_main": reactive.Value(None),
    
    # Decomposed components
    "X_data": reactive.Value(None),      # Expression matrix
    "obs_data": reactive.Value(None),    # Cell annotations
    "var_data": reactive.Value(None),    # Feature info
    
    # Metadata
    "obs_names": reactive.Value(None),   # Annotation columns
    "var_names": reactive.Value(None),   # Feature names
    
    # Output data for downloads
    "df_boxplot": reactive.Value(None),
    "df_nn": reactive.Value(None),
    # ... etc
}
```

---

## Core Reactive Patterns

### Pattern 1: Auto-Update When Data Changes

```python
# In effect_update_server.py
@reactive.Effect
def update_all_feature_dropdowns():
    """Automatically runs when var_names changes"""
    choices = shared['var_names'].get()  # Subscribe to changes
    
    # Update ALL feature dropdowns at once
    ui.update_select("boxplot_features", choices=choices)
    ui.update_select("umap_features", choices=choices)
    ui.update_select("spatial_features", choices=choices)
```

**Result**: Upload data once → All 12 feature dropdowns update automatically!

---

### Pattern 2: Button-Triggered Analysis

```python
# In boxplot_server.py
@output
@render_widget
@reactive.event(input.go_btn, ignore_none=True)  # Wait for button!
def create_boxplot():
    """Only runs when user clicks 'Generate' button"""
    
    # Get data from shared state
    adata = reconstruct_adata_from_shared(shared)
    
    # Generate plot
    fig, df = spac.visualization.boxplot_interactive(adata, ...)
    
    # Store data for download
    shared['df_boxplot'].set(df)
    
    return fig  # Display to user
```

**Result**: Expensive computations only run when needed!

---

### Pattern 3: Cached Computations

```python
@reactive.calc  # Magic: result is cached!
def expensive_transform():
    """Runs once, result reused everywhere"""
    data = shared['X_data'].get()
    return costly_computation(data)  # Only computed once

@output
@render.plot
def plot1():
    result = expensive_transform()  # Uses cache
    return plot_a(result)

@output  
@render.plot
def plot2():
    result = expensive_transform()  # Same cache, no recompute!
    return plot_b(result)
```

**Result**: 10x faster rendering when multiple plots use same data!

---

## Real Example: What Happens When You Upload Data?

**Step-by-step flow**:

1. **User uploads** `experiment.h5ad` file

2. **data_input_server** activates:
   ```python
   shared['adata_main'].set(load(file))  # One line!
   ```

3. **Automatic decomposition** happens:
   ```python
   shared['X_data'].set(adata.X)
   shared['obs_data'].set(adata.obs)
   shared['var_data'].set(adata.var)
   shared['obs_names'].set(list(adata.obs.keys()))
   # ... etc
   ```

4. **effect_update_server** detects changes and **automatically updates**:
   - 15+ feature dropdowns
   - 20+ annotation dropdowns
   - 5+ layer dropdowns
   - Nearest neighbor phenotype choices
   
   **All with just 4 effect functions!**

5. **Feature servers** now have access to data via `shared.get()`

6. **User sees** all dropdowns populated, ready to analyze

**Total lines of code for this magic: ~50 lines**  
**Without reactive: Would need 200+ lines of manual UI updates**

---

## Adding a New Feature: 3 Simple Steps

### Step 1: Create UI Module

```python
# ui/my_feature_ui.py
def my_feature_ui():
    return ui.nav_panel(
        "My Feature",
        ui.input_select("mf_annotation", "Select:", choices=[]),
        ui.input_action_button("go_mf", "Generate"),
        ui.output_plot("mf_plot")
    )
```

### Step 2: Create Server Module

```python
# server/my_feature_server.py
def my_feature_server(input, output, session, shared):
    @output
    @render.plot
    @reactive.event(input.go_mf, ignore_none=True)
    def mf_plot():
        adata = get_from_shared(shared)
        return create_plot(adata)
```

### Step 3: Add Effect (in effect_update_server.py)

```python
@reactive.Effect
def update_my_feature():
    choices = shared['obs_names'].get()
    ui.update_select("mf_annotation", choices=choices)
```

**Done!** New feature integrated with full reactive support.

---

## Key Benefits for Our Team

### 1. **Rapid Development**
- Add new plots in ~30 minutes
- No need to manually wire up data flow
- Templates make it consistent

### 2. **Easy Debugging**
- Data flow is explicit and traceable
- One place to check: shared state
- Print statements show reactive chain

### 3. **Great Performance**
- `@reactive.calc` caches expensive operations
- Only recompute when dependencies change
- UI updates are batched automatically

### 4. **Maintainable Code**
- Modules are independent
- Clear separation of concerns
- Easy to test individual components

### 5. **Scalable Architecture**
- Currently 12 feature modules
- Can easily add 20 more
- No spaghetti code!

---

## Common Questions

### Q: What if I need data from another module?

**A**: Don't call modules directly! Use shared state:

```python
# ❌ DON'T: Try to call other modules
result = spatial_server.get_data()  # Doesn't work!

# ✅ DO: Use shared state
shared['spatial_result'].set(data)  # One module writes
result = shared['spatial_result'].get()  # Another reads
```

### Q: Why are there two decorators (Effect vs effect)?

**A**: They're identical! Use either:
- `@reactive.Effect` (uppercase) - Traditional style
- `@reactive.effect` (lowercase) - New style

Both work the same way.

### Q: When should I use @reactive.calc vs @reactive.effect?

**A**: 
- **`@reactive.calc`**: For **computing values** (returns something)
  ```python
  @reactive.calc
  def compute_score():
      return calculate(data)
  ```

- **`@reactive.effect`**: For **side effects** (returns nothing)
  ```python
  @reactive.effect
  def update_ui():
      ui.update_select(...)  # No return
  ```

### Q: How do I prevent infinite loops?

**A**: Use `@reactive.event` to break the loop:

```python
# ❌ BAD: Infinite loop
@reactive.effect
def bad():
    x = shared['counter'].get()
    shared['counter'].set(x + 1)  # Triggers itself forever!

# ✅ GOOD: Only on button click
@reactive.effect
@reactive.event(input.increment_btn)
def good():
    x = shared['counter'].get()
    shared['counter'].set(x + 1)  # Only when button clicked
```

---

## Demo: Live Code Walkthrough

### Example 1: See Auto-Update in Action

1. Open `app.py` and find `shared` dictionary
2. Open `data_input_server.py` - see how data is loaded
3. Open `effect_update_server.py` - see how UI updates
4. Run app and upload data
5. Watch console logs show reactive chain executing

### Example 2: Add a Simple Feature

Let's add a "Data Summary" feature that shows basic stats:

1. Create `ui/summary_ui.py`
2. Create `server/summary_server.py`
3. Add effect in `effect_update_server.py`
4. Register in `app.py`
5. Test it!

**Time estimate**: 15 minutes

---

## Resources & Documentation

### 📚 Full Documentation

- **[REACTIVE_ARCHITECTURE.md](REACTIVE_ARCHITECTURE.md)** - Comprehensive guide (50 pages)
- **[docs/REACTIVE_PATTERNS_QUICK_REF.md](docs/REACTIVE_PATTERNS_QUICK_REF.md)** - Quick reference
- **[docs/architecture_diagram.mmd](docs/architecture_diagram.mmd)** - Visual architecture
- **[docs/data_flow_diagrams.mmd](docs/data_flow_diagrams.mmd)** - Sequence diagrams

### 🎓 Learning Path

1. **Day 1**: Read REACTIVE_ARCHITECTURE.md overview & core patterns
2. **Day 2**: Follow developer guide, add a simple feature
3. **Day 3**: Study existing modules, add a plot module
4. **Ongoing**: Use quick reference for daily development

### 🔗 External Resources

- [Shiny for Python Docs](https://shiny.rstudio.com/py/)
- [Reactive Programming Guide](https://shiny.rstudio.com/py/docs/reactive-programming.html)

---

## Summary: Why This Architecture Rocks

| Traditional Approach | Our Reactive Approach |
|---------------------|----------------------|
| Manual UI updates | Automatic updates |
| Tight coupling between modules | Loose coupling via shared state |
| Hard to add features | Easy plugin architecture |
| Duplicate code everywhere | Reusable patterns |
| Difficult debugging | Clear data flow |
| Poor performance | Optimized with caching |

---

## Next Steps

### For Demo Viewers
1. ✅ Review this document
2. ✅ Look at architecture diagram
3. ✅ Explore one module in detail
4. ✅ Try adding a simple feature

### For New Contributors
1. ✅ Read full REACTIVE_ARCHITECTURE.md
2. ✅ Set up development environment
3. ✅ Follow developer guide step-by-step
4. ✅ Submit your first feature PR!

### For Existing Team
1. ✅ Refactor old code to use patterns
2. ✅ Add missing docstrings
3. ✅ Optimize with @reactive.calc
4. ✅ Share knowledge with team

---

## Questions?

**Contact**: SPAC Development Team  
**Documentation**: See `docs/` folder  
**Contributing**: See `CONTRIBUTING.md`

---

**Thank you for your interest in SPAC Shiny!** 🎉

We hope this reactive architecture makes your development experience smooth and enjoyable.

---

*Last Updated: November 21, 2025*  
*Version: 1.0*
