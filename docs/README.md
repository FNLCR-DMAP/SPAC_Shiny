# SPAC Shiny App - Documentation

This directory contains comprehensive documentation for understanding and contributing to the SPAC Shiny for Python application.

## 📚 Documentation Index

### For All Users

- **[REACTIVE_ARCHITECTURE.md](REACTIVE_ARCHITECTURE.md)** - Complete guide to the app's reactive architecture
  - Overview of reactive programming concepts
  - Detailed architecture diagrams
  - Core reactive components explanation
  - Data flow documentation
  - Developer guide for adding features
  - Best practices and anti-patterns

### For Quick Reference

- **[REACTIVE_PATTERNS_QUICK_REF.md](REACTIVE_PATTERNS_QUICK_REF.md)** - Condensed reference guide
  - Decorator cheat sheet
  - Common code patterns
  - Shared state API reference
  - Quick troubleshooting guide
  - Minimal module templates

### Visual Diagrams

- **[architecture_diagram.mmd](architecture_diagram.mmd)** - High-level system architecture
  - Component relationships
  - Data flow paths
  - Module interactions
  - State management visualization

- **[data_flow_diagrams.mmd](data_flow_diagrams.mmd)** - Sequence diagrams for key operations
  - Data loading flow
  - Plot generation flow
  - Download flow
  - Dynamic UI insertion
  - Cached computation flow
  - Effect chain flow
  - Error prevention flow

## 🎯 Where to Start

### **For Demo/Presentation:**
1. Start with [architecture_diagram.mmd](architecture_diagram.mmd) - render this as an image
2. Review key sections in [REACTIVE_ARCHITECTURE.md](../REACTIVE_ARCHITECTURE.md):
   - Overview
   - Architecture Diagram (text-based)
   - Core Reactive Components
   - Reactive Patterns
3. Use [data_flow_diagrams.mmd](data_flow_diagrams.mmd) to show specific flows

### **For New Contributors:**
1. Read [REACTIVE_ARCHITECTURE.md](../REACTIVE_ARCHITECTURE.md) thoroughly
2. Keep [REACTIVE_PATTERNS_QUICK_REF.md](REACTIVE_PATTERNS_QUICK_REF.md) open while coding
3. Follow the "Developer Guide" section step-by-step
4. Refer to [data_flow_diagrams.mmd](data_flow_diagrams.mmd) to understand specific operations

### **For Experienced Developers:**
1. Use [REACTIVE_PATTERNS_QUICK_REF.md](REACTIVE_PATTERNS_QUICK_REF.md) as your main reference
2. Consult [REACTIVE_ARCHITECTURE.md](../REACTIVE_ARCHITECTURE.md) for architectural decisions
3. Review [data_flow_diagrams.mmd](data_flow_diagrams.mmd) when debugging complex flows

## 🖼️ Viewing Mermaid Diagrams

Mermaid diagrams (`.mmd` files) can be viewed in several ways:

### Option 1: GitHub/GitLab
- Push to GitHub/GitLab and view directly in the web interface
- Both platforms render Mermaid diagrams automatically

### Option 2: VS Code
1. Install "Markdown Preview Mermaid Support" extension
2. Open the `.mmd` file
3. Use preview: `Cmd+Shift+V` (Mac) or `Ctrl+Shift+V` (Windows/Linux)

### Option 3: Online Editors
- [Mermaid Live Editor](https://mermaid.live/) - Official editor
- Copy/paste diagram code to view and export as PNG/SVG

### Option 4: Command Line (Export as Image)
```bash
# Install mermaid-cli
npm install -g @mermaid-js/mermaid-cli

# Generate PNG
mmdc -i architecture_diagram.mmd -o architecture_diagram.png

# Generate SVG (better quality)
mmdc -i architecture_diagram.mmd -o architecture_diagram.svg
```

## 📖 Documentation Structure

```
docs/
├── README.md                          (This file)
├── REACTIVE_ARCHITECTURE.md           (Main comprehensive guide)
├── QUICK_START_DOCS.md                (Quick start guide)
├── REACTIVE_PATTERNS_QUICK_REF.md     (Quick reference guide)
├── DEMO_PRESENTATION.md               (Demo/presentation material)
├── VISUAL_CHEAT_SHEET.md              (Printable cheat sheet)
├── DOCUMENTATION_SUMMARY.md           (What was created)
├── architecture_diagram.mmd           (System architecture diagram)
├── data_flow_diagrams.mmd             (Sequence diagrams)
└── documentation_navigation_map.mmd   (Documentation map)

../
├── .github/
│   └── copilot-instructions.md        (Development guidelines)
├── CONTRIBUTING.md                    (Contribution guidelines)
└── README.md                          (Project overview)
```

## 🔑 Key Concepts Summary

### Reactive Programming in SPAC
The app uses **Shiny for Python's reactive framework** where:
- Changes to reactive values automatically trigger dependent computations
- UI updates happen automatically when data changes
- Expensive operations are cached and reused
- Button clicks can explicitly trigger operations

### Shared State Pattern
- Central `shared` dictionary contains all reactive values
- Modules communicate only through shared state (no direct calls)
- Data flows: Input → State → Effects → UI → Feature Servers → Outputs

### Module Organization
- **UI modules** (`ui/*.py`): Define user interface components
- **Server modules** (`server/*.py`): Implement business logic
- **Effect updates** (`server/effect_update_server.py`): Synchronize UI
- **Utilities** (`utils/*.py`): Reusable helper functions

## 🎓 Learning Path

### Beginner
1. **Read**: [REACTIVE_ARCHITECTURE.md](REACTIVE_ARCHITECTURE.md) - Overview & Core Components sections
2. **Study**: Simple patterns (Pattern 1 & 2) in the same document
3. **View**: [data_flow_diagrams.mmd](data_flow_diagrams.mmd) - Diagram 1 & 2
4. **Practice**: Add a simple UI element that displays data from shared state

### Intermediate
1. **Read**: [REACTIVE_ARCHITECTURE.md](REACTIVE_ARCHITECTURE.md) - All reactive patterns
2. **Study**: Complete Developer Guide section
3. **Reference**: [REACTIVE_PATTERNS_QUICK_REF.md](REACTIVE_PATTERNS_QUICK_REF.md)
4. **Practice**: Add a new feature module with plot generation

### Advanced
1. **Master**: All patterns in [REACTIVE_ARCHITECTURE.md](REACTIVE_ARCHITECTURE.md)
2. **Study**: [data_flow_diagrams.mmd](data_flow_diagrams.mmd) - All diagrams
3. **Optimize**: Use cached computations and advanced reactive patterns
4. **Practice**: Refactor existing modules, optimize performance

## 🛠️ Development Workflow

When adding a new feature:

1. **Plan** - Review Developer Guide in [REACTIVE_ARCHITECTURE.md](../REACTIVE_ARCHITECTURE.md)
2. **Reference** - Keep [REACTIVE_PATTERNS_QUICK_REF.md](REACTIVE_PATTERNS_QUICK_REF.md) open
3. **Code** - Follow the step-by-step template
4. **Test** - Verify reactive flow works correctly
5. **Document** - Add comments using NumPy-style docstrings
6. **Review** - Check against best practices checklist

## 📝 Example Use Cases

### Use Case 1: Understanding Data Flow
**Goal**: Understand what happens when user uploads a file

**Documents to read**:
1. [REACTIVE_ARCHITECTURE.md](REACTIVE_ARCHITECTURE.md) - "Data Flow" section
2. [data_flow_diagrams.mmd](data_flow_diagrams.mmd) - Diagram 1: "Data Loading Flow"

### Use Case 2: Adding a New Plot Type
**Goal**: Add a new visualization module

**Documents to follow**:
1. [REACTIVE_ARCHITECTURE.md](REACTIVE_ARCHITECTURE.md) - "Developer Guide" section
2. [REACTIVE_PATTERNS_QUICK_REF.md](REACTIVE_PATTERNS_QUICK_REF.md) - "Quick Module Template"
3. [data_flow_diagrams.mmd](data_flow_diagrams.mmd) - Diagram 2: "Boxplot Generation Flow" (as reference)

### Use Case 3: Debugging Dropdown Not Updating
**Goal**: Fix UI synchronization issue

**Documents to consult**:
1. [REACTIVE_PATTERNS_QUICK_REF.md](REACTIVE_PATTERNS_QUICK_REF.md) - "Quick Troubleshooting"
2. [REACTIVE_ARCHITECTURE.md](REACTIVE_ARCHITECTURE.md) - "Pattern 2: Effect-Based UI Updates"
3. [data_flow_diagrams.mmd](data_flow_diagrams.mmd) - Diagram 6: "Effect Chain Flow"

### Use Case 4: Optimizing Performance
**Goal**: Reduce unnecessary recomputations

**Documents to study**:
1. [REACTIVE_PATTERNS_QUICK_REF.md](REACTIVE_PATTERNS_QUICK_REF.md) - "Pattern 4: Cached Computation"
2. [REACTIVE_ARCHITECTURE.md](REACTIVE_ARCHITECTURE.md) - "Best Practices" - Performance section
3. [data_flow_diagrams.mmd](data_flow_diagrams.mmd) - Diagram 5: "Cached Computation Flow"

## 🤝 Contributing to Documentation

If you find errors or want to improve the documentation:

1. Check [CONTRIBUTING.md](../SCSAWorkflow/CONTRIBUTING.md) for contribution guidelines
2. Follow NumPy-style docstrings for code documentation
3. Use conventional commit messages (see [.github/copilot-instructions.md](../.github/copilot-instructions.md))
4. Test that Mermaid diagrams render correctly before committing

## 📞 Questions?

- For general questions: See [README.md](../README.md)
- For development guidelines: See [.github/copilot-instructions.md](../.github/copilot-instructions.md)
- For contribution process: See [CONTRIBUTING.md](../SCSAWorkflow/CONTRIBUTING.md)

## 📅 Documentation Versioning

- **Current Version**: 1.0
- **Last Updated**: November 21, 2025
- **Maintained By**: SPAC Development Team

---

**Quick Links**:
- 📘 [Main Architecture Doc](REACTIVE_ARCHITECTURE.md)
- 🚀 [Quick Start Guide](QUICK_START_DOCS.md)
- ⚡ [Quick Reference](REACTIVE_PATTERNS_QUICK_REF.md)
- 🗺️ [Architecture Diagram](architecture_diagram.mmd)
- 📊 [Flow Diagrams](data_flow_diagrams.mmd)
