# Welcome to SPAC - Spatial Analysis of Cellular Data

## Overview
SPAC_Shiny is an interactive dashboard designed to empower researchers with intuitive, dynamic, and customizable visualization tools for spatial single-cell datasets derived from cancerous tumors.

## Quick Start Guide

### 1. Load Your Data
- Navigate to the **Data Input** tab
- Upload your spatial single-cell data file
- Supported formats: `.h5ad`, `.pickle`
- Sample data is pre-loaded for demonstration

### 2. Basic Workflow
1. **Data Input**: Load and inspect your dataset
2. **Annotations**: Explore cell type annotations and metadata
3. **Features**: Analyze gene expression patterns
4. **Visualizations**: Create plots and spatial maps
5. **Analysis**: Perform statistical analyses

### 3. Key Features

#### 🔬 **Spatial Analysis**
- Spatial distribution plots
- Neighborhood analysis
- Ripley's L function for spatial statistics

#### 📊 **Statistical Visualizations**
- Boxplots with optimized performance (10x faster)
- Histograms and distribution plots
- Feature vs annotation comparisons

#### 🗺️ **Dimensional Reduction**
- UMAP embeddings
- Interactive scatter plots
- Custom color schemes

#### 🎯 **Data Exploration**
- Interactive filtering and subsetting
- Real-time parameter adjustment
- Export capabilities

## Data Requirements

### Expected Format
- **AnnData objects** (`.h5ad` files)
- **Pickled data** (`.pickle` files)
- Spatial coordinates in `obsm['spatial']`
- Cell type annotations in `obs`

### Anndata Structure
```
adata.obs: Cell metadata and annotations
adata.obsm['spatial']: X,Y coordinates
adata.var: Gene/feature information
adata.X: Expression matrix
```

## Navigation Tips

### Recommended Analysis Order
1. Start with **Data Input** to load your dataset
2. Explore **Annotations** to understand your cell types
3. Check **Features** for gene expression overview
4. Use **Spatial** tab for spatial analysis
5. Create and download plots in visualization tabs

### Performance Notes
- Large datasets (>1M cells) are optimized for speed
- Use filtering options to focus on specific regions
- Boxplot and histogram functions are 10x faster than standard implementations

## Need Help?
- Each tab contains specific instructions
- Hover over input elements for tooltips
- Sample data is provided for testing features
