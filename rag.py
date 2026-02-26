import ollama
import numpy as np
from scipy.spatial.distance import cosine

PAPER_CHUNKS = [
    {
        "id": "abstract",
        "text": """SPAC (analysis of SPAtial single-Cell datasets) is a scalable web-based 
        platform for single-cell spatial analysis using a four-tier architecture: 
        foundational dependencies, SPAC Python package, interactive analysis layer, 
        and real-time visualization layer (Shiny for Python). Benchmarked on 2.6 million 
        cells from a 4T1 breast tumor model, achieving 20x GPU speedup over CPU."""
    },
    {
        "id": "architecture_browser",
        "text": """SPAC's browser-based access allows researchers to execute analyses 
        directly in a web browser, eliminating local software installations. Bench 
        biologists can launch workflows, tune parameters, and inspect results 
        interactively without command-line expertise."""
    },
    {
        "id": "architecture_hpc",
        "text": """SPAC's scalable HPC/cloud back-end connects to on-premises or 
        cloud-based HPC resources with 150GB+ memory and GPU accelerators. Automated 
        load balancing, parallel computing, and containerized execution allow processing 
        of datasets with millions of cells. Deployed on NIH Biowulf (Slurm). The NIDAP 
        HPC Connector dispatches jobs to CPU/GPU nodes and returns results automatically 
        without requiring HPC credentials."""
    },
    {
        "id": "data_format",
        "text": """SPAC converts CSV input files (from HALO, MCMICRO, QuPath, Visiopharm) 
        into AnnData (H5AD) format. AnnData stores: X (expression matrix, cells x features), 
        obs (per-cell metadata/annotations), var (feature metadata), obsm (embeddings: 
        spatial coordinates, UMAP), layers (raw and normalized), uns (unstructured metadata). 
        Both raw and normalized layers are retained to preserve data lineage."""
    },
    {
        "id": "features_annotations",
        "text": """In SPAC, Features are measurable per-cell attributes such as biomarker 
        intensities (e.g., Hif1a, NOS2, COX2, E-cadherin, vimentin, Ki67, aSMA). 
        Annotations are categorical labels assigned to cells: cell type, spatial region 
        (e.g., hypoxic, normal stroma), or computational clusters (e.g., PhenoGraph). 
        Phenotype codes combine multiple markers: e.g., CD4+CD25+FOXP3+ = regulatory T-cell."""
    },
    {
        "id": "preprocessing",
        "text": """SPAC preprocessing includes arcsinh transformation, quantile scaling, 
        and batch correction. The interface provides per-feature distributions, missingness 
        checks, and outlier detection. Supports cell- and feature-level filtering. 
        Both raw and normalized layers are retained in the AnnData object."""
    },
    {
        "id": "clustering_phenograph",
        "text": """SPAC supports unsupervised clustering via PhenoGraph: 
        (1) normalize marker intensity data, (2) construct k-nearest-neighbor (KNN) graph 
        using Euclidean distance, (3) refine edge weights by shared neighbor count, 
        (4) apply Louvain or Leiden community detection. Key parameters: k (number of 
        neighbors, controls graph connectivity) and resolution (controls cluster granularity, 
        higher = more clusters). Optimal params for 2.6M cell benchmark: k=35, resolution=0.6, 
        yielding 16 biologically meaningful clusters."""
    },
    {
        "id": "gpu_performance",
        "text": """GPU acceleration benchmarks: PhenoGraph clustering on 2.6 million cells 
        with 9 biomarkers took ~2.5 hours on AMD EPYC 7543 CPU vs ~7 minutes on NVIDIA A100 
        GPU (~20x speedup) using the Grapheno implementation. Near-linear scalability 
        demonstrated up to 5 million cells. Running multiple clustering resolutions adds 
        minimal overhead due to GPU parallelism — all k/resolution combinations complete 
        under 30 minutes."""
    },
    {
        "id": "umap_dimensionality",
        "text": """SPAC uses UMAP (Uniform Manifold Approximation and Projection) for 
        dimensionality reduction, projecting high-dimensional cell data into 2D for 
        visualization. Helps identify clusters and subpopulation structure. A flat/blob UMAP 
        may indicate insufficient clustering, parameter tuning needed, or genuinely low 
        subpopulation structure. Also supports t-SNE. Spatial UMAP groups cells by 
        microenvironment similarity."""
    },
    {
        "id": "spatial_nearest_neighbor",
        "text": """Nearest Neighbor analysis in SPAC quantifies the proximity of a specified 
        source phenotype relative to other cell types, highlighting patterns of local 
        adjacency or avoidance. Useful for quantifying immune cell infiltration and 
        cell-cell proximity patterns in tissue."""
    },
    {
        "id": "spatial_ripley",
        "text": """Ripley's L statistic measures whether two phenotypes exhibit spatial 
        aggregation, dispersion, or random distribution across varying spatial scales. 
        Random distribution is modeled using a Poisson point process assuming homogeneous 
        tissue distribution. Positive L value at a given radius = clustering; negative = 
        dispersion. Particularly useful for comparing case vs. control conditions."""
    },
    {
        "id": "spatial_neighborhood",
        "text": """SPAC builds Neighborhood Graphs using Squidpy via KNN or radius-based 
        criteria with customized edge correction. Computes a Cluster Interaction Matrix 
        tallying edges between distinct phenotypes. Applies permutation-based Neighborhood 
        Enrichment scores to determine if phenotype pairs co-locate more or less than 
        expected by chance. Results can be stratified by user-defined annotations 
        (region, timepoint)."""
    },
    {
        "id": "visualization_tools",
        "text": """SPAC visualization tools include: hierarchical heatmaps with dendrograms 
        (z-score normalized marker expressions), Sankey plots (proportional flow between 
        annotations), relational heatmaps (frequency distribution between categorical 
        annotations), spatial plots (static and interactive via Plotly Express with 
        zoom/hover/pinned colors), UMAP plots, boxplots, scatterplots. Interactive spatial 
        plot supports real-time exploration with multiple simultaneous annotation layers."""
    },
    {
        "id": "case_study_results",
        "text": """Case study: 4T1 triple-negative breast tumor mouse model (BALB/c mice, 
        8 biomarkers: Hif1a, NOS2, COX2, β-catenin, vimentin, E-cadherin, Ki67, aSMA + PIMO). 
        Three tumor regions: normoxia (6.8% O2), hypoxia (1.3% O2), necrosis. 
        Key findings: Epithelial cells (E-cadherin+/β-catenin+) localize in normoxic regions 
        (16.1% of normoxic cells). Hypoxic_PIMO-dim cells enriched in hypoxic/necrotic regions 
        (53.2% of hypoxic cells). 312,000 Hypoxic_PIMO-dim cells from hypoxic regions."""
    },
    {
        "id": "knowledge_based_phenotyping",
        "text": """Knowledge-based (manual gating) phenotyping sets intensity thresholds 
        to categorize marker expressions as binary states (0 or 1). Multiple markers 
        combine into composite phenotype_codes (e.g., CD4+CD25+FOXP3+ for regulatory T-cells). 
        Domain experts specify phenotype_name and phenotype_code without custom programming. 
        Ensures consistency across multi-slide datasets."""
    },
    {
        "id": "availability",
        "text": """SPAC availability: 
        GitHub (Python package): https://github.com/FNLCR-DMAP/SCSAWorkflow
        GitHub (Shiny dashboard): https://github.com/FNLCR-DMAP/SPAC_Shiny
        Hosted app: https://appshare.cancer.gov/spac-interactive-visualization/
        Galaxy tools: suite_spac_tools on Galaxy Main ToolShed
        License: MIT. Python requirements: 3.9.13, NumPy 1.26.4, pandas 1.5, 
        AnnData 0.10, Scanpy 1.9, Matplotlib 3.9.2."""
    },
]
client = ollama.Client(host="http://host.docker.internal:11434")

def embed(text: str) -> np.ndarray:
    response = client.embeddings(model="nomic-embed-text", prompt=text)
    return np.array(response["embedding"])

texts = [chunk["text"] for chunk in PAPER_CHUNKS]
ids = [chunk["id"] for chunk in PAPER_CHUNKS]
embeddings = [embed(text) for text in texts]

def retrieve_context(query: str, n_results: int = 3) -> str:
    query_embedding = embed(query)
    scores = [1 - cosine(query_embedding, emb) for emb in embeddings]
    top_indices = np.argsort(scores)[::-1][:n_results]
    chunks = [texts[i] for i in top_indices]
    return "\n\n---\n\n".join(chunks)