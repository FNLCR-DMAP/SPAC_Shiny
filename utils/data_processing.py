"""
SPAC Shiny App - Data Processing Module
This module contains functions for loading and processing data for the SPAC Shiny app.
"""

from pathlib import Path
import logging
import pickle

import anndata as ad
import pandas as pd

logger = logging.getLogger(__name__)

# Simple in-memory cache for loaded datasets
_data_cache = {}


def cached_load_data(file_path):
    """
    Load AnnData with caching to avoid repeated file I/O.

    This function caches loaded datasets in system RAM, providing
    significant performance improvements when switching between
    analysis modules or reloading the same dataset.

    Parameters
    ----------
    file_path : str or Path
        Path to the data file (.h5ad, .pickle, or other AnnData format)

    Returns
    -------
    anndata.AnnData
        Loaded AnnData object (cached if previously loaded)

    Notes
    -----
    Cache is stored in system RAM and persists until app restart.
    Subsequent loads of the same file are nearly instant.

    Examples
    --------
    >>> adata = cached_load_data("data/sample.h5ad")
    Loading data from data/sample.h5ad
    >>> # Second call with same file - instant from cache
    >>> adata = cached_load_data("data/sample.h5ad")
    Retrieved data from cache: data/sample.h5ad
    """
    # Convert to absolute path for consistent cache keys
    cache_key = str(Path(file_path).resolve())

    # Check cache first
    if cache_key in _data_cache:
        logger.info(f"Retrieved data from cache: {file_path}")
        return _data_cache[cache_key]

    # Load data from file
    logger.info(f"Loading data from {file_path}")

    with open(file_path, 'rb') as file:
        if file_path.endswith('.pickle'):
            adata = pickle.load(file)
        elif file_path.endswith('.h5ad'):
            adata = ad.read_h5ad(file_path)
        else:
            adata = ad.read(file_path)

    # Store in cache
    _data_cache[cache_key] = adata
    logger.info(
        f"Cached data: {adata.n_obs} cells, "
        f"{adata.n_vars} genes"
    )

    return adata


def clear_data_cache():
    """
    Clear the data cache to free system RAM.

    This removes all cached datasets from memory. Useful if you need
    to free up RAM or reload data from disk.

    Examples
    --------
    >>> clear_data_cache()
    Cleared data cache
    """
    global _data_cache
    _data_cache.clear()
    logger.info("Cleared data cache")


def get_cache_info():
    """
    Get information about currently cached datasets.

    Returns
    -------
    dict
        Dictionary with cache statistics including number of cached
        files and total memory usage

    Examples
    --------
    >>> info = get_cache_info()
    >>> print(f"Cached files: {info['num_files']}")
    Cached files: 2
    """
    import sys

    total_size = sum(
        sys.getsizeof(adata)
        for adata in _data_cache.values()
    )

    return {
        'num_files': len(_data_cache),
        'cached_files': list(_data_cache.keys()),
        'total_size_mb': total_size / (1024 * 1024)
    }


def read_html_file(filepath):
    """Read HTML file content"""
    try:
        with open(filepath, 'r', encoding='utf-8') as file:
            return file.read()
    except FileNotFoundError:
        return ""


def read_markdown_file(filepath):
    """Read Markdown file content and convert to HTML"""
    try:
        with open(filepath, 'r', encoding='utf-8') as file:
            content = file.read()
            
            def process_inline_formatting(text):
                """Process bold, italic, and code formatting"""
                # Handle bold text **text**
                import re
                text = re.sub(r'\*\*(.*?)\*\*', r'<strong>\1</strong>', text)
                # Handle italic text *text*
                text = re.sub(r'\*(.*?)\*', r'<em>\1</em>', text)
                # Handle inline code `code`
                text = re.sub(r'`(.*?)`', r'<code>\1</code>', text)
                return text
            
            lines = content.split('\n')
            html_lines = []
            in_ordered_list = False
            in_unordered_list = False
            
            for line in lines:
                original_line = line
                line = line.strip()
                
                # Handle headers
                if line.startswith('#### '):
                    if in_ordered_list:
                        html_lines.append('</ol>')
                        in_ordered_list = False
                    if in_unordered_list:
                        html_lines.append('</ul>')
                        in_unordered_list = False
                    header_text = process_inline_formatting(line[5:])
                    html_lines.append(f'<h4>{header_text}</h4>')
                elif line.startswith('### '):
                    if in_ordered_list:
                        html_lines.append('</ol>')
                        in_ordered_list = False
                    if in_unordered_list:
                        html_lines.append('</ul>')
                        in_unordered_list = False
                    header_text = process_inline_formatting(line[4:])
                    html_lines.append(f'<h3>{header_text}</h3>')
                elif line.startswith('## '):
                    if in_ordered_list:
                        html_lines.append('</ol>')
                        in_ordered_list = False
                    if in_unordered_list:
                        html_lines.append('</ul>')
                        in_unordered_list = False
                    header_text = process_inline_formatting(line[3:])
                    html_lines.append(f'<h2>{header_text}</h2>')
                elif line.startswith('# '):
                    if in_ordered_list:
                        html_lines.append('</ol>')
                        in_ordered_list = False
                    if in_unordered_list:
                        html_lines.append('</ul>')
                        in_unordered_list = False
                    header_text = process_inline_formatting(line[2:])
                    html_lines.append(f'<h1>{header_text}</h1>')
                
                # Handle numbered lists (1. 2. 3. etc.)
                elif line and line[0].isdigit() and '. ' in line:
                    if in_unordered_list:
                        html_lines.append('</ul>')
                        in_unordered_list = False
                    if not in_ordered_list:
                        html_lines.append('<ol>')
                        in_ordered_list = True
                    list_text = line.split('. ', 1)[1]
                    formatted_text = process_inline_formatting(list_text)
                    html_lines.append(f'<li>{formatted_text}</li>')
                
                # Handle bullet lists (- or *)
                elif line.startswith('- ') or line.startswith('* '):
                    if in_ordered_list:
                        html_lines.append('</ol>')
                        in_ordered_list = False
                    if not in_unordered_list:
                        html_lines.append('<ul>')
                        in_unordered_list = True
                    list_text = line[2:]
                    formatted_text = process_inline_formatting(list_text)
                    html_lines.append(f'<li>{formatted_text}</li>')
                
                # Handle empty lines
                elif not line:
                    if in_ordered_list:
                        html_lines.append('</ol>')
                        in_ordered_list = False
                    if in_unordered_list:
                        html_lines.append('</ul>')
                        in_unordered_list = False
                    html_lines.append('<br>')
                
                # Handle regular paragraphs
                else:
                    if in_ordered_list:
                        html_lines.append('</ol>')
                        in_ordered_list = False
                    if in_unordered_list:
                        html_lines.append('</ul>')
                        in_unordered_list = False
                    formatted_text = process_inline_formatting(line)
                    html_lines.append(f'<p>{formatted_text}</p>')
            
            # Close any remaining lists
            if in_ordered_list:
                html_lines.append('</ol>')
            if in_unordered_list:
                html_lines.append('</ul>')
            
            return '\n'.join(html_lines)
            
    except FileNotFoundError:
        return "<div class='alert alert-warning'>Content not found.</div>"
    except Exception as e:
        return f"<div class='alert alert-danger'>Error: {str(e)}</div>"

def load_data(file_path):
    """
    Load data from a specified file path. Supports .h5ad and .pickle formats.
    
    Parameters:
        file_path (str): The path to the data file.

    Returns:
        adata: Loaded AnnData object or None if loading fails.
    """
    try:
        with open(file_path, 'rb') as file:
            if file_path.endswith('.h5ad'):
                try:
                    adata = ad.read_h5ad(file)
                    return adata
                except Exception as e:
                    print(f"Error loading .h5ad file: {e}")
                    return None
            elif file_path.endswith('.pickle'):
                try:
                    adata = pd.read_pickle(file)
                    return adata
                except Exception as e:
                    print(f"Error loading .pickle file: {e}")
                    return None
            else:
                print("Unsupported file format. Please provide a .h5ad or .pickle file.")
                return None
    except FileNotFoundError:
        print(
            f"Preloaded data file {file_path} not found.",
            "Proceeding without preloaded data."
        )


def get_annotation_label_counts(adata: ad.AnnData):
    """
    Return a dictionary of every annotation (column in adata.obs),
    where the value is a dict of {label: cell_count}.

    Parameters:
        adata (AnnData): AnnData object containing the data.

    Returns:
        dict: A dictionary where keys are annotation labels and values 
        are dictionaries of {label: cell_count} for each annotation.

    Example structure:
      {
        "cell_type": {"T-cell": 100, "B-cell": 80, ...},
        "condition": {"disease": 120, "healthy": 60, ...},
        ...
      }
    """
    if adata is None or not hasattr(adata, "obs") or adata.obs.empty:
        return {}

    annotation_counts = {}
    for col in adata.obs.columns:
        # value_counts returns a Series of {label: count}
        vc = adata.obs[col].value_counts(dropna=False)
        annotation_counts[col] = vc.to_dict()

    return annotation_counts


def get_annotation_top_labels(
    adata: ad.AnnData, annotation: str, top_n: int = 10
):
    """
    Return the top labels and their counts for a specific annotation column.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix.
    annotation : str
        Column name in ``adata.obs`` to extract labels from.
    top_n : int or None, optional
        Number of top labels to return. If None, return all labels.
        Default is 10.

    Returns
    -------
    list of (label, count)
        A list of tuples sorted by count descending. Returns an empty list
        if the annotation is not present or no data available.
    """
    if adata is None or not hasattr(adata, "obs"):
        return []

    if annotation not in adata.obs.columns:
        return []

    vc = adata.obs[annotation].value_counts(dropna=False)
    items = list(vc.items()) if hasattr(vc, 'items') else list(vc.iteritems())
    # items is list of (label, count)
    # sort descending by count (value_counts is already sorted but ensure it)
    items_sorted = sorted(items, key=lambda x: x[1], reverse=True)

    if top_n is None:
        return items_sorted
    return items_sorted[:top_n]


def get_rl_pairs(adata: ad.AnnData):
    """
    Extract Ripley phenotype pairs from ``adata.uns['ripley_l']``.

    The function expects ``adata.uns['ripley_l']`` to be a DataFrame-like
    structure with columns ``center_phenotype`` and ``neighbor_phenotype``.
    It returns a list of strings formatted as "CENTER -> NEIGHBOR".

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix.

    Returns
    -------
    list[str]
        List of unique pair strings. Empty list if nothing found.
    """
    if adata is None or not hasattr(adata, "uns"):
        return []

    ripley_results = None
    try:
        ripley_results = adata.uns.get('ripley_l')
    except Exception:
        ripley_results = None

    if ripley_results is None:
        return []

    try:
        unique_df = (
            ripley_results[["center_phenotype", "neighbor_phenotype"]]
            .drop_duplicates()
        )
        choices = [
            f"{str(row[0])} -> {str(row[1])}"
            for _, row in unique_df.iterrows()
        ]
    except Exception:
        choices = []

    return choices


def get_spatial_distance_columns(adata: ad.AnnData):
    """
    Return the column names for a spatial distance matrix stored in the
    AnnData object. Checks both ``adata.obsm['spatial_distance']`` and
    ``adata.uns['spatial_distance']``.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix.

    Returns
    -------
    list[str] or None
        List of column names if available, otherwise None.
    """
    if adata is None:
        return None

    # Prefer obsm
    try:
        if hasattr(adata, 'obsm') and 'spatial_distance' in adata.obsm:
            distance_df = adata.obsm['spatial_distance']
            if hasattr(distance_df, 'columns'):
                return list(distance_df.columns)
    except Exception:
        pass

    try:
        if hasattr(adata, 'uns') and 'spatial_distance' in adata.uns:
            distance_df = adata.uns['spatial_distance']
            if hasattr(distance_df, 'columns'):
                return list(distance_df.columns)
    except Exception:
        pass

    return None
