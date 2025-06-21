"""
SPAC Shiny App - Data Processing Module
This module contains functions for loading and processing data for the SPAC Shiny app.
"""

import anndata as ad
import pandas as pd

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
