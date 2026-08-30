"""Utilities for consistent download filename generation."""

from __future__ import annotations

from typing import Mapping


MIME_TO_EXTENSION: Mapping[str, str] = {
    "text/csv": "csv",
    "text/html": "html",
    "text/plain": "txt",
    "text/tab-separated-values": "tsv",
    "application/json": "json",
    "application/pdf": "pdf",
    "application/xml": "xml",
    "application/zip": "zip",
    "image/png": "png",
    "image/jpeg": "jpg",
    "image/webp": "webp",
}


def infer_extension(mime_type: str | None, default_extension: str = "bin") -> str:
    """Infer file extension from MIME type."""
    if not mime_type:
        return default_extension
    return MIME_TO_EXTENSION.get(mime_type.strip().lower(), default_extension)


def build_download_filename(
    shared,
    graph_type: str,
    extension: str | None = None,
    mime_type: str | None = None,
    default_extension: str = "bin",
) -> str:
    """Build filenames in the format dataset_graphtype.ext.

    For future file types, callers can pass only `mime_type` and avoid hardcoding
    extensions in each server module.
    """
    dataset_name = shared.get("input_filename")
    if hasattr(dataset_name, "get"):
        dataset_name = dataset_name.get()

    if not dataset_name:
        dataset_name = "dataset"

    clean_graph_type = str(graph_type).strip().replace(" ", "_").lower()
    resolved_extension = extension or infer_extension(mime_type, default_extension)
    clean_extension = str(resolved_extension).strip().lstrip(".").lower()
    return f"{dataset_name}_{clean_graph_type}.{clean_extension}"
