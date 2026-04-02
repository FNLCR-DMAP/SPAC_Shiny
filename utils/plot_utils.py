"""
Plotting utility functions for SPAC Shiny application.

This module provides reusable helper functions for customizing matplotlib
plots, such as axis label formatting and styling.
"""

import io
from typing import List

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.text import Text


def fig_to_png_bytes(fig) -> bytes:
    """
    Render a matplotlib Figure (or seaborn ClusterGrid) to PNG bytes and
    immediately close the underlying Figure to free memory.

    Parameters
    ----------
    fig : matplotlib.figure.Figure or seaborn.matrix.ClusterGrid
        The figure to serialise.

    Returns
    -------
    bytes
        Raw PNG data.
    """
    buf = io.BytesIO()
    underlying = fig.fig if hasattr(fig, 'fig') else fig
    underlying.savefig(buf, format='png', bbox_inches='tight')
    plt.close(underlying)
    buf.seek(0)
    return buf.read()


def png_bytes_to_figure(png_bytes: bytes):
    """
    Reconstruct a display-ready matplotlib Figure from cached PNG bytes.

    The returned figure wraps the image in a plain imshow axes with
    tight_layout disabled, so it cannot accumulate layout mutations across
    repeated Shiny renders.

    Parameters
    ----------
    png_bytes : bytes
        Raw PNG data produced by :func:`fig_to_png_bytes`.

    Returns
    -------
    matplotlib.figure.Figure
    """
    img = mpimg.imread(io.BytesIO(png_bytes))
    height_px, width_px = img.shape[:2]
    dpi = 100
    fig, ax = plt.subplots(figsize=(width_px / dpi, height_px / dpi), dpi=dpi)
    ax.imshow(img)
    ax.axis('off')
    fig.subplots_adjust(0, 0, 1, 1)
    fig.set_tight_layout(False)
    return fig


def abbreviate_labels(labels: List[Text], limit: int) -> List[str]:
    """
    Abbreviate a list of matplotlib Text objects.

    Parameters
    ----------
    labels : List[matplotlib.text.Text]
        List of label objects to abbreviate.
    limit : int
        Maximum number of characters allowed.

    Returns
    -------
    List[str]
        List of abbreviated label strings.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> ax.set_xticklabels(["LongLabel1", "LongLabel2", "Short"])
    >>> abbreviated = abbreviate_labels(ax.get_xticklabels(), limit=5)
    >>> print(abbreviated)
    ['LongL', 'LongL', 'Short']
    """
    return [label.get_text()[:limit] if label.get_text() else ""
            for label in labels]


def apply_axis_style(labels: List[Text], fontsize: int,
                     fontfamily: str = "DejaVu Sans") -> None:
    """
    Apply font size and family to axis labels.

    Parameters
    ----------
    labels : List[matplotlib.text.Text]
        List of label objects to style.
    fontsize : int
        Font size to apply.
    fontfamily : str, optional
        Font family to apply, by default "DejaVu Sans".

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> ax.set_xlabel("X Label")
    >>> apply_axis_style(ax.get_xticklabels(), fontsize=12)
    """
    for label in labels:
        label.set_fontsize(fontsize)
        label.set_fontfamily(fontfamily)