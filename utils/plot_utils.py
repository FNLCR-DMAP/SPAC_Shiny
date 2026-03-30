"""
Plotting utility functions for SPAC Shiny application.

This module provides reusable helper functions for customizing matplotlib
plots, such as axis label formatting and styling.
"""

from typing import List
from matplotlib.text import Text


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