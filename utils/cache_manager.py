"""
Per-session visualization result cache with LRU eviction.

Each Shiny session creates its own VisualizationCache instance (stored in
shared['cache']).  Results are keyed by a composite of:

    (dataset_version, viz_name, normalized_params)

Because dataset_version is part of the key, replacing or subsetting the
AnnData object automatically makes all previous cache entries unreachable
without requiring an explicit flush — though an explicit invalidate() is
also called on dataset change to reclaim memory immediately.

Usage (inside a @render.* function)
-------------------------------------
    cache   = shared['cache']
    version = shared['dataset_version'].get()

    params = {
        'annotation': input.bp_anno(),
        'layer':      input.bp_layer(),
        'features':   tuple(sorted(input.bp_features())),
        # ... all inputs that affect the plot
    }

    def compute():
        fig, df = expensive_visualization(...)
        return fig, df

    fig, df = cache.get_or_compute('boxplot', version, params, compute)
    shared['df_boxplot'].set(df)
    return fig
"""

from __future__ import annotations

import logging
import time
from collections import OrderedDict
from typing import Any, Callable, Dict, Hashable, Tuple

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Parameter normalisation
# ---------------------------------------------------------------------------

def normalize_params(value: Any) -> Hashable:
    """
    Recursively convert an arbitrary value to a hashable representation.

    Rules
    -----
    * dict  → sorted tuple of (key, normalized_value) pairs
    * list  → tuple of normalized values
    * set   → sorted tuple of normalized values
    * tuple → tuple of normalized values
    * everything else is returned as-is (assumed already hashable)
    """
    if isinstance(value, dict):
        return tuple(
            sorted((k, normalize_params(v)) for k, v in value.items())
        )
    if isinstance(value, (list, tuple)):
        return tuple(normalize_params(v) for v in value)
    if isinstance(value, set):
        try:
            return tuple(sorted(normalize_params(v) for v in value))
        except TypeError:
            return tuple(normalize_params(v) for v in value)
    return value


# ---------------------------------------------------------------------------
# Cache class
# ---------------------------------------------------------------------------

class VisualizationCache:
    """
    LRU cache for visualization outputs.

    Parameters
    ----------
    max_size : int
        Maximum number of (fig, df) entries to hold before the oldest
        entry is evicted.  Defaults to 50.

    Notes
    -----
    * Only (fig, df) pairs are cached — AnnData objects are never stored.
    * Thread-safety: Shiny for Python runs in a single-threaded async loop
      per session, so no locking is required here.
    """

    def __init__(self, max_size: int = 50) -> None:
        self._cache: OrderedDict = OrderedDict()
        self.max_size = max_size

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _make_key(
        self,
        viz_name: str,
        dataset_version: int,
        params: Dict,
    ) -> Tuple:
        return (dataset_version, viz_name, normalize_params(params))

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def get_or_compute(
        self,
        viz_name: str,
        dataset_version: int,
        params: Dict,
        compute_fn: Callable[[], Any],
    ) -> Any:
        """
        Return cached result or compute, store, and return a fresh one.

        Parameters
        ----------
        viz_name : str
            Unique name for this visualization (e.g. ``'boxplot'``).
        dataset_version : int
            Current value of ``shared['dataset_version']``.  Changes
            whenever the dataset is replaced or subsetted.
        params : dict
            All inputs that affect the visualization output.  Mutable
            types (lists, dicts, sets) are normalized internally.
        compute_fn : callable
            Zero-argument callable that produces and returns the result
            to cache.  It is called only on a cache miss.

        Returns
        -------
        object
            The cached or freshly computed result.
        """
        t0 = time.perf_counter()
        key = self._make_key(viz_name, dataset_version, params)

        if key in self._cache:
            # Promote to most-recently-used position
            self._cache.move_to_end(key)
            elapsed = time.perf_counter() - t0
            logger.debug(
                "Cache HIT  [%s] (version=%s, cache_size=%d)",
                viz_name, dataset_version, len(self._cache),
            )
            logger.info(
                "Time taken to retrieve %s from cache: %.6f seconds",
                viz_name, elapsed,
            )
            return self._cache[key]

        logger.debug(
            "Cache MISS [%s] (version=%s) — computing ...",
            viz_name, dataset_version,
        )
        result = compute_fn()

        self._cache[key] = result
        self._cache.move_to_end(key)

        # Evict oldest entries if over limit
        while len(self._cache) > self.max_size:
            evicted_key, _ = self._cache.popitem(last=False)
            logger.debug("Cache EVICT [%s]", evicted_key[1])

        return result

    def invalidate(self) -> None:
        """
        Clear the entire cache.

        Call this whenever the underlying dataset is replaced so that
        stale figures and DataFrames are released from memory immediately.
        """
        count = len(self._cache)
        self._cache.clear()
        logger.debug("Cache cleared (%d entries released)", count)

    @property
    def size(self) -> int:
        """Number of entries currently held in the cache."""
        return len(self._cache)

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"VisualizationCache(size={self.size}/{self.max_size})"
        )
