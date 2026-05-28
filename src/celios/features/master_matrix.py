"""Master matrix assembly helpers.

This module assembles the wide `<SIDM>__<source>` master matrix from
already-prepared source matrices. It intentionally does not load files or
resolve identifiers; callers must provide precomputed inputs.
"""
from __future__ import annotations

from typing import Dict, Iterable, List, Mapping, Optional, Sequence

import numpy as np
import pandas as pd


SOURCE_PRIORITY = ("mutations", "cnv", "TF", "expression")


def source_column_name(sidm: str, source: str) -> str:
    """Return the canonical master-matrix column name for a source block."""
    return f"{sidm}__{source}"


def build_symbol_column(node_dict: Mapping[str, Sequence[object]], node_index: Iterable[str]) -> pd.Series:
    """Build the node-level `symbol` column from node_dict values."""
    node_symbols = {str(node): syms for node, syms in (node_dict or {}).items()}
    index = list(node_index)
    return pd.Series([node_symbols.get(node) for node in index], index=index, name="symbol")


def _source_block_from_node_matrix(
    node_matrix: Optional[pd.DataFrame],
    sidm_list: Sequence[str],
    node_index: Sequence[str],
    source: str,
    fill_value=np.nan,
) -> pd.DataFrame:
    """Convert a node x SIDM matrix into a `<SIDM>__<source>` block."""
    if node_matrix is None:
        block = pd.DataFrame(index=node_index)
        for sidm in sidm_list:
            block[source_column_name(sidm, source)] = fill_value
        return block

    aligned = node_matrix.reindex(node_index)
    data = {}
    for sidm in sidm_list:
        if sidm in aligned.columns:
            data[source_column_name(sidm, source)] = aligned[sidm].values
        else:
            data[source_column_name(sidm, source)] = fill_value
    return pd.DataFrame(data, index=node_index)


def assemble_master_matrix(
    node_dict: Mapping[str, Sequence[object]],
    sidm_list: Sequence[str],
    source_matrices: Optional[Mapping[str, Optional[pd.DataFrame]]] = None,
    node_index: Optional[Sequence[str]] = None,
    include_symbol: bool = True,
    source_order: Optional[Sequence[str]] = None,
    selected_sources: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """Assemble the master matrix from precomputed node-level source matrices.

    Args:
        node_dict: mapping of node -> symbols
        sidm_list: ordered list of SIDM identifiers
        source_matrices: mapping of source name -> node x SIDM DataFrame
        node_index: optional node ordering. Defaults to node_dict keys.
        include_symbol: whether to include a `symbol` column
        source_order: explicit order of sources; defaults to SOURCE_PRIORITY
        selected_sources: optional list of sources to include; if provided, only these sources are assembled
    """
    if not node_dict:
        raise ValueError("node_dict must be provided to assemble the master matrix")
    if not sidm_list:
        raise ValueError("sidm_list must be provided to assemble the master matrix")

    source_matrices = source_matrices or {}
    source_order = tuple(source_order or SOURCE_PRIORITY)
    node_index = list(node_index or node_dict.keys())
    
    # If selected_sources is provided, only include those sources
    if selected_sources is not None:
        selected_sources = set(selected_sources)
        source_order = tuple(s for s in source_order if s in selected_sources)

    master = pd.DataFrame(index=node_index)
    if include_symbol:
        master["symbol"] = build_symbol_column(node_dict, node_index)

    for source in source_order:
        block = _source_block_from_node_matrix(
            source_matrices.get(source),
            sidm_list=sidm_list,
            node_index=node_index,
            source=source,
        )
        master = pd.concat([master, block], axis=1)

    return master


def add_source_block(
    master: pd.DataFrame,
    node_matrix: Optional[pd.DataFrame],
    sidm_list: Sequence[str],
    source: str,
    node_index: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """Append a source block to an existing master matrix.

    This keeps repeated source-column code out of callers that still build
    the master matrix incrementally.
    """
    if master is None:
        raise ValueError("master must be provided")
    node_index = list(node_index or master.index)
    block = _source_block_from_node_matrix(node_matrix, sidm_list, node_index, source)
    return pd.concat([master, block], axis=1)


__all__ = [
    "SOURCE_PRIORITY",
    "add_source_block",
    "assemble_master_matrix",
    "build_symbol_column",
    "source_column_name",
]
