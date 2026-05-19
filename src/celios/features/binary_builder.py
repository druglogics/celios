"""Binary matrix loader and node aggregation helpers.

Responsibilities:
- Load binary mutation/CNV files using `binary_parser` with auto-detection
- Normalize sample IDs to SIDM using `celios.utils.normalize_df_sample_ids`
- Provide node-level aggregation from gene-level binary matrices
"""
from typing import Dict, Optional, Tuple

import pandas as pd

from celios.features.binary_parser import FormatDetector, get_binary_parser
from celios.utils import normalize_df_sample_ids


def load_binary_matrix(
    filepath: str,
    format_override: Optional[str] = None,
    model_registry: Optional[str] = None,
    alias_map: Optional[Dict[str, str]] = None,
    sidm_list: Optional[list] = None,
    collapse_method: str = "max",
    verbose: bool = False,
) -> Tuple[pd.DataFrame, Dict]:
    """Load a binary matrix file and normalize sample IDs to SIDM.

    Returns (genes_df, metadata) where genes_df is genes x SIDM DataFrame.
    """
    # Detect format
    fmt = FormatDetector.detect(filepath, format_override=format_override)
    parser = get_binary_parser(fmt, verbose=verbose)

    genes_df, metadata = parser.load(filepath, model_registry=model_registry)

    # At this point, ModelID-based formats are already mapped to SIDM by the parser.
    # For OLD-style matrices we may need to normalize sample identifiers (aliases -> SIDM).
    df_norm, unmapped = normalize_df_sample_ids(genes_df, alias_map=alias_map, sidm_list=sidm_list, collapse_method=collapse_method, verbose=verbose)

    # If normalization returned a DataFrame with SIDM columns, prefer it
    if df_norm is not None and not df_norm.empty:
        genes_df = df_norm

    # ensure gene index name consistency
    try:
        genes_df.index.name = "gene_symbol"
    except Exception:
        pass

    return genes_df, metadata


def aggregate_genes_to_nodes(genes_df: pd.DataFrame, node_dict: Dict[str, list], agg_method: str = "max") -> pd.DataFrame:
    """Aggregate gene-level binary matrix to node-level.

    Aggregation uses `max` (logical OR) by default so a node is active
    if any mapped gene is active for that SIDM.
    """
    if genes_df is None or genes_df.empty:
        return genes_df
    if not node_dict:
        raise ValueError("node_dict must be provided to aggregate genes to nodes")

    node_rows = {}
    # Ensure genes uppercase for matching
    gene_index = genes_df.index.str.upper()
    df_upper = genes_df.copy()
    df_upper.index = gene_index

    for node, symbols in node_dict.items():
        mapped = [s.upper() for s in symbols if s]
        present = [g for g in mapped if g in df_upper.index]
        if not present:
            # create all-NaN row
            node_rows[node] = pd.Series(index=df_upper.columns, dtype=float)
            continue
        sub = df_upper.loc[present]
        if agg_method == "max":
            node_rows[node] = sub.max(axis=0, skipna=True)
        elif agg_method == "mean":
            node_rows[node] = sub.mean(axis=0, skipna=True)
        else:
            node_rows[node] = sub.iloc[0]

    node_df = pd.DataFrame(node_rows).T
    node_df.index.name = "node"
    return node_df
