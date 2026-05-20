"""Binary matrix loader and node aggregation helpers.

Responsibilities:
- Load binary mutation/CNV files using `binary_parser` with auto-detection
- Normalize sample IDs to SIDM using the pre-resolved alias map
- Provide node-level aggregation from gene-level binary matrices
"""
from typing import Dict, List, Optional, Tuple

import pandas as pd

from celios.features.binary_parser import FormatDetector, get_binary_parser


def _normalize_binary_samples_locally(
    genes_df: pd.DataFrame,
    alias_map: Optional[Dict[str, str]] = None,
    sidm_list: Optional[List[str]] = None,
    collapse_method: str = "max",
    sample_axis: str = "columns",
    verbose: bool = False,
) -> Tuple[pd.DataFrame, List[str], List[str]]:
    """Keep only user-relevant samples and rename them to SIDM locally."""
    if genes_df is None or genes_df.empty:
        return genes_df, [], []

    if sample_axis not in {"columns", "rows"}:
        raise ValueError("sample_axis must be either 'columns' or 'rows'")

    alias_map = alias_map or {}
    sidm_set = set(sidm_list or [])

    mapped_sources: Dict[str, List[str]] = {}
    unmapped = []

    sample_labels = genes_df.columns if sample_axis == "columns" else genes_df.index

    for sample_label in sample_labels:
        key = str(sample_label).strip()
        sidm = None
        if key in sidm_set:
            sidm = key
        elif key in alias_map:
            sidm = alias_map[key]
        else:
            lowered = key.lower()
            for alias, mapped in alias_map.items():
                if alias.lower() == lowered:
                    sidm = mapped
                    break

        if sidm and (not sidm_set or sidm in sidm_set):
            mapped_sources.setdefault(sidm, []).append(sample_label)
        else:
            unmapped.append(key)

    if not mapped_sources:
        if sample_axis == "columns":
            return genes_df.iloc[:, 0:0].copy(), unmapped, []
        return genes_df.iloc[0:0, :].copy(), unmapped, []

    new_columns = {}
    ordered_sidm = sidm_list or list(mapped_sources.keys())
    for sidm in ordered_sidm:
        source_labels = mapped_sources.get(sidm)
        if not source_labels:
            continue
        if sample_axis == "columns":
            if len(source_labels) == 1:
                new_columns[sidm] = genes_df.loc[:, source_labels[0]]
            elif collapse_method == "mean":
                new_columns[sidm] = genes_df.loc[:, source_labels].mean(axis=1, skipna=True)
            elif collapse_method == "first":
                new_columns[sidm] = genes_df.loc[:, source_labels[0]]
            else:
                new_columns[sidm] = genes_df.loc[:, source_labels].max(axis=1, skipna=True)
        else:
            if len(source_labels) == 1:
                new_columns[sidm] = genes_df.loc[source_labels[0]]
            elif collapse_method == "mean":
                new_columns[sidm] = genes_df.loc[source_labels].mean(axis=0, skipna=True)
            elif collapse_method == "first":
                new_columns[sidm] = genes_df.loc[source_labels[0]]
            else:
                new_columns[sidm] = genes_df.loc[source_labels].max(axis=0, skipna=True)

    if sample_axis == "columns":
        selected = pd.DataFrame(new_columns, index=genes_df.index)
        matched = list(selected.columns)
    else:
        selected = pd.DataFrame(new_columns).T
        matched = list(selected.index)

    if verbose:
        print(f"[MAPPING] Selected SIDM columns: {len(matched)}")
        if sample_axis == "rows":
            matched_model_ids = [str(label) for labels in mapped_sources.values() for label in labels]
            print(f"[MAPPING] first 10 mutation ModelIDs matched: {matched_model_ids[:10]}")
        if unmapped:
            print(f"[MAPPING] Unmapped binary sample IDs: {unmapped[:10]}")

    return selected, unmapped, matched


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

    genes_df, metadata = parser.load(
        filepath,
        model_registry=model_registry,
        resolve_model_ids=False,
    )

    sample_axis = (metadata or {}).get("sample_axis", "columns")
    raw_sample_ids = (metadata or {}).get("raw_sample_ids")
    input_shape = (metadata or {}).get("input_shape", genes_df.shape)

    if verbose:
        alias_map = alias_map or {}
        ach_aliases = sorted(alias for alias in alias_map if str(alias).upper().startswith("ACH-"))
        mutation_sample_ids = [str(sample_id).strip() for sample_id in (raw_sample_ids or [])]
        alias_key_set = {str(alias).strip() for alias in alias_map.keys()}
        intersection = sorted(set(mutation_sample_ids) & alias_key_set)
        print(f"[BINARY] alias_map total size: {len(alias_map)}")
        print(f"[BINARY] ACH alias count: {len(ach_aliases)}")
        print(f"[BINARY] first 20 ACH aliases: {ach_aliases[:20]}")
        print(f"[BINARY] first 20 sample IDs from mutation file: {mutation_sample_ids[:20]}")
        print(f"[BINARY] mutation sample ID / alias_map key intersection: {intersection[:20]}")

    if verbose:
        print(f"[BINARY] detected binary format: {fmt}")
        print(f"[BINARY] input shape: {input_shape}")
        print(f"[BINARY] detected sample axis: {sample_axis}")
        if raw_sample_ids is None:
            raw_sample_ids = list(genes_df.columns if sample_axis == "columns" else genes_df.index)
        print(f"[BINARY] first 10 sample IDs: {list(raw_sample_ids)[:10]}")

    genes_df, unmapped, matched_sample_ids = _normalize_binary_samples_locally(
        genes_df,
        alias_map=alias_map,
        sidm_list=sidm_list,
        collapse_method=collapse_method,
        sample_axis=sample_axis,
        verbose=verbose,
    )

    if sample_axis == "rows":
        genes_df = genes_df.T

    try:
        genes_df.index.name = "gene_symbol"
    except Exception:
        pass

    if isinstance(metadata, dict):
        metadata = dict(metadata)
        metadata["unmapped_sample_ids"] = unmapped
        metadata["selected_sample_ids"] = matched_sample_ids
        metadata["output_shape"] = tuple(genes_df.shape)
        metadata["sample_axis"] = sample_axis

    if verbose:
        print(f"[BINARY] matched sample IDs: {matched_sample_ids[:10]}")
        print(f"[BINARY] output shape after SIDM filtering: {genes_df.shape}")

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
