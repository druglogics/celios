"""Helpers to normalize DataFrame sample identifiers to SIDM.

This module is a focused, importable version of the normalization logic
previously embedded inside `features.training.ActivityMatrix`.
"""
from typing import Tuple, List, Dict, Optional

import pandas as pd


def normalize_df_sample_ids(
    df: pd.DataFrame,
    alias_map: Optional[Dict[str, str]] = None,
    sidm_list: Optional[List[str]] = None,
    collapse_method: str = "max",
    verbose: bool = False,
) -> Tuple[pd.DataFrame, List[str]]:
    """Normalize DataFrame sample identifiers (columns) to SIDM.

    Args:
        df: input DataFrame whose columns are sample identifiers
        alias_map: mapping alias -> SIDM (case-sensitive keys expected)
        sidm_list: list of allowed SIDM identifiers
        collapse_method: how to collapse multiple columns mapping to same SIDM ('max'|'mean'|'first')

    Returns:
        (normalized_df, unmapped_list)
    """
    if df is None or df.empty:
        return df, []

    alias_map = alias_map or {}
    mapped_cols = {}
    unmapped = []

    for col in list(df.columns):
        key = str(col).strip()
        # direct SIDM
        if sidm_list and key in sidm_list:
            mapped_cols.setdefault(key, []).append(col)
            continue
        # alias exact match
        if key in alias_map:
            sidm = alias_map[key]
            mapped_cols.setdefault(sidm, []).append(col)
            continue
        # case-insensitive alias match
        low = key.lower()
        found = None
        for a, s in alias_map.items():
            if a.lower() == low:
                found = s
                break
        if found:
            mapped_cols.setdefault(found, []).append(col)
            continue
        # not mapped
        unmapped.append(key)

    if not mapped_cols:
        if verbose:
            print("[MAPPING] No columns mapped to SIDM")
        return df, unmapped

    new_cols = {}
    for sidm, orig_cols in mapped_cols.items():
        if len(orig_cols) == 1:
            new_cols[sidm] = df[orig_cols[0]]
        else:
            if collapse_method == "max":
                new_cols[sidm] = df[orig_cols].max(axis=1, skipna=True)
            elif collapse_method == "mean":
                new_cols[sidm] = df[orig_cols].mean(axis=1, skipna=True)
            else:
                new_cols[sidm] = df[orig_cols].iloc[:, 0]

    new_df = pd.DataFrame(new_cols)

    if verbose:
        duplicates_collapsed = sum(1 for _, cols in mapped_cols.items() if len(cols) > 1)
        print(f"[MAPPING] Mapped SIDM columns: {len(new_df.columns)}")
        print(f"[MAPPING] Unmapped columns: {len(unmapped)}")
        print(f"[MAPPING] SIDMs with duplicate source columns collapsed: {duplicates_collapsed}")
        if unmapped:
            print(f"[MAPPING] Unmapped examples: {unmapped[:10]}")

    return new_df, unmapped
