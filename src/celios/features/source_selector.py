"""Final source selection helpers.

This module applies the source priority rules over a master matrix and can
rename SIDM columns back to user-facing cell-line names when available.
"""
from __future__ import annotations

from typing import Dict, Iterable, List, Mapping, Optional, Sequence

import pandas as pd

from celios.features.master_matrix import SOURCE_PRIORITY, source_column_name


def _rename_sidm_columns_to_cell_line_names(
    df: pd.DataFrame,
    sidm_dict: Optional[Mapping[str, str]] = None,
) -> pd.DataFrame:
    """Rename SIDM columns to cell-line names when a mapping exists."""
    if df is None or not sidm_dict:
        return df

    rename_map = {
        sidm: sidm_dict.get(sidm, sidm)
        for sidm in df.columns
        if sidm in sidm_dict
    }
    if not rename_map:
        return df
    return df.rename(columns=rename_map)


def select_from_master(
    master: pd.DataFrame,
    sidm_list: Sequence[str],
    selected_sources: Optional[Sequence[str]] = None,
    sidm_dict: Optional[Mapping[str, str]] = None,
    include_symbol: bool = True,
    source_priority: Optional[Sequence[str]] = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """Apply source priority to a master matrix and return the final matrix.

    Priority order defaults to mutations > cnv > TF > expression.
    """
    if master is None:
        raise ValueError("master DataFrame must be provided")
    if not sidm_list:
        raise ValueError("sidm_list must be provided")

    selected_sources = list(selected_sources or SOURCE_PRIORITY)
    source_priority = tuple(source_priority or SOURCE_PRIORITY)
    nodes = master.index

    result = pd.DataFrame(index=nodes, columns=list(sidm_list), dtype=float)

    for sidm in sidm_list:
        out = pd.Series(index=nodes, dtype=float)

        for source in source_priority:
            if source not in selected_sources:
                continue
            col = source_column_name(sidm, source)
            if col not in master.columns:
                continue

            values = master[col]
            if source in ("mutations", "cnv"):
                mask = values.isin([0, 1])
                need = out.isna()
                assign = mask & need
                out[assign] = values[assign].astype(float)
            else:
                mask = values.notna()
                need = out.isna()
                assign = mask & need
                out[assign] = values[assign].astype(float)

        result[sidm] = out

    if include_symbol and "symbol" in master.columns:
        result = result.reset_index().rename(columns={"index": "node_name"})
        result.insert(1, "symbol", master["symbol"].reindex(result["node_name"]).values)
        result = result.set_index("node_name")

    result = _rename_sidm_columns_to_cell_line_names(result, sidm_dict=sidm_dict)

    if verbose:
        print(f"[SOURCE] Selected final activity matrix shape: {result.shape}")

    return result


__all__ = [
    "select_from_master",
    "_rename_sidm_columns_to_cell_line_names",
]
