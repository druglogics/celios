"""Standalone cell-line resolution helper (light shim).

This module provides a small wrapper around the existing resolver in
`celios.utils.cell_line_resolver` so callers in `features` can import a
single entrypoint while we incrementally migrate logic.
"""
from typing import Optional, Tuple, Dict, List

import pandas as pd

from celios.utils import load_csv_file
from celios.utils.cell_line_resolver import (
    resolve_sidm_from_dataframe,
    normalize_identifier,
    detect_identifier_type,
)


def load_cell_line_file(cell_line_file: str, verbose: bool = False):
    """Load a cell-line file with delimiter fallback (csv/tsv).

    Returns a pandas DataFrame.
    """
    if verbose:
        print(f"[SIDM] Loading cell line file: {cell_line_file}")
    try:
        df = load_csv_file(cell_line_file, sep=None, engine="python")
    except Exception:
        df = load_csv_file(cell_line_file)
    if verbose:
        try:
            print(f"[SIDM] Columns: {list(df.columns)}")
            print(f"[SIDM] First rows:\n{df.head(5).to_string(index=False)}")
        except Exception:
            pass
    return df


def normalize_id(value: object) -> str:
    """Wrapper to normalize a single identifier using resolver logic."""
    return normalize_identifier(value)


def detect_id_type(value: object) -> str:
    """Wrapper to detect an identifier type (sidm, model_id, rrid, cvcl, name)."""
    return detect_identifier_type(value)


def resolve_cell_lines(cell_line_file: Optional[str] = None, df: Optional[pd.DataFrame] = None, verbose: bool = False) -> Tuple[List[str], Dict[str, str], Dict[str, str], Dict[str, int]]:
    """Resolve SIDMs from a user-provided cell-line file or DataFrame.

    Returns: (sidm_list, sidm_dict, alias_to_sidm, resolution_report)
    """
    if df is None:
        if not cell_line_file:
            raise ValueError("Either 'cell_line_file' or 'df' must be provided")
        df = load_cell_line_file(cell_line_file, verbose=verbose)

    if verbose:
        print("[SIDM] Beginning resolution of provided cell-line table")
    # Pre-normalize identifier-like columns where helpful (light touch)
    # Keep heavy lifting to resolve_sidm_from_dataframe
    sidm_dict, not_found, resolution_report = resolve_sidm_from_dataframe(df, verbose=verbose)
    if verbose:
        total = resolution_report.get('total_rows') if isinstance(resolution_report, dict) else None
        resolved = resolution_report.get('resolved') if isinstance(resolution_report, dict) else None
        unresolved = resolution_report.get('unresolved') if isinstance(resolution_report, dict) else None
        alias_map = resolution_report.get('alias_to_sidm') if isinstance(resolution_report, dict) else None
        print(f"[SIDM] Resolution summary: total={total}, resolved={resolved}, unresolved={unresolved}")
        print(f"[SIDM] Alias map size: {len(alias_map) if alias_map else 0}")
    alias_map = resolution_report.get("alias_to_sidm", {}) if isinstance(resolution_report, dict) else {}
    sidm_list = list(sidm_dict.keys())
    return sidm_list, sidm_dict, alias_map, resolution_report
