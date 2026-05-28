"""Expression matrix loader and normalization helpers.

Responsibilities:
- Reuse `activity_parser` format detection and parsing
- Normalize sample IDs to SIDM before filtering
- Preserve row-wise log + min-max normalization
"""
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

from celios.features.activity_parser import FormatDetector, get_parser
from celios.utils import normalize_df_sample_ids


def load_expression_matrix(
    activity_file: str,
    format_override: Optional[str] = None,
    verbose: bool = False,
    deep_debug: bool = False,
) -> Tuple[pd.DataFrame, Dict]:
    """Load and parse an expression matrix using the existing parser stack."""
    if verbose and deep_debug:
        print(f"[EXPRESSION] Loading activity file: {activity_file}")

    detected_format = FormatDetector.detect(activity_file, format_override=format_override)
    if verbose and deep_debug:
        print(f"[EXPRESSION] Detected format: {detected_format}")

    parser = get_parser(detected_format, verbose=verbose)
    df, metadata = parser.load(activity_file)
    metadata = metadata or {}
    metadata["format"] = metadata.get("format", detected_format)

    if verbose and deep_debug:
        print(f"[EXPRESSION] Parsed raw matrix shape: {df.shape}")

    return df, metadata


def prepare_expression_matrix(
    activity_file: Optional[str] = None,
    activity_df: Optional[pd.DataFrame] = None,
    sidm_list: Optional[list] = None,
    alias_map: Optional[Dict[str, str]] = None,
    format_override: Optional[str] = None,
    verbose: bool = False,
    deep_debug: bool = False,
) -> Tuple[pd.DataFrame, Dict]:
    """Load expression data, normalize sample IDs, and retain SIDM columns only."""
    if activity_df is None:
        if not activity_file:
            raise ValueError("Either activity_file or activity_df must be provided")
        df, metadata = load_expression_matrix(
            activity_file,
            format_override=format_override,
            verbose=verbose,
            deep_debug=deep_debug,
        )
    else:
        df = activity_df.copy()
        metadata = {"format": format_override or "unknown"}

    if df.index.name != "symbol":
        df.index.name = "symbol"
    df.index = df.index.astype(str).str.upper()

    df_norm, unmapped = normalize_df_sample_ids(
        df,
        alias_map=alias_map,
        sidm_list=sidm_list,
        collapse_method="mean",
        verbose=verbose,
    )

    if unmapped and verbose and deep_debug:
        print(f"[EXPRESSION] Unmapped sample IDs: {unmapped[:20]}")

    if df_norm is not None and not df_norm.empty:
        df = df_norm

    if sidm_list is not None:
        cols = [c for c in sidm_list if c in df.columns]
        df = df.loc[:, cols]

    if verbose and deep_debug:
        print(f"[EXPRESSION] Prepared matrix shape: {df.shape}")

    return df, metadata


def normalize_expression_matrix(df: pd.DataFrame, verbose: bool = False, deep_debug: bool = False) -> pd.DataFrame:
    """Apply log transform and row-wise min-max normalization to expression data."""
    if df is None:
        return None

    out = df.copy(deep=True)
    offset = 0.01
    numeric_cols = out.select_dtypes(include=[float, int]).columns
    out.loc[:, numeric_cols] = out.loc[:, numeric_cols].apply(lambda x: np.log(x + offset))
    out.loc[:, numeric_cols] = out.loc[:, numeric_cols].apply(
        lambda x: (x - x.min()) / (x.max() - x.min()),
        axis=1,
    )

    if verbose and deep_debug:
        print(f"[EXPRESSION] Normalized matrix shape: {out.shape}")

    return out
