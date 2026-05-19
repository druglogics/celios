"""TF activity matrix helper.

Responsibilities:
- Load TF activity data
- Preserve p-value filtering
- Map conditions to SIDM using resolved cell-line metadata
- Map TF sources to nodes using node dictionaries
- Preserve binary threshold conversion
"""
from typing import Dict, List, Optional

import numpy as np
import pandas as pd


def _map_condition_to_sidm(value: object, sidm_dict: Optional[Dict[str, str]] = None, alias_map: Optional[Dict[str, str]] = None):
    """Map a TF condition value to SIDM using display-name and alias maps."""
    if value is None:
        return None

    text = str(value)
    if sidm_dict:
        reverse = {display_name: sidm for sidm, display_name in sidm_dict.items()}
        if text in reverse:
            return reverse[text]

    alias_map = alias_map or {}
    if text in alias_map:
        return alias_map[text]

    lowered = text.lower()
    for alias, sidm in alias_map.items():
        if str(alias).lower() == lowered:
            return sidm

    return None


def load_tf_matrix(
    tf_activity_file: str,
    node_dict_reversed: Dict[str, List[str]],
    sidm_list: List[str],
    sidm_dict: Optional[Dict[str, str]] = None,
    alias_map: Optional[Dict[str, str]] = None,
    p_value_threshold: float = 0.05,
    binary_threshold: float = 0.0,
    verbose: bool = False,
) -> pd.DataFrame:
    """Load TF activity and aggregate it to a node x SIDM matrix.

    Returns a DataFrame indexed by node_name, columns = SIDM IDs.
    """
    if not tf_activity_file:
        return None
    if not node_dict_reversed:
        raise ValueError("node_dict_reversed must be provided to map TF sources to nodes")
    if not sidm_list:
        raise ValueError("sidm_list must be provided to build the TF matrix")

    if verbose:
        print(f"[TF] Loading TF activity file: {tf_activity_file}")

    tf = pd.read_csv(tf_activity_file, sep=r"\s+")

    if "p_value" in tf.columns:
        tf = tf[tf["p_value"] < p_value_threshold]
        if verbose:
            print(f"[TF] Rows after p-value filtering (< {p_value_threshold}): {len(tf)}")

    if "condition" not in tf.columns or "source" not in tf.columns or "score" not in tf.columns:
        raise ValueError("TF activity file must contain at least 'condition', 'source', and 'score' columns")

    tf["sidm"] = tf["condition"].map(lambda x: _map_condition_to_sidm(x, sidm_dict=sidm_dict, alias_map=alias_map))
    unmapped_conditions = tf[tf["sidm"].isna()]["condition"].unique().tolist()
    if unmapped_conditions and verbose:
        print(f"[TF] Unmapped conditions: {unmapped_conditions[:20]}")

    tf = tf.dropna(subset=["sidm"])

    # map TF source -> node_name using node_dict_reversed
    tf["node_name"] = tf["source"].map(node_dict_reversed)
    tf = tf.explode("node_name").dropna(subset=["node_name"])

    pivot = tf.pivot_table(index="node_name", columns="sidm", values="score", aggfunc="max")

    for sidm in sidm_list:
        if sidm not in pivot.columns:
            pivot[sidm] = np.nan
    pivot = pivot.loc[:, sidm_list]

    thresh = binary_threshold
    binary = (pivot > thresh).astype(float)
    node_activity = pivot.where(pivot.isna(), binary)

    if verbose:
        print(f"[TF] Node activity matrix shape: {node_activity.shape}")

    return node_activity
