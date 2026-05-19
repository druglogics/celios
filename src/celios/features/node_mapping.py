"""Stateless helpers for node dictionary loading and symbol mapping.

This module extracts the reusable pieces around node dictionaries:
- loading a node_dict from a dict, DataFrame, or file path
- reversing node_dict to a gene/symbol -> node mapping
- normalizing gene symbols for consistent matching
- looking up nodes for a gene symbol or list of symbols
"""
from __future__ import annotations

import ast
import re
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

import pandas as pd

from celios.utils.io import load_node_dict_from_csv


def normalize_symbol(value: object) -> str:
    """Normalize a gene or symbol string for matching.

    Rules are intentionally light-touch:
    - convert to string
    - NFKC-normalize through Python string handling where possible
    - trim whitespace
    - collapse internal whitespace to a single space
    - uppercase
    """
    if value is None:
        return ""
    text = str(value).strip()
    if not text:
        return ""
    text = re.sub(r"\s+", " ", text)
    return text.upper()


def normalize_symbol_list(values: Iterable[object]) -> List[str]:
    """Normalize a sequence of symbols and drop empty entries."""
    return [symbol for symbol in (normalize_symbol(value) for value in values) if symbol]


def load_node_dict(node_dict_input: Optional[object], verbose: bool = False) -> Dict[str, List[str]]:
    """Load a node dictionary from a path, dict, DataFrame, or iterable-like input.

    Accepted forms:
    - dict: returned after normalizing values to lists of strings
    - pandas.DataFrame: first two columns are treated as node and symbol columns
    - str/path-like: loaded via `load_node_dict_from_csv`

    Returns a mapping of node_name -> list[str].
    """
    if node_dict_input is None:
        raise ValueError("node_dict must be provided to map genes to nodes")

    if isinstance(node_dict_input, dict):
        cleaned: Dict[str, List[str]] = {}
        for node, symbols in node_dict_input.items():
            if symbols is None:
                cleaned[str(node)] = []
            elif isinstance(symbols, (list, tuple, set)):
                cleaned[str(node)] = normalize_symbol_list(symbols)
            else:
                cleaned[str(node)] = normalize_symbol_list([symbols])
        return cleaned

    if isinstance(node_dict_input, pd.DataFrame):
        df = node_dict_input.copy()
        if df.empty:
            return {}
        if df.shape[1] < 2:
            raise ValueError("Node dictionary DataFrame must contain at least two columns")

        node_col = df.columns[0]
        symbol_col = df.columns[1]
        cleaned = {}
        for _, row in df.iterrows():
            node = str(row[node_col]).strip()
            val = row[symbol_col]
            if pd.isna(val):
                cleaned[node] = []
                continue
            if isinstance(val, str) and val.strip().startswith("["):
                try:
                    parsed = ast.literal_eval(val)
                    cleaned[node] = normalize_symbol_list(parsed)
                    continue
                except Exception:
                    pass
            if isinstance(val, str) and "," in val:
                cleaned[node] = normalize_symbol_list([part.strip().strip("'\"") for part in val.split(",")])
            else:
                cleaned[node] = normalize_symbol_list([val])
        return cleaned

    # File path / path-like input
    path = str(node_dict_input)
    try:
        df = pd.read_csv(path)
        if df.shape[1] >= 2:
            return load_node_dict(df, verbose=verbose)
    except Exception:
        pass

    # Fall back to the shared CSV loader for historical node_dict.csv layouts.
    try:
        loaded = load_node_dict_from_csv(path, verbose=verbose)
    except Exception as exc:
        raise ValueError(f"Failed to load node_dict from {path}: {exc}") from exc

    cleaned = {}
    for node, symbols in loaded.items():
        if isinstance(symbols, (list, tuple, set)):
            cleaned[str(node)] = normalize_symbol_list(symbols)
        else:
            cleaned[str(node)] = normalize_symbol_list([symbols])
    return cleaned


def reverse_node_dict(node_dict: Mapping[str, Sequence[object]], normalize: bool = True) -> Dict[str, List[str]]:
    """Reverse node_dict into symbol -> [node_name, ...]."""
    reverse: Dict[str, List[str]] = {}
    for node, symbols in (node_dict or {}).items():
        for symbol in symbols or []:
            key = normalize_symbol(symbol) if normalize else str(symbol).strip()
            if not key:
                continue
            reverse.setdefault(key, []).append(str(node))
    return reverse


def gene_to_node_map(node_dict: Mapping[str, Sequence[object]], normalize: bool = True) -> Dict[str, List[str]]:
    """Alias for reverse_node_dict for readability at the call site."""
    return reverse_node_dict(node_dict, normalize=normalize)


def node_for_symbol(symbol: object, node_dict: Mapping[str, Sequence[object]], normalize: bool = True) -> List[str]:
    """Return all nodes mapped to a single symbol."""
    reverse = reverse_node_dict(node_dict, normalize=normalize)
    key = normalize_symbol(symbol) if normalize else str(symbol).strip()
    return reverse.get(key, [])


def map_symbols_to_nodes(symbols: Iterable[object], node_dict: Mapping[str, Sequence[object]], normalize: bool = True) -> Dict[str, List[str]]:
    """Return a symbol -> nodes mapping for a collection of symbols."""
    reverse = reverse_node_dict(node_dict, normalize=normalize)
    out: Dict[str, List[str]] = {}
    for symbol in symbols:
        key = normalize_symbol(symbol) if normalize else str(symbol).strip()
        if key:
            out[key] = reverse.get(key, [])
    return out


__all__ = [
    "gene_to_node_map",
    "load_node_dict",
    "map_symbols_to_nodes",
    "node_for_symbol",
    "normalize_symbol",
    "normalize_symbol_list",
    "reverse_node_dict",
]
