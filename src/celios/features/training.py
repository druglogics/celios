"""
Orchestration layer for building activity matrices.

This file provides a minimal orchestration surface that delegates
all heavy-lifting to the feature helper modules. Public API is kept
stable: `extract_omics` and `resolve_cell_lines`.
"""
import logging
from typing import List, Dict, Optional

import os
import pandas as pd
import sys
import time

from celios.utils import activitymatrix_report, save_file, add_resolution_report
from celios.features.cline import resolve_cell_lines as cline_resolve
from celios.features.node_mapping import load_node_dict as load_node_dict_helper, reverse_node_dict
from celios.features.expression_builder import (
    load_expression_matrix,
    prepare_expression_matrix,
    normalize_expression_matrix,
)
from celios.features.binary_builder import load_binary_matrix, aggregate_genes_to_nodes
from celios.features.tf_builder import load_tf_matrix
from celios.features.master_matrix import assemble_master_matrix
from celios.features.source_selector import select_from_master

logger = logging.getLogger(__name__)


class ActivityMatrix:
    """Lightweight container for pipeline configuration and orchestration.

    Stores input paths and options; all processing is delegated to helper
    modules. This preserves the previous public API while keeping the file
    compact and easy to review.
    """

    def __init__(
        self,
        activity_file: Optional[str] = None,
        cell_line_file: Optional[str] = None,
        node_dict: Optional[Dict[str, List[str]]] = None,
        tf_activity_file: Optional[str] = None,
        mutations_file: Optional[str] = None,
        cnv_file: Optional[str] = None,
        directory_output: Optional[str] = None,
        verbose: bool = False,
        data_sources: List[str] = None,
        format_override: Optional[str] = None,
        mutations_format_override: Optional[str] = None,
        cnv_format_override: Optional[str] = None,
        pre_resolved_cell_lines: Optional[Dict[str, object]] = None,
    ):
        # inputs
        self.activity_file = activity_file
        self.cell_line_file = cell_line_file
        self.node_dict_input = node_dict
        self.node_dict = node_dict if isinstance(node_dict, dict) else {}
        self.node_dict_reversed = None

        self.tf_activity_file = tf_activity_file
        self.mutations_file = mutations_file
        self.cnv_file = cnv_file

        # options
        self.directory_output = directory_output
        self.verbose = verbose
        self.format_override = format_override
        self.mutations_format_override = mutations_format_override
        self.cnv_format_override = cnv_format_override
        self.data_sources = data_sources or ["mutations", "cnv", "TF", "expression"]

        # runtime state
        self.activity_raw_df: Optional[pd.DataFrame] = None
        self.activity_df: Optional[pd.DataFrame] = None
        self.activity_normalized: Optional[pd.DataFrame] = None

        self.sidm_list: Optional[List[str]] = None
        self.sidm_dict: Optional[Dict[str, str]] = None
        self.alias_to_sidm: Optional[Dict[str, str]] = None

        self.activity_matrix: Optional[pd.DataFrame] = None
        self.format_metadata: Optional[Dict] = None

        if isinstance(pre_resolved_cell_lines, dict):
            self.sidm_list = pre_resolved_cell_lines.get("sidm_list")
            self.sidm_dict = pre_resolved_cell_lines.get("sidm_dict")
            self.alias_to_sidm = pre_resolved_cell_lines.get("alias_to_sidm")

        # small tunables kept here for backward compatibility
        self.p_value_tf_threshold = 0.05

    # --- logging helper -------------------------------------------------
    def _log(self, *args, level="info"):
        if not self.verbose:
            return
        if level == "info":
            logger.info(*args)
        else:
            logger.debug(*args)

    # --- cell-line resolution (thin wrapper around features.cline) -------
    def _ensure_sidm(self):
        """Ensure SIDM mapping exists; delegate to features.cline.resolve_cell_lines."""
        if self.sidm_list is not None and self.sidm_dict is not None:
            if self.verbose:
                print(f"[STEP 2] Cell line ID resolution already available: {len(self.sidm_list)} SIDMs")
            return self.sidm_list, self.sidm_dict

        if not self.cell_line_file:
            raise ValueError("cell_line_file must be provided to extract SIDM list")

        sidm_list, sidm_dict, alias_map, resolution_report = cline_resolve(cell_line_file=self.cell_line_file, verbose=self.verbose)
        self.sidm_list = sidm_list
        self.sidm_dict = sidm_dict
        self.alias_to_sidm = alias_map or {}
        add_resolution_report(resolution_report)
        if self.verbose:
            print(f"[STEP 2] Resolved SIDMs: {len(self.sidm_list)}")
        return self.sidm_list, self.sidm_dict

    def resolve_cell_line_resolution(self):
        return self._ensure_sidm()

    # --- node dict helpers (delegates to node_mapping) ------------------
    def load_node_dict(self, node_dict_input):
        return load_node_dict_helper(node_dict_input, verbose=self.verbose)

    def _ensure_node_dict(self):
        if self.node_dict:
            return self.node_dict
        if not hasattr(self, "node_dict_input") or self.node_dict_input is None:
            raise ValueError("node_dict must be provided to map genes to nodes")
        self.node_dict = self.load_node_dict(self.node_dict_input)
        return self.node_dict

    def _reverse_node_dict(self):
        if self.node_dict_reversed is not None:
            return self.node_dict_reversed
        self.node_dict_reversed = reverse_node_dict(self.node_dict)
        return self.node_dict_reversed

    # --- expression helpers (thin shims) --------------------------------
    def _load_activity_raw(self):
        if not self.activity_file:
            raise ValueError("No activity_file provided")
        df, metadata = load_expression_matrix(self.activity_file, format_override=self.format_override, verbose=self.verbose)
        self.format_metadata = metadata
        self.activity_raw_df = df
        return self.activity_raw_df

    def _prepare_activity_df(self):
        if not self.activity_file:
            raise ValueError("activity_file is required when 'expression' is in data_sources")
        self._ensure_sidm()
        df, metadata = prepare_expression_matrix(
            activity_file=self.activity_file,
            activity_df=self.activity_raw_df,
            sidm_list=self.sidm_list,
            alias_map=self.alias_to_sidm,
            format_override=self.format_override,
            verbose=self.verbose,
        )
        if self.activity_raw_df is None:
            self.activity_raw_df = df.copy()
        self.format_metadata = metadata or self.format_metadata
        self.activity_df = df
        return self.activity_df

    def _normalize_expression(self):
        if self.activity_df is None:
            self._prepare_activity_df()
        df = normalize_expression_matrix(self.activity_df, verbose=self.verbose)
        self.activity_normalized = df
        return self.activity_normalized

    # --- binary / TF helpers (delegation) --------------------------------
    def _load_binary_table(self, filepath: str, format_override: Optional[str] = None) -> pd.DataFrame:
        if not filepath:
            return None
        if self.sidm_list is None:
            self._ensure_sidm()
        df, metadata = load_binary_matrix(
            filepath,
            format_override=format_override,
            model_registry=None,
            alias_map=self.alias_to_sidm,
            sidm_list=self.sidm_list,
            collapse_method='max',
            verbose=self.verbose,
        )
        cols = [c for c in (df.columns if df is not None else []) if c in (self.sidm_list or [])]
        if df is not None and cols:
            df = df.loc[:, cols]
        return df

    def _tf_node_matrix(self):
        if not self.tf_activity_file:
            return None
        if self.sidm_list is None:
            self._ensure_sidm()
        if self.node_dict_reversed is None:
            self._reverse_node_dict()
        node_tf = load_tf_matrix(
            tf_activity_file=self.tf_activity_file,
            node_dict_reversed=self.node_dict_reversed,
            sidm_list=self.sidm_list,
            sidm_dict=self.sidm_dict,
            alias_map=self.alias_to_sidm,
            p_value_threshold=self.p_value_tf_threshold,
            binary_threshold=getattr(self, 'negative_tf_activity_threshold', 0),
            verbose=self.verbose,
        )
        return node_tf

    # --- master assembly & selection (orchestration) ---------------------
    def _node_expression_matrix(self):
        if self.activity_normalized is None:
            self._normalize_expression()
        if not self.node_dict:
            raise ValueError("node_dict must be provided to map genes to nodes")
        node_expr = aggregate_genes_to_nodes(self.activity_normalized, node_dict=self.node_dict, agg_method="mean")
        return node_expr

    def _build_master_matrix(self, data_sources: List[str] = None) -> pd.DataFrame:
        if data_sources is not None:
            self.data_sources = data_sources
        # ensure prerequisites
        self._ensure_sidm()
        self._ensure_node_dict()
        self._reverse_node_dict()

        source_matrices = {}
        if 'expression' in self.data_sources:
            source_matrices['expression'] = self._node_expression_matrix()
        if 'TF' in self.data_sources:
            source_matrices['TF'] = self._tf_node_matrix()
        if 'mutations' in self.data_sources:
            source_matrices['mutations'] = None
            if self.mutations_file:
                muts = self._load_binary_table(self.mutations_file, format_override=self.mutations_format_override)
                source_matrices['mutations'] = aggregate_genes_to_nodes(muts, node_dict=self.node_dict, agg_method="max")
        if 'cnv' in self.data_sources:
            source_matrices['cnv'] = None
            if self.cnv_file:
                cnv = self._load_binary_table(self.cnv_file, format_override=self.cnv_format_override)
                source_matrices['cnv'] = aggregate_genes_to_nodes(cnv, node_dict=self.node_dict, agg_method="max")

        master = assemble_master_matrix(
            node_dict=self.node_dict,
            sidm_list=self.sidm_list,
            source_matrices=source_matrices,
            node_index=list(self.node_dict.keys()),
            include_symbol=True,
        )
        self.activity_matrix = master
        return master

    def _select_from_master(self, master: pd.DataFrame = None, selected_sources: List[str] = None) -> pd.DataFrame:
        if master is None:
            master = self.activity_matrix
        return select_from_master(
            master=master,
            sidm_list=self.sidm_list,
            selected_sources=selected_sources or self.data_sources,
            sidm_dict=self.sidm_dict,
            include_symbol=True,
            verbose=self.verbose,
        )

    # --- pipeline orchestration -----------------------------------------
    def _run_pipeline(self, directory_output: Optional[str] = None, selected_sources: List[str] = None, save_master: bool = True, make_report: bool = True) -> pd.DataFrame:
        # directory override
        if directory_output is not None:
            self.directory_output = directory_output

        if self.directory_output is None and (save_master or make_report):
            raise ValueError("directory_output must be provided to save outputs or reports")
        if self.directory_output is not None:
            os.makedirs(self.directory_output, exist_ok=True)

        start_time = time.perf_counter()
        cmd = ' '.join(sys.argv)

        # ensure SIDM and node_dict
        try:
            self._ensure_sidm()
        except Exception:
            logger.exception('Failed to ensure SIDM mapping')
        try:
            self._ensure_node_dict()
        except Exception:
            logger.exception('Failed to ensure node dictionary')

        # build requested node-level matrices
        node_expr = node_tf = node_muts = node_cnv = None

        if 'expression' in (selected_sources or self.data_sources):
            try:
                self._load_activity_raw()
                self._prepare_activity_df()
                self._normalize_expression()
                self._reverse_node_dict()
                node_expr = self._node_expression_matrix()
            except Exception:
                logger.exception('Error preparing expression data')

        if 'mutations' in (selected_sources or self.data_sources) and self.mutations_file:
            try:
                muts = self._load_binary_table(self.mutations_file, format_override=self.mutations_format_override)
                node_muts = aggregate_genes_to_nodes(muts, node_dict=self.node_dict, agg_method="max")
            except Exception:
                logger.exception('Error loading mutations')

        if 'cnv' in (selected_sources or self.data_sources) and self.cnv_file:
            try:
                cnv = self._load_binary_table(self.cnv_file, format_override=self.cnv_format_override)
                node_cnv = aggregate_genes_to_nodes(cnv, node_dict=self.node_dict, agg_method="max")
            except Exception:
                logger.exception('Error loading CNV')

        if 'TF' in (selected_sources or self.data_sources) and self.tf_activity_file:
            try:
                node_tf = self._tf_node_matrix()
            except Exception:
                logger.exception('Error creating TF matrix')

        source_matrices = {
            'expression': node_expr,
            'TF': node_tf,
            'mutations': node_muts,
            'cnv': node_cnv,
        }

        master = assemble_master_matrix(
            node_dict=self.node_dict,
            sidm_list=self.sidm_list,
            source_matrices=source_matrices,
            node_index=list(self.node_dict.keys()),
            include_symbol=True,
        )

        master_fp = None
        if save_master:
            master_fp = os.path.join(self.directory_output, 'activity_master_matrix.csv')
            save_file(master, self.directory_output, 'activity_master_matrix.csv', index=True)
            if self.verbose:
                logger.info('Master matrix saved: %s', master_fp)

        final = select_from_master(
            master=master,
            sidm_list=self.sidm_list,
            selected_sources=selected_sources or self.data_sources,
            sidm_dict=self.sidm_dict,
            include_symbol=True,
            verbose=self.verbose,
        )

        final_fp = os.path.join(self.directory_output, 'activity_from_master.csv') if self.directory_output else None
        if final_fp:
            save_file(final, self.directory_output, 'activity_from_master.csv', index=True)
            if self.verbose:
                logger.info('Final activity saved: %s', final_fp)

        end_time = time.perf_counter()
        runtime = end_time - start_time

        if make_report:
            report_data = {
                'selected_sources': selected_sources or self.data_sources,
                'command': cmd,
                'runtime_seconds': runtime,
                'node_dict': self.node_dict,
                'missing_hgnc_nodes': [n for n, syms in (self.node_dict or {}).items() if not syms],
                'sidm_list': self.sidm_list,
                'sidm_dict': self.sidm_dict,
                'raw_shape': getattr(self.activity_raw_df, 'shape', None),
                'prep_shape': getattr(self.activity_df, 'shape', None),
                'norm_shape': getattr(self.activity_normalized, 'shape', None),
                'node_expr_shape': getattr(node_expr, 'shape', None),
                'muts_shape': getattr(node_muts, 'shape', None),
                'cnv_shape': getattr(node_cnv, 'shape', None),
                'tfm_shape': getattr(node_tf, 'shape', None),
                'master_fp': master_fp,
                'master_shape': getattr(master, 'shape', None),
                'final_fp': final_fp,
                'final_shape': getattr(final, 'shape', None),
                'nodes_with_all_missing': [],
                'format_metadata': self.format_metadata,
            }
            activitymatrix_report(self.directory_output, report_data, verbose=self.verbose)

        return final


def extract_omics(
    activity_file: Optional[str] = None,
    cell_line_file: Optional[str] = None,
    node_dict: Optional[Dict[str, List[str]]] = None,
    tf_activity_file: Optional[str] = None,
    mutations_file: Optional[str] = None,
    cnv_file: Optional[str] = None,
    directory_output: Optional[str] = None,
    verbose: bool = False,
    data_sources: List[str] = None,
    save_master: bool = True,
    make_report: bool = True,
    format_override: Optional[str] = None,
    mutations_format_override: Optional[str] = None,
    cnv_format_override: Optional[str] = None,
    pre_resolved_cell_lines: Optional[Dict[str, object]] = None,
):
    am = ActivityMatrix(
        activity_file=activity_file,
        cell_line_file=cell_line_file,
        node_dict=node_dict,
        tf_activity_file=tf_activity_file,
        mutations_file=mutations_file,
        cnv_file=cnv_file,
        directory_output=directory_output,
        verbose=verbose,
        data_sources=data_sources,
        format_override=format_override,
        mutations_format_override=mutations_format_override,
        cnv_format_override=cnv_format_override,
        pre_resolved_cell_lines=pre_resolved_cell_lines,
    )

    # try to load node_dict if path-like
    if isinstance(node_dict, (str,)):
        try:
            am.load_node_dict(node_dict)
        except Exception:
            if verbose:
                logger.info("Failed to auto-load node_dict from path: %s", node_dict)

    return am._run_pipeline(directory_output=directory_output, selected_sources=data_sources, save_master=save_master, make_report=make_report)


def resolve_cell_lines(
    cell_line_file: Optional[str] = None,
    verbose: bool = False,
) -> Dict[str, object]:
    am = ActivityMatrix(cell_line_file=cell_line_file, verbose=verbose)
    sidm_list, sidm_dict = am.resolve_cell_line_resolution()
    if verbose:
        print(f"[STEP 2] Resolution complete: {len(sidm_list)} SIDMs, {len(am.alias_to_sidm or {})} aliases")
    return {
        "sidm_list": sidm_list,
        "sidm_dict": sidm_dict,
        "alias_to_sidm": am.alias_to_sidm or {},
    }


__all__ = ["extract_omics", "resolve_cell_lines"]
