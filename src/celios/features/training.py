"""
training.py

Omics Activity builder that produces a master matrix with rows as model nodes
and columns for each combination of cell-line (SIDM) and data source (expression, TF, mutations, cnv).

Columns are named with the pattern: <SIDM>__<source> (e.g. "SIDM123__expression").

This keeps all processed data in a single DataFrame and lets users select which data
sources to use downstream by selecting the relevant columns.

This file aims to be a minimal, clear implementation to try the different architecture
requested by the user.
"""

import logging
from typing import List, Dict, Optional

import os
import pandas as pd
import ast
import numpy as np
import sys
import time
import re

from celios.utils import activitymatrix_report, save_file, load_csv_file, add_resolution_report
from celios.utils.cell_line_resolver import resolve_sidm_from_dataframe
from celios.features.activity_parser import FormatDetector, get_parser
from celios.features.binary_parser import (
    FormatDetector as BinaryFormatDetector,
    get_binary_parser,
)

logger = logging.getLogger(__name__)


class ActivityMatrix:
    """Builds a master node x (SIDM x source) matrix.

    Args:
        activity_file: path to gene-level activity CSV (multi-header like original repo)
        cell_line_file: path to cell line file (contains 'SIDM' and 'cell_line_name')
        node_dict: mapping node_name -> list of HGNC symbols
        tf_activity_file, mutations_file, cnv_file: optional file paths
        verbose: if True, emits logger.info messages
        data_sources: list of sources to include by default
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
    ):
        self.activity_file = activity_file
        self.cell_line_file = cell_line_file
        self.node_dict = node_dict or {}
        self.node_dict_reversed = None

        self.tf_activity_file = tf_activity_file
        self.mutations_file = mutations_file
        self.cnv_file = cnv_file

        self.directory_output = directory_output
        self.verbose = verbose
        self.format_override = format_override
        self.mutations_format_override = mutations_format_override
        self.cnv_format_override = cnv_format_override

        self.activity_raw_df: Optional[pd.DataFrame] = None
        self.activity_df: Optional[pd.DataFrame] = None  # gene-level with single header
        self.activity_normalized: Optional[pd.DataFrame] = None

        self.sidm_list: Optional[List[str]] = None
        self.sidm_dict: Optional[Dict[str, str]] = None

        self.mutations_data: Optional[pd.DataFrame] = None
        self.cnv_data: Optional[pd.DataFrame] = None
        self.activity_matrix: Optional[pd.DataFrame] = None

        self.data_sources = data_sources or ["mutations", "cnv", "TF", "expression"]

        # thresholds for TF filtering (kept small & configurable)
        self.p_value_tf_threshold = 0.05

        # Store format metadata from parser
        self.format_metadata: Optional[Dict] = None

    # ------------------------------------------------------------------
    # Helpers to load basic inputs
    # ------------------------------------------------------------------
    def _log(self, *args, level="info"):
        if not self.verbose:
            return
        if level == "info":
            logger.info(*args)
        else:
            logger.debug(*args)

    def _load_activity_raw(self):
        """Load raw activity file using format auto-detection and pluggable parser.

        Supports multiple formats (old, 26Q1, etc.) with automatic detection.
        Stores the parsed DataFrame and format metadata.
        
        Returns:
            DataFrame with gene symbols as index and SIDM IDs as columns.
        """
        if not self.activity_file:
            raise ValueError("No activity_file provided")
        
        self._log("Loading activity file: %s", self.activity_file)
        
        # Auto-detect or validate format
        detected_format = FormatDetector.detect(self.activity_file, format_override=self.format_override)
        self._log("Activity file format: %s", detected_format)
        
        # Get appropriate parser and load file
        parser = get_parser(detected_format, verbose=self.verbose)
        df, format_metadata = parser.load(self.activity_file)
        
        # Store metadata
        self.format_metadata = format_metadata
        
        # Store the parsed DataFrame
        self.activity_raw_df = df
        return self.activity_raw_df

    def load_node_dict(self, node_dict_input):
        """Load or coerce a node dictionary.

        Accepts either a dict already, or a path (string/Path) to a CSV file
        with at least two columns (node_name, HGNC_symbol).
        Sets `self.node_dict` to a mapping node -> list[str].
        """
        if node_dict_input is None:
            return None

        # If already a dict, accept as-is
        if isinstance(node_dict_input, dict):
            self.node_dict = node_dict_input
            return self.node_dict

        # Otherwise assume a path
        path = str(node_dict_input)
        try:
            df = load_csv_file(path)
        except Exception:
            df = pd.read_csv(path)

        # heuristics for columns
        if 'node_name' in df.columns and ('HGNC_symbol' in df.columns or 'symbol' in df.columns):
            node_col = 'node_name'
            sym_col = 'HGNC_symbol' if 'HGNC_symbol' in df.columns else 'symbol'
        else:
            node_col = df.columns[0]
            sym_col = df.columns[1] if len(df.columns) > 1 else None

        if sym_col is None:
            raise ValueError('Node dictionary file must contain at least two columns: node and symbol(s)')

        cleaned = {}
        for _, row in df.iterrows():
            node = row[node_col]
            val = row[sym_col]
            if pd.isna(val):
                cleaned[node] = []
                continue
            if isinstance(val, list):
                cleaned[node] = [str(v).strip() for v in val if v is not None]
                continue
            if isinstance(val, str) and val.strip().startswith('['):
                try:
                    parsed = ast.literal_eval(val)
                    cleaned[node] = [str(v).strip() for v in parsed if v is not None]
                    continue
                except Exception:
                    pass
            if isinstance(val, str) and ',' in val:
                cleaned[node] = [s.strip().strip("'\"") for s in val.split(',') if s.strip()]
            else:
                cleaned[node] = [str(val).strip()]

        self.node_dict = cleaned
        return self.node_dict

    def _ensure_sidm(self):
        """Ensure `self.sidm_list` and `self.sidm_dict` are populated from `cell_line_file`.
        
          Resolution strategy:
             1. Read user table (csv/tsv auto-detected)
             2. Resolve identifiers row-wise to SIDM via local Model.csv resolver
                 supporting SIDM, ModelID(ACH), RRID/CVCL, CCLE_ID, and names.
        
        Returns:
            tuple: (sidm_list, sidm_dict) where sidm_dict maps SIDM -> cell_line_name
        
        Raises:
            ValueError: If cell_line_file is not provided or no SIDM mappings can be found.
        """
        if self.sidm_list is not None and self.sidm_dict is not None:
            return self.sidm_list, self.sidm_dict
        if not self.cell_line_file:
            raise ValueError("cell_line_file must be provided to extract SIDM list")
        self._log("Loading cell line file: %s", self.cell_line_file)

        # Auto-detect delimiter to support both .csv and .tsv cell-line files.
        try:
            df = load_csv_file(self.cell_line_file, sep=None, engine="python")
        except Exception:
            df = load_csv_file(self.cell_line_file)

        try:
            sidm_dict, not_found, resolution_report = resolve_sidm_from_dataframe(df)
            # Log the resolution report to run logs
            add_resolution_report(resolution_report)
        except Exception as e:
            raise ValueError(
                f"Failed to resolve SIDM mapping from cell_line_file: {e}. "
                "Provide at least one supported identifier column (SIDM, ModelID, RRID/CVCL, CCLE_ID, or cell_line_name)."
            ) from e

        if not sidm_dict:
            raise ValueError(
                "No cell lines from cell_line_file could be resolved to SIDM. "
                "Check identifier values or include SIDM column explicitly."
            )

        if not_found:
            warning_msg = (
                "WARNING: Could not resolve %d cell line identifier(s): %s\n"
                "These rows will be excluded from the analysis."
                % (len(not_found), ", ".join(not_found))
            )
            print("=" * 80)
            print(warning_msg)
            print("=" * 80)
            if self.verbose:
                logger.warning(warning_msg)

        self.sidm_list = list(sidm_dict.keys())
        self.sidm_dict = sidm_dict
        self._log("Resolved %s SIDM entries (excluded %s)", len(self.sidm_list), len(not_found))
        return self.sidm_list, self.sidm_dict

    # ------------------------------------------------------------------
    # Expression processing
    # ------------------------------------------------------------------
    def _prepare_activity_df(self):
        """Create a simplified gene-level DataFrame with index 'symbol' and columns=SIDM IDs.

        This method adapts the parser output (which may have different column structures
        depending on format) to a normalized form with SIDM columns.

        Sets `self.activity_df`.
        
        Raises:
            ValueError: If activity file is not available (required for expression processing).
        """
        if not self.activity_file:
            raise ValueError(
                "activity_file is required when 'expression' is in data_sources. "
                "Provide activity_file or remove 'expression' from data_sources."
            )
        
        if self.activity_raw_df is None:
            self._load_activity_raw()

        df = self.activity_raw_df.copy()

        # Parser output already has symbol as index
        # Ensure index is named 'symbol' and values are normalized (uppercase)
        if df.index.name != "symbol":
            df.index.name = "symbol"
        df.index = df.index.astype(str).str.upper()

        # ensure we have sidm_list for filtering
        self._ensure_sidm()
        
        # Filter columns for SIDMs present
        # For old format: columns are already SIDM IDs
        # For 26Q1 format: columns may be ModelID or SequencingID, but we have SIDM from cell_line_file
        cols = [c for c in df.columns if c in self.sidm_list]
        
        if not cols:
            # If no direct SIDM column match, columns might be in a different format
            # Log warning and keep all columns for now
            self._log(
                "Warning: No columns matched SIDM list. Columns in DataFrame: %s. "
                "SIDM list: %s. Keeping all columns.",
                list(df.columns)[:5],
                self.sidm_list[:5],
            )
        else:
            df = df.loc[:, cols]

        self.activity_df = df
        self._log("Prepared activity_df shape: %s", df.shape)
        return self.activity_df

    def _normalize_expression(self):
        """Log-transform and min-max normalize expression per gene (row-wise).

        Result stored in `self.activity_normalized`.
        """
        if self.activity_df is None:
            self._prepare_activity_df()
        df = self.activity_df.copy(deep=True)
        offset = 0.01
        numeric_cols = df.select_dtypes(include=[float, int]).columns
        df.loc[:, numeric_cols] = df.loc[:, numeric_cols].apply(lambda x: np.log(x + offset))
        # min-max per row
        df.loc[:, numeric_cols] = df.loc[:, numeric_cols].apply(lambda x: (x - x.min()) / (x.max() - x.min()), axis=1)
        self.activity_normalized = df
        self._log("Normalized activity shape: %s", df.shape)
        return self.activity_normalized

    # ------------------------------------------------------------------
    # Node aggregation
    # ------------------------------------------------------------------
    def _reverse_node_dict(self):
        if self.node_dict_reversed is not None:
            return self.node_dict_reversed
        rev = {}
        for node, symbols in (self.node_dict or {}).items():
            for s in symbols:
                rev.setdefault(s, []).append(node)
        self.node_dict_reversed = rev
        return rev

    def _node_expression_matrix(self):
        """Aggregate gene-level normalized expression to node-level (mean across mapped genes).

        Returns DataFrame indexed by node_name, columns = SIDM IDs.
        """
        if self.activity_normalized is None:
            self._normalize_expression()
        if not self.node_dict:
            raise ValueError("node_dict must be provided to map genes to nodes")

        # prepare a dataframe with symbol as column and sidm columns
        df = self.activity_normalized.copy()
        # map each gene (symbol) to node_name(s) and explode
        mapping = pd.Series({s: self.node_dict_reversed.get(s, []) for s in df.index})
        # Turn mapping into rows: for each symbol, list of nodes
        # Create a DataFrame with index symbol and column node_name (exploded)
        exploded = (
            mapping.to_frame('node_list')
            .reset_index()
            .rename(columns={'index': 'symbol'})
            .explode('node_list')
            .dropna(subset=['node_list'])
        )
        if exploded.empty:
            raise ValueError("No gene symbols in activity map to provided node_dict")

        exploded = exploded.set_index('symbol')
        # join expression values to node mapping and aggregate mean by node
        joined = df.join(exploded['node_list'])
        joined = joined.reset_index().rename(columns={'node_list': 'node_name'})
        # now group by node_name and take mean across numeric columns only
        node_expr = joined.groupby('node_name').mean(numeric_only=True)
        self._log("Node expression matrix shape: %s", node_expr.shape)
        return node_expr

    # ------------------------------------------------------------------
    # Mutations / CNV processing
    # ------------------------------------------------------------------
    def _load_binary_table(
        self,
        filepath: str,
        format_override: Optional[str] = None,
        index_col: str = "gene_symbol",
    ) -> pd.DataFrame:
        """Load a binary matrix (mutations, CNV) with format auto-detection.

        Supports two formats:
        - OLD: genes × SIDM
        - 26Q1: ModelID × genes (transposed + mapped to SIDM)

        Args:
            filepath: Path to binary matrix file
            format_override: Force format ("old" | "26q1" | None for auto-detect)
            index_col: Column name for index (typically 'gene_symbol')

        Returns:
            DataFrame with genes as index, SIDM columns, binary values
        """
        if not filepath:
            return None

        # Auto-detect format
        try:
            detected_format = BinaryFormatDetector.detect(filepath, format_override)
        except ValueError as e:
            # Fall back to old behavior if detection fails
            self._log(f"Binary format detection failed: {e}. Using old format.", level="info")
            df = load_csv_file(filepath, index_col=0)
            if "gene_symbol" in df.columns:
                df = df.set_index("gene_symbol")
            if self.sidm_list is None:
                self._ensure_sidm()
            cols = [c for c in df.columns if c in self.sidm_list]
            df = df.loc[:, cols]
            return df

        # Use binary parser factory
        parser = get_binary_parser(detected_format, verbose=self.verbose)
        df, metadata = parser.load(filepath)

        # Filter to SIDM list
        if self.sidm_list is None:
            self._ensure_sidm()
        cols = [c for c in df.columns if c in self.sidm_list]
        df = df.loc[:, cols]

        self._log(f"Loaded {detected_format} binary matrix: {metadata['shape']} genes×samples")

        return df

    def _node_binary_matrix(self, bin_df: pd.DataFrame) -> pd.DataFrame:
        """Aggregate gene-level binary matrix to node-level by taking column-wise max across genes mapped to a node.

        Returns DataFrame indexed by node_name, columns = SIDMs, values are 0/1/NaN.
        """
        if bin_df is None:
            return None
        if not self.node_dict:
            raise ValueError("node_dict required for node-level binary aggregation")

        rows = []
        for node, genes in self.node_dict.items():
            # select rows for genes intersection
            sub = bin_df.loc[bin_df.index.isin(genes)]
            if sub.empty:
                vals = pd.Series(index=bin_df.columns, dtype=float)
            else:
                # take maximum across genes to indicate presence
                vals = sub.max(axis=0, skipna=True)
            vals.name = node
            rows.append(vals)
        node_bin = pd.DataFrame(rows)
        node_bin.index = list(self.node_dict.keys())
        return node_bin

    # ------------------------------------------------------------------
    # TF processing (simple mapping)
    # ------------------------------------------------------------------
    def _tf_node_matrix(self):
        """Load TF activity file and aggregate to node-level. Expects a file with columns
        including at least ['source', 'condition', 'p_value', 'score'] where 'condition'
        corresponds to the cell line name (we map to SIDM using sidm_dict).
        """
        if not self.tf_activity_file:
            return None
        if self.sidm_dict is None:
            self._ensure_sidm()
        # load TF activity (use pandas to be permissive)
        tf = pd.read_csv(self.tf_activity_file, sep=r'\s+')
        # filter by p_value
        if 'p_value' in tf.columns:
            tf = tf[tf['p_value'] < self.p_value_tf_threshold]
        # map condition (cell line name) to sidm
        # create reversed mapping cell_line_name -> SIDM
        rev = {v: k for k, v in self.sidm_dict.items()}
        tf = tf[tf['condition'].isin(rev.keys())]
        tf['sidm'] = tf['condition'].map(rev)

        # map TF source -> node_name using node_dict_reversed
        if self.node_dict_reversed is None:
            self._reverse_node_dict()
        tf['node_name'] = tf['source'].map(self.node_dict_reversed)
        # explode node_name if lists
        tf = tf.explode('node_name').dropna(subset=['node_name'])
        # compute node activity as max score per node x sidm
        pivot = tf.pivot_table(index='node_name', columns='sidm', values='score', aggfunc='max')
        # ensure columns are full SIDM set
        for sidm in self.sidm_list:
            if sidm not in pivot.columns:
                pivot[sidm] = np.nan
        pivot = pivot.loc[:, self.sidm_list]
        # Convert TF scores to binary activity (1 if score > threshold, 0 if score <= threshold)
        # Preserve NaN where no TF score is available.
        thresh = getattr(self, 'negative_tf_activity_threshold', 0)
        def to_binary(v):
            if pd.isna(v):
                return np.nan
            return 1.0 if v > thresh else 0.0

        # Vectorized conversion to binary while preserving NaNs (avoid deprecated applymap)
        binary = (pivot > thresh).astype(float)
        node_activity = pivot.where(pivot.isna(), binary)
        return node_activity

    # ------------------------------------------------------------------
    # Master matrix builder
    # ------------------------------------------------------------------
    def _build_master_matrix(self, data_sources: List[str] = None) -> pd.DataFrame:
        """Build and return the master DataFrame.

        Columns are named `<SIDM>__<source>`. The DataFrame includes a 'symbol' column
        listing the node's mapped HGNC symbols.
        """
        if data_sources is not None:
            self.data_sources = data_sources
        if self.verbose:
            self._log("Building master activity matrix with sources: %s", self.data_sources)

        # ensure sidm and only prepare activity if expression is needed
        if 'expression' in self.data_sources:
            self._prepare_activity_df()
        
        self._reverse_node_dict()

        # node-level symbol list column
        node_symbols = {node: syms for node, syms in (self.node_dict or {}).items()}

        # prepare containers
        node_index = list(self.node_dict.keys())
        master = pd.DataFrame(index=node_index)
        master['symbol'] = master.index.map(node_symbols)

        # build expression node matrix if requested
        if 'expression' in self.data_sources:
            node_expr = self._node_expression_matrix()
            # ensure SIDM columns present
            expr_data = {}
            for sidm in self.sidm_list:
                colname = f"{sidm}__expression"
                # node_expr may not have all nodes (if no mapped genes) -> reindex
                vals = node_expr.reindex(master.index).get(sidm)
                expr_data[colname] = vals.values
            master = pd.concat([master, pd.DataFrame(expr_data, index=master.index)], axis=1)

        # build TF node matrix if requested
        if 'TF' in self.data_sources:
            node_tf = self._tf_node_matrix()
            if node_tf is not None:
                tf_data = {}
                for sidm in self.sidm_list:
                    colname = f"{sidm}__TF"
                    vals = node_tf.reindex(master.index).get(sidm)
                    tf_data[colname] = vals.values
                master = pd.concat([master, pd.DataFrame(tf_data, index=master.index)], axis=1)
            else:
                tf_data = {f"{sidm}__TF": np.nan for sidm in self.sidm_list}
                master = pd.concat([master, pd.DataFrame(tf_data, index=master.index)], axis=1)

        # build mutations/node binary
        if 'mutations' in self.data_sources:
            if self.mutations_file:
                muts = self._load_binary_table(
                    self.mutations_file, format_override=self.mutations_format_override
                )
                node_muts = self._node_binary_matrix(muts)
                mutations_data = {}
                for sidm in self.sidm_list:
                    colname = f"{sidm}__mutations"
                    vals = node_muts.reindex(master.index).get(sidm) if node_muts is not None else None
                    mutations_data[colname] = vals.values if vals is not None else np.nan
                master = pd.concat([master, pd.DataFrame(mutations_data, index=master.index)], axis=1)
            else:
                mutations_data = {f"{sidm}__mutations": np.nan for sidm in self.sidm_list}
                master = pd.concat([master, pd.DataFrame(mutations_data, index=master.index)], axis=1)

        # build CNV node binary
        if 'cnv' in self.data_sources:
            if self.cnv_file:
                cnv = self._load_binary_table(
                    self.cnv_file, format_override=self.cnv_format_override
                )
                node_cnv = self._node_binary_matrix(cnv)
                cnv_data = {}
                for sidm in self.sidm_list:
                    colname = f"{sidm}__cnv"
                    vals = node_cnv.reindex(master.index).get(sidm) if node_cnv is not None else None
                    cnv_data[colname] = vals.values if vals is not None else np.nan
                master = pd.concat([master, pd.DataFrame(cnv_data, index=master.index)], axis=1)
            else:
                cnv_data = {f"{sidm}__cnv": np.nan for sidm in self.sidm_list}
                master = pd.concat([master, pd.DataFrame(cnv_data, index=master.index)], axis=1)

        # optional: save
        self.activity_matrix = master
        self._log("Master matrix built with shape %s", master.shape)
        return master

    def _select_from_master(self, master: pd.DataFrame = None, selected_sources: List[str] = None) -> pd.DataFrame:
        """Apply selection logic over a master matrix to produce final node x SIDM activity.

        Priority per node×SIDM: mutations -> cnv -> TF -> expression. The `master` DataFrame
        is expected to have columns named `<SIDM>__<source>` (e.g. `SIDM001__mutations`).

        `selected_sources` filters which sources are considered (subset of ['mutations','cnv','TF','expression']).
        Returns a DataFrame indexed by node_name with columns = SIDM IDs.
        """
        if master is None:
            if self.activity_matrix is None:
                raise ValueError("master DataFrame must be provided or built via _build_master_matrix()")
            master = self.activity_matrix

        if selected_sources is None:
            selected_sources = self.data_sources

        # result DataFrame: index = nodes, columns = SIDMs
        nodes = master.index
        sidms = self.sidm_list
        result = pd.DataFrame(index=nodes, columns=sidms, dtype=float)

        for sidm in sidms:
            # column names
            col_mut = f"{sidm}__mutations"
            col_cnv = f"{sidm}__cnv"
            col_tf = f"{sidm}__TF"
            col_expr = f"{sidm}__expression"

            # start with NaNs
            out = pd.Series(index=nodes, dtype=float)

            # mutations
            if 'mutations' in selected_sources and col_mut in master.columns:
                muts = master[col_mut]
                mask_m = muts.isin([0, 1])
                out[mask_m] = muts[mask_m].astype(float)

            # cnv
            if 'cnv' in selected_sources and col_cnv in master.columns:
                cnv = master[col_cnv]
                mask_c = cnv.isin([0, 1])
                need = out.isna()
                assign = mask_c & need
                out[assign] = cnv[assign].astype(float)

            # TF
            if 'TF' in selected_sources and col_tf in master.columns:
                tf = master[col_tf]
                mask_t = tf.notna()
                need = out.isna()
                assign = mask_t & need
                out[assign] = tf[assign].astype(float)

            # expression
            if 'expression' in selected_sources and col_expr in master.columns:
                expr = master[col_expr]
                mask_e = expr.notna()
                need = out.isna()
                assign = mask_e & need
                out[assign] = expr[assign].astype(float)

            result[sidm] = out

        # attach symbol column from master if present
        if 'symbol' in master.columns:
            result = result.reset_index().rename(columns={'index': 'node_name'})
            result.insert(1, 'symbol', master['symbol'].reindex(result['node_name']).values)
            result = result.set_index('node_name')

        # optional save
        self._log("Selected final activity matrix shape: %s", result.shape)

        # Map SIDM column names to cell line names if mapping available
        if self.sidm_dict is None:
            try:
                self._ensure_sidm()
            except Exception:
                # if we can't build sidm mapping, leave as SIDM ids
                pass

        if self.sidm_dict is not None:
            # build mapping for columns that are SIDM IDs
            rename_map = {sidm: self.sidm_dict.get(sidm, sidm) for sidm in result.columns if sidm in self.sidm_dict}
            if rename_map:
                result = result.rename(columns=rename_map)

        return result

    def _run_pipeline(self, directory_output: Optional[str] = None, selected_sources: List[str] = None, save_master: bool = True, make_report: bool = True) -> pd.DataFrame:
        """Run the full master-matrix pipeline and produce outputs + a run log.

        This orchestrates: prepare inputs, build master matrix (optionally save it),
        select the final activity matrix using `selected_sources`, save the results,
        and write a short run log into `directory_output/run_log.txt`.
        Returns the final selected DataFrame.
        """
        # allow override of directory output
        if directory_output is not None:
            self.directory_output = directory_output

        # If no directory provided, operate in-memory and skip saving.
        # Only require a directory when saving outputs or reports is requested.
        if self.directory_output is None:
            if save_master or make_report:
                raise ValueError("directory_output must be provided to save outputs or reports")
        else:
            os.makedirs(self.directory_output, exist_ok=True)

        start_time = time.perf_counter()
        cmd = ' '.join(sys.argv)

        # simple preparation stage
        try:
            self._ensure_sidm()
        except Exception:
            if self.verbose:
                logger.warning('Failed to ensure SIDM mapping')

        raw = prep = norm = node_expr = muts = cnv = tfm = None

        # load/prepare activity only if expression is needed
        if 'expression' in (selected_sources or self.data_sources):
            try:
                raw = self._load_activity_raw()
                prep = self._prepare_activity_df()
                norm = self._normalize_expression()
                self._reverse_node_dict()
                node_expr = self._node_expression_matrix()
            except Exception as e:
                if self.verbose:
                    logger.exception('Error preparing expression data: %s', e)

        # load binary tables and TF matrix if requested
        try:
            if 'mutations' in (selected_sources or self.data_sources) and self.mutations_file:
                muts = self._load_binary_table(
                    self.mutations_file, format_override=self.mutations_format_override
                )
        except Exception as e:
            if self.verbose:
                logger.exception('Error loading mutations: %s', e)

        try:
            if 'cnv' in (selected_sources or self.data_sources) and self.cnv_file:
                cnv = self._load_binary_table(
                    self.cnv_file, format_override=self.cnv_format_override
                )
        except Exception as e:
            if self.verbose:
                logger.exception('Error loading CNV: %s', e)

        try:
            if 'TF' in (selected_sources or self.data_sources) and self.tf_activity_file:
                tfm = self._tf_node_matrix()
        except Exception as e:
            if self.verbose:
                logger.exception('Error creating TF matrix: %s', e)

        # build master
        master = self._build_master_matrix(data_sources=selected_sources)
        master_fp = None
        if save_master:
            master_fp = os.path.join(self.directory_output, 'activity_master_matrix.csv')
            save_file(master, self.directory_output, 'activity_master_matrix.csv', index=True)
            if self.verbose:
                logger.info('Master matrix saved: %s', master_fp)

        # select final
        final = self._select_from_master(master, selected_sources=selected_sources)
        final_fp = os.path.join(self.directory_output, 'activity_from_master.csv')
        save_file(final, self.directory_output, 'activity_from_master.csv', index=True)
        if self.verbose:
            logger.info('Final activity saved: %s', final_fp)

        # diagnostics
        nan_counts = final.isna().sum(axis=1)
        nodes_with_all_missing = nan_counts[nan_counts == final.shape[1]].index.tolist()

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
                'raw_shape': getattr(raw, 'shape', None),
                'prep_shape': getattr(prep, 'shape', None),
                'norm_shape': getattr(norm, 'shape', None),
                'node_expr_shape': getattr(node_expr, 'shape', None),
                'muts_shape': getattr(muts, 'shape', None),
                'cnv_shape': getattr(cnv, 'shape', None),
                'tfm_shape': getattr(tfm, 'shape', None),
                'master_fp': master_fp,
                'master_shape': getattr(master, 'shape', None),
                'final_fp': final_fp,
                'final_shape': getattr(final, 'shape', None),
                'nodes_with_all_missing': nodes_with_all_missing,
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
):
    """Public convenience function that runs the ActivityMatrix pipeline.

    This is the only symbol exported by this module (see `__all__`).
    It constructs an `ActivityMatrix` and runs the pipeline, returning the final
    selected activity DataFrame.

    Args:
        format_override: Optional format hint for activity ('old' or '26Q1'). If None, auto-detects.
        mutations_format_override: Optional format hint for mutations ('old' or '26q1'). If None, auto-detects.
        cnv_format_override: Optional format hint for CNV ('old' or '26q1'). If None, auto-detects.
    """
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
    )

    # If node_dict was provided as a path, attempt to load it via helper
    if isinstance(node_dict, (str,)):
        try:
            am.load_node_dict(node_dict)
        except Exception:
            # fall back silently; run_pipeline will report missing mapping
            if verbose:
                logger.info("Failed to auto-load node_dict from path: %s", node_dict)

    # run internal pipeline; keep API limited to this convenience function
    return am._run_pipeline(directory_output=directory_output, selected_sources=data_sources, save_master=save_master, make_report=make_report)


__all__ = ["extract_omics"]
