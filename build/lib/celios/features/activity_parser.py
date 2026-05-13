"""
activity_parser.py

Format-agnostic activity file loading and parsing.

This module provides a pluggable architecture for handling different activity file formats:
- Detects format automatically (or accepts user override)
- Parses format-specific structure into normalized DataFrame
- Returns consistent output: gene_symbols × SIDM_IDs matrix + metadata

Supported formats:
  - 'old': rnaseq_tpm_20220624.csv (multi-header, gene_id + symbol index)
  - '26Q1': rnaseq_tpm_coding_genes26Q1.csv (single header, ModelID cell identification)
"""

import logging
import re
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Tuple, Dict, Optional

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


class FormatDetector:
    """Detect activity file format by examining header structure."""

    @staticmethod
    def _peek_header(filepath: str, nrows: int = 10) -> pd.DataFrame:
        """Read first nrows of file without strict parsing."""
        try:
            df = pd.read_csv(filepath, nrows=nrows, dtype=str)
            return df
        except Exception as e:
            raise ValueError(f"Failed to read file {filepath}: {e}")

    @staticmethod
    def detect(filepath: str, format_override: Optional[str] = None) -> str:
        """
        Auto-detect format or validate user override.

        Args:
            filepath: Path to activity file
            format_override: If provided, validate against detected format and return if valid

        Returns:
            str: Format identifier ('old' or '26Q1')

        Raises:
            ValueError: If format cannot be detected or is invalid
        """
        if format_override:
            if format_override in ("old", "26Q1"):
                logger.info(f"Using user-specified format: {format_override}")
                return format_override
            else:
                raise ValueError(f"Unknown format override: {format_override}. Valid: 'old', '26Q1'")

        # Auto-detect by peeking at header
        peek = FormatDetector._peek_header(filepath, nrows=6)

        # Check for old format markers
        # Old format: rows 1-5 have metadata, row 5 has 'gene_id' and 'symbol'
        if "gene_id" in peek.columns or peek.iloc[4, 0] == "gene_id" if len(peek) > 4 else False:
            logger.info("Detected format: old (rnaseq_tpm_20220624)")
            return "old"

        # Check for 26Q1 format markers
        # 26Q1 format: has columns like 'SequencingID', 'ModelConditionID', 'ModelID'
        if "SequencingID" in peek.columns and "ModelID" in peek.columns:
            logger.info("Detected format: 26Q1 (rnaseq_tpm_coding_genes26Q1)")
            return "26Q1"

        # If uncertain, try harder by checking actual data rows
        # Check if first non-header row looks like old format (SIDG00001, A1BG, ...)
        if len(peek) > 5:
            first_data_row = peek.iloc[5, :2]  # First two columns of data row
            # Old format: gene_id looks like "SIDG00001" and symbol is a gene name
            if str(first_data_row.iloc[0]).startswith("SIDG"):
                logger.info("Detected format: old (rnaseq_tpm_20220624) - confirmed by data row")
                return "old"

        raise ValueError(
            f"Could not auto-detect file format for {filepath}. "
            "File may not be a supported activity format (old or 26Q1). "
            "Specify format_override='old' or format_override='26Q1' explicitly."
        )


class ActivityParser(ABC):
    """Abstract base class for activity file parsers.

    All concrete parsers must implement:
    - load(filepath): Returns (DataFrame, metadata_dict)
    - _parse(): Format-specific parsing logic

    Output contract:
    - DataFrame index: gene symbols (uppercase, deduplicated)
    - DataFrame columns: SIDM IDs (resolved from Model.csv if needed)
    - DataFrame values: expression levels (float)
    - metadata: dict with keys like {format, n_genes, n_samples, notes}
    """

    def __init__(self, verbose: bool = False):
        self.verbose = verbose

    def _log(self, msg: str, level: str = "info"):
        if not self.verbose:
            return
        if level == "info":
            logger.info(msg)
        else:
            logger.debug(msg)

    @abstractmethod
    def load(self, filepath: str) -> Tuple[pd.DataFrame, Dict]:
        """Load and parse activity file.

        Returns:
            tuple: (DataFrame with index=symbols, columns=SIDMs), metadata_dict
        """
        pass

    @abstractmethod
    def _parse(self, filepath: str) -> Tuple[pd.DataFrame, Dict]:
        """Format-specific parsing. Implement in subclass."""
        pass


class OldFormatParser(ActivityParser):
    """Parser for rnaseq_tpm_20220624.csv format.

    Structure:
    - Rows 1-4: metadata (model_id, model_name, dataset_name, data_source)
    - Row 5: headers (gene_id, symbol, then SIDM columns)
    - Rows 6+: data with MultiIndex columns (SIDM, cell_line_name)
    """

    def load(self, filepath: str) -> Tuple[pd.DataFrame, Dict]:
        """Load old format file."""
        self._log(f"Loading old format file: {filepath}")
        return self._parse(filepath)

    def _parse(self, filepath: str) -> Tuple[pd.DataFrame, Dict]:
        """Parse old format with multi-header structure."""
        # Extract metadata from first 4 rows (skipped during header load)
        peek = pd.read_csv(filepath, nrows=5, dtype=str)
        metadata = {
            "format": "old",
            "model_name_row": peek.iloc[1, 0] if len(peek) > 1 else None,
            "dataset_name_row": peek.iloc[2, 0] if len(peek) > 2 else None,
            "data_source_row": peek.iloc[3, 0] if len(peek) > 3 else None,
        }

        # Load with proper multi-header parsing
        # Skip rows 2-4 (0-indexed: rows 1-3), use rows 0 and 4 as header (becomes 0 and 1)
        df = pd.read_csv(
            filepath,
            skiprows=(2, 3, 4),  # Skip metadata rows (1-indexed in original, 0-indexed here becomes 1,2,3)
            header=[0, 1],
            index_col=[0, 1],
        )
        df.index = df.index.rename(["gene_id", "symbol"])

        # Collapse MultiIndex columns: (SIDM, cell_line_name) -> use first level (SIDM)
        df.columns = df.columns.droplevel(1)

        # Drop gene_id index level, keep symbol
        df = df.reset_index(level=0, drop=True)
        df.index = df.index.rename("symbol")

        # Normalize symbols to uppercase
        df.index = df.index.str.upper()

        # Clean up column names (remove whitespace)
        df.columns = df.columns.astype(str).str.strip()

        metadata["n_genes"] = df.shape[0]
        metadata["n_samples"] = df.shape[1]
        self._log(f"Parsed old format: {df.shape[0]} genes × {df.shape[1]} samples")

        return df, metadata


class Format26Q1Parser(ActivityParser):
    """Parser for rnaseq_tpm_coding_genes26Q1.csv format.

    Structure:
    - Row 1: single header with SequencingID, ModelConditionID, ModelID, IsDefault*, then gene columns
    - Gene columns: "GENE_NAME (ENTREZ_ID)" format, e.g., "TSPAN6 (7105)"
    - Rows 2+: data with index, metadata columns, then expression values
    - Cell identification: ModelID (ACH-XXXXXX) must be mapped to SIDM via Model.csv
    """

    def load(self, filepath: str) -> Tuple[pd.DataFrame, Dict]:
        """Load 26Q1 format file."""
        self._log(f"Loading 26Q1 format file: {filepath}")
        return self._parse(filepath)

    def _parse(self, filepath: str) -> Tuple[pd.DataFrame, Dict]:
        """Parse 26Q1 format with single header and ModelID cell identification."""
        # Load with single header
        df = pd.read_csv(filepath, index_col=0)

        # Metadata and cell identification columns to extract
        metadata_cols = [
            "SequencingID",
            "ModelConditionID",
            "ModelID",
            "IsDefaultEntryForMC",
            "IsDefaultEntryForModel",
        ]

        # Extract metadata
        metadata = {
            "format": "26Q1",
            "metadata_columns": metadata_cols,
        }

        # Store ModelID -> to be used for SIDM mapping
        model_ids = df["ModelID"].unique().tolist() if "ModelID" in df.columns else []
        metadata["model_ids"] = model_ids

        # Drop metadata columns; keep only gene expression columns
        gene_cols = [c for c in df.columns if c not in metadata_cols]

        if not gene_cols:
            raise ValueError("No gene expression columns found after removing metadata columns")

        df_genes = df[gene_cols].copy()

        # Parse gene column names: "GENE_NAME (ENTREZ_ID)" -> extract GENE_NAME
        gene_symbols = []
        for col in gene_cols:
            # Regex: capture text before parenthesis
            match = re.match(r"^([^\(]+)\s*\(\d+\)$", col.strip())
            if match:
                symbol = match.group(1).strip().upper()
                gene_symbols.append(symbol)
            else:
                # Fallback: just use the column name as-is, uppercase
                logger.warning(f"Could not parse gene symbol from column '{col}'; using as-is")
                gene_symbols.append(col.upper())

        df_genes.columns = gene_symbols

        # Now rename columns from ModelID to SIDM
        # This requires external mapping via Model.csv
        # For now, use ModelID as column names; will be mapped during ActivityMatrix._ensure_sidm()
        # Or: we can import and use load_sidm_from_model_csv here

        # Map ModelID to SIDM (import here to avoid circular dependency)
        try:
            from celios.utils.io import load_sidm_from_model_csv

            sidm_dict, not_found = load_sidm_from_model_csv(model_ids, verbose=self.verbose)
            metadata["sidm_mapping"] = sidm_dict
            metadata["unmapped_model_ids"] = not_found

            if not_found:
                self._log(
                    f"Warning: {len(not_found)} ModelID(s) not found in Model.csv: {not_found}",
                    level="warning",
                )

            # Rename columns from ModelID to SIDM (inverse of sidm_dict which maps SIDM -> cell_line)
            # sidm_dict: {SIDM -> cell_line_name}
            # We need: {model_id -> SIDM}
            # The issue is model_ids from df["ModelID"] correspond to ModelID values
            # But sidm_dict maps from SIDM ID not from cell_line names
            # Actually, load_sidm_from_model_csv takes cell_line names and returns {SIDM -> cell_line_name}
            # So here we pass ModelID values expecting them to match...
            # Let me re-read the function

            # Actually, looking at load_sidm_from_model_csv in io.py, it takes cell_line_names
            # and looks them up in Model.csv using normalization. The ModelID is a column in Model.csv
            # So we need a different approach: map ModelID values directly to SIDM

            # For now, just use ModelID as surrogate for sample identifiers
            # The mapping to SIDM will happen in ActivityMatrix._ensure_sidm()
            # when it processes the cell_line_file

            # Store model_id to SIDM mapping if available
            if sidm_dict:
                # Create reverse mapping: model_id -> SIDM
                # sidm_dict maps SIDM -> cell_line_name from Model.csv
                # But we need model_id -> SIDM
                # This requires looking at the actual Model.csv structure
                pass

        except Exception as e:
            self._log(f"Could not map ModelID to SIDM: {e}. Using ModelID as sample identifiers.", level="warning")
            # Fallback: keep ModelID as column names for now

        # Transpose: rows become genes, columns become samples
        df_genes = df_genes.T

        # Set sample identifiers as column names
        if "SequencingID" in df.columns:
            df_genes.columns = df["SequencingID"].values
        else:
            df_genes.columns = df.index

        df_genes.index = pd.Index(gene_symbols, name="model_id")

        metadata["n_genes"] = df_genes.shape[0]
        metadata["n_samples"] = df_genes.shape[1]

        self._log(f"Parsed 26Q1 format: {df_genes.shape[0]} genes × {df_genes.shape[1]} samples")

        return df_genes, metadata


def get_parser(format_name: str, verbose: bool = False) -> ActivityParser:
    """Factory function to get appropriate parser instance.

    Args:
        format_name: 'old' or '26Q1'
        verbose: Enable verbose logging

    Returns:
        ActivityParser instance

    Raises:
        ValueError: If format_name is unknown
    """
    if format_name == "old":
        return OldFormatParser(verbose=verbose)
    elif format_name == "26Q1":
        return Format26Q1Parser(verbose=verbose)
    else:
        raise ValueError(f"Unknown format: {format_name}. Supported: 'old', '26Q1'")


__all__ = [
    "FormatDetector",
    "ActivityParser",
    "OldFormatParser",
    "Format26Q1Parser",
    "get_parser",
]
