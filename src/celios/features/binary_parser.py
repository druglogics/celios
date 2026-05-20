"""
Format-agnostic binary matrix parsing (mutations, CNV) with auto-detection.

Supports two formats:
1. OLD (CCLE_*_binary.csv): genes as rows, SIDM as columns
2. 26Q1 (ModelID-based): ModelID as rows, genes as columns

Uses Strategy Pattern for pluggable parsers with auto-detection and override support.
"""

import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, Optional, Tuple

import pandas as pd

from celios.utils.io import load_csv_file, load_sidm_from_modelid

logger = logging.getLogger(__name__)

OLD_GENE_HEADER_NAMES = {
    "gene_symbol",
    "symbol",
    "gene",
    "hugo_symbol",
    "hugosymbol",
}


def _normalize_header_token(token: str) -> str:
    return str(token).strip().strip('"').strip("'").lower()


def _detect_old_gene_column(columns) -> str:
    normalized_columns = [_normalize_header_token(column) for column in columns]

    if normalized_columns and normalized_columns[0] in OLD_GENE_HEADER_NAMES:
        return columns[0]

    if len(columns) > 1 and normalized_columns[0].startswith("unnamed") and normalized_columns[1] in OLD_GENE_HEADER_NAMES:
        return columns[1]

    raise ValueError(
        "Cannot detect old binary matrix gene identifier column. "
        f"Expected one of {sorted(OLD_GENE_HEADER_NAMES)} in the first column."
    )


class FormatDetector:
    """Detect binary matrix format by inspecting file headers."""

    @staticmethod
    def detect(filepath: str, format_override: Optional[str] = None) -> str:
        """
        Detect binary matrix format.

        Args:
            filepath: Path to binary matrix file
            format_override: Force a format ("old" | "26q1" | None for auto-detect)

        Returns:
            Format name: "old" or "26q1"

        Raises:
            ValueError: If format cannot be determined or is invalid
        """
        if format_override:
            if format_override.lower() not in ("old", "26q1"):
                raise ValueError(
                    f"Invalid format_override: {format_override}. "
                    "Must be 'old' or '26q1'."
                )
            logger.info(f"Binary matrix format override: {format_override}")
            return format_override.lower()

        # Read first line to inspect headers
        try:
            with open(filepath, "r") as f:
                header = f.readline().strip()
        except Exception as e:
            raise ValueError(f"Cannot read file {filepath}: {e}")

        # Check for format indicators in the header tokens
        header_tokens = header.split(",")
        normalized_tokens = [_normalize_header_token(token) for token in header_tokens]

        if normalized_tokens and normalized_tokens[0] in OLD_GENE_HEADER_NAMES:
            logger.info(f"Auto-detected OLD binary format from {Path(filepath).name}")
            return "old"
        if len(normalized_tokens) > 1 and normalized_tokens[0] in {"", "unnamed: 0", "unnamed:0"} and normalized_tokens[1] in OLD_GENE_HEADER_NAMES:
            logger.info(f"Auto-detected OLD binary format from {Path(filepath).name}")
            return "old"
        if normalized_tokens and normalized_tokens[0] == "modelid":
            logger.info(f"Auto-detected 26Q1 binary format from {Path(filepath).name}")
            return "26q1"

        raise ValueError(
            f"Cannot auto-detect binary matrix format from {filepath}. "
            "Expected an old-style gene column (e.g. 'symbol') or 'ModelID' in header."
        )


class BinaryMatrixParser(ABC):
    """Abstract base class for binary matrix parsers."""

    def __init__(self, verbose: bool = False):
        """Initialize parser."""
        self.verbose = verbose

    def load(
        self,
        filepath: str,
        model_registry: Optional[str] = None,
        resolve_model_ids: bool = True,
    ) -> Tuple[pd.DataFrame, Dict]:
        """
        Load and parse binary matrix file.

        Args:
            filepath: Path to binary matrix file
            model_registry: Path to Model.csv registry (for ModelID→SIDM mapping)

        Returns:
            Tuple of:
            - DataFrame: genes as index, SIDM as columns, binary values
            - metadata_dict: {'format': format_name, 'shape': (n_genes, n_samples), ...}
        """
        df, metadata = self._parse(filepath, model_registry, resolve_model_ids=resolve_model_ids)
        return df, metadata

    @abstractmethod
    def _parse(
        self,
        filepath: str,
        model_registry: Optional[str] = None,
        resolve_model_ids: bool = True,
    ) -> Tuple[pd.DataFrame, Dict]:
        """Parse file and return normalized DataFrame. Implemented by subclasses."""
        pass


class OldBinaryMatrixParser(BinaryMatrixParser):
    """Parser for OLD binary format (genes × SIDM)."""

    def _parse(
        self,
        filepath: str,
        model_registry: Optional[str] = None,
        resolve_model_ids: bool = True,
    ) -> Tuple[pd.DataFrame, Dict]:
        """
        Parse OLD format: gene_symbol index, SIDM columns.

        Returns:
            DataFrame with genes as index, SIDM as columns, metadata dict
        """
        logger.info(f"Parsing OLD binary format from {Path(filepath).name}")

        df = load_csv_file(filepath)

        gene_col = _detect_old_gene_column(df.columns)
        df = df.set_index(gene_col)
        df.index.name = "gene_symbol"

        # Drop the unnamed index column if present
        df = df.loc[:, ~df.columns.str.contains("^Unnamed")]

        # Normalize gene symbols to uppercase
        df.index = df.index.str.upper()

        # Coerce values to numeric (preserve NaN for missing entries).
        # Empty cells in binary matrices should NOT be interpreted as inactive (0).
        # Keep NaN so downstream training can choose whether to integrate or ignore.
        df = df.apply(pd.to_numeric, errors="coerce")

        # Validate values: allow only 0.0, 1.0, or NaN. Log a warning if others
        unique_vals = pd.unique(df.values.ravel())
        unexpected = [v for v in unique_vals if not (pd.isna(v) or v in (0, 0.0, 1, 1.0))]
        if unexpected:
            logger.warning(
                f"Unexpected non-binary values in OLD format matrix: {unexpected}. "
                "These will be coerced to NaN."
            )
            # Coerce unexpected values to NaN
            df = df.where(df.isin([0, 0.0, 1, 1.0]))

        n_genes, n_samples = df.shape
        metadata = {
            "format": "old",
            "shape": (n_genes, n_samples),
            "input_shape": (n_genes, n_samples),
            "sample_axis": "columns",
            "raw_sample_ids": list(df.columns),
            "n_genes": n_genes,
            "n_samples": n_samples,
        }

        if self.verbose:
            logger.info(
                f"OLD binary matrix: {n_genes} genes × {n_samples} SIDM samples"
            )

        return df, metadata


class ModelIDBinaryMatrixParser(BinaryMatrixParser):
    """Parser for 26Q1 binary format (ModelID × genes, transposed to genes × SIDM)."""

    def _parse(
        self,
        filepath: str,
        model_registry: Optional[str] = None,
        resolve_model_ids: bool = True,
    ) -> Tuple[pd.DataFrame, Dict]:
        """
        Parse 26Q1 format: ModelID index, gene columns.

        Returns:
            DataFrame with genes as index, SIDM as columns (transposed + mapped)
            metadata dict
        """
        logger.info(f"Parsing 26Q1 binary format from {Path(filepath).name}")

        # Read with ModelID as index
        df = pd.read_csv(filepath, index_col=0)
        input_shape = df.shape

        # Column 0 becomes the index (ModelID)
        model_ids = df.index.tolist()

        if self.verbose:
            logger.info(f"Loaded {len(model_ids)} ModelIDs from file")

        model_to_sidm = {}
        not_found = []

        # Normalize gene symbols to uppercase (column names are genes)
        df.columns = df.columns.str.upper()

        # Coerce values to numeric (preserve NaN for missing entries).
        # 26Q1 files may use empty cells for missing data which should NOT be
        # interpreted as inactive (0). Keep NaN so downstream training can
        # choose whether to integrate or ignore missing measurements.
        df = df.apply(pd.to_numeric, errors="coerce")

        # Validate values: allow only 0.0, 1.0, or NaN. Log a warning if others
        unique_vals = pd.unique(df.values.ravel())
        unexpected = [v for v in unique_vals if not (pd.isna(v) or v in (0, 0.0, 1, 1.0))]
        if unexpected:
            logger.warning(
                f"Unexpected non-binary values in 26Q1 matrix: {unexpected}. "
                "These will be coerced to NaN."
            )
            # Coerce unexpected values to NaN
            df = df.where(df.isin([0, 0.0, 1, 1.0]))

        if resolve_model_ids:
            # Backward-compatible path: map ModelID → SIDM during parsing.
            model_to_sidm, not_found = load_sidm_from_modelid(
                model_ids, model_registry=model_registry, verbose=self.verbose
            )

            if not_found:
                logger.warning(f"ModelIDs not found in registry: {not_found}")

            # Map index from ModelID to SIDM and drop any rows that could not be resolved.
            df.index = df.index.map(model_to_sidm)
            df = df[df.index.notna()]
            df.index.name = "SIDM"

            # Transpose to genes × SIDM
            df = df.T
            df.index.name = "gene_symbol"

            n_genes, n_samples = df.shape
            metadata = {
                "format": "26q1",
                "shape": (n_genes, n_samples),
                "input_shape": input_shape,
                "sample_axis": "columns",
                "raw_sample_ids": model_ids,
                "n_genes": n_genes,
                "n_samples": n_samples,
                "model_ids_mapped": len(model_to_sidm),
                "model_ids_not_found": len(not_found),
                "resolve_model_ids": resolve_model_ids,
            }

            if self.verbose:
                logger.info(
                    f"26Q1 binary matrix: {n_genes} genes × {n_samples} SIDM samples "
                    f"(mapped {len(model_to_sidm)} ModelIDs)"
                )

            return df, metadata

        # Keep raw sample identifiers for local normalization in the builder.
        df.index.name = "ModelID"

        n_samples, n_genes = df.shape
        metadata = {
            "format": "26q1",
            "shape": (n_samples, n_genes),
            "input_shape": input_shape,
            "sample_axis": "rows",
            "raw_sample_ids": model_ids,
            "n_genes": n_genes,
            "n_samples": n_samples,
            "model_ids_mapped": 0,
            "model_ids_not_found": 0,
            "resolve_model_ids": resolve_model_ids,
        }

        if self.verbose:
            logger.info(
                f"26Q1 binary matrix: {n_samples} ModelID samples × {n_genes} genes "
                "(local SIDM resolution disabled)"
            )

        return df, metadata


def get_binary_parser(format_name: str, verbose: bool = False) -> BinaryMatrixParser:
    """
    Factory function to get appropriate binary matrix parser.

    Args:
        format_name: "old" or "26q1"
        verbose: Enable verbose logging

    Returns:
        Parser instance

    Raises:
        ValueError: If format_name is invalid
    """
    format_name = format_name.lower()

    if format_name == "old":
        return OldBinaryMatrixParser(verbose=verbose)
    elif format_name == "26q1":
        return ModelIDBinaryMatrixParser(verbose=verbose)
    else:
        raise ValueError(
            f"Unknown binary matrix format: {format_name}. "
            "Must be 'old' or '26q1'."
        )
