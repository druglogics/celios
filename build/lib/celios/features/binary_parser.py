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

        # Check for format indicators in header
        if "gene_symbol" in header:
            logger.info(f"Auto-detected OLD binary format from {Path(filepath).name}")
            return "old"
        elif "ModelID" in header or "model_id" in header.lower():
            logger.info(f"Auto-detected 26Q1 binary format from {Path(filepath).name}")
            return "26q1"
        else:
            raise ValueError(
                f"Cannot auto-detect binary matrix format from {filepath}. "
                "Expected 'gene_symbol' (old) or 'ModelID' (26Q1) in header."
            )


class BinaryMatrixParser(ABC):
    """Abstract base class for binary matrix parsers."""

    def __init__(self, verbose: bool = False):
        """Initialize parser."""
        self.verbose = verbose

    def load(
        self, filepath: str, model_registry: Optional[str] = None
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
        df, metadata = self._parse(filepath, model_registry)
        return df, metadata

    @abstractmethod
    def _parse(
        self, filepath: str, model_registry: Optional[str] = None
    ) -> Tuple[pd.DataFrame, Dict]:
        """Parse file and return normalized DataFrame. Implemented by subclasses."""
        pass


class OldBinaryMatrixParser(BinaryMatrixParser):
    """Parser for OLD binary format (genes × SIDM)."""

    def _parse(
        self, filepath: str, model_registry: Optional[str] = None
    ) -> Tuple[pd.DataFrame, Dict]:
        """
        Parse OLD format: gene_symbol index, SIDM columns.

        Returns:
            DataFrame with genes as index, SIDM as columns, metadata dict
        """
        logger.info(f"Parsing OLD binary format from {Path(filepath).name}")

        # Read with gene_symbol as index (column 1, skip the unnamed first column)
        df = load_csv_file(filepath, index_col=1)  # Column 1 is 'gene_symbol'
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
        self, filepath: str, model_registry: Optional[str] = None
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

        # Column 0 becomes the index (ModelID)
        model_ids = df.index.tolist()

        if self.verbose:
            logger.info(f"Loaded {len(model_ids)} ModelIDs from file")

        # Map ModelID → SIDM
        sidm_mapping, not_found = load_sidm_from_modelid(
            model_ids, model_registry=model_registry, verbose=self.verbose
        )

        if not_found:
            logger.warning(f"ModelIDs not found in registry: {not_found}")

        # Map index from ModelID to SIDM
        df.index = df.index.map(sidm_mapping)
        df.index.name = "SIDM"

        # Normalize gene symbols to uppercase (column names are genes)
        df.columns = df.columns.str.upper()

        # Transpose to genes × SIDM
        df = df.T
        df.index.name = "gene_symbol"

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

        n_genes, n_samples = df.shape
        metadata = {
            "format": "26q1",
            "shape": (n_genes, n_samples),
            "n_genes": n_genes,
            "n_samples": n_samples,
            "model_ids_mapped": len(sidm_mapping),
            "model_ids_not_found": len(not_found),
        }

        if self.verbose:
            logger.info(
                f"26Q1 binary matrix: {n_genes} genes × {n_samples} SIDM samples "
                f"(mapped {len(sidm_mapping)} ModelIDs)"
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
