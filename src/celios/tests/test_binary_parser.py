"""
Test suite for binary_parser module (mutations/CNV format detection and parsing).

Tests auto-detection, format-specific parsers, and end-to-end integration.
"""

import os
import pytest
import pandas as pd
import numpy as np
from pathlib import Path

from celios.features.binary_parser import (
    FormatDetector,
    OldBinaryMatrixParser,
    ModelIDBinaryMatrixParser,
    get_binary_parser,
)


# ============================================================================
# Test Data Setup
# ============================================================================

# Path to activity input files (used as reference for real file existence)
DATA_DIR = Path(__file__).parent.parent.parent.parent / "data" / "activity_input"
CCLE_MUTS_FILE = DATA_DIR / "CCLE_muts_binary.csv"


@pytest.fixture(scope="session")
def data_files_exist():
    """Check if real data files are available."""
    return {
        "ccle_muts": CCLE_MUTS_FILE.exists(),
    }


# ============================================================================
# Test FormatDetector
# ============================================================================

class TestFormatDetector:
    """Test format auto-detection with real and synthetic data."""

    def test_detect_old_format_real_file(self, data_files_exist):
        """Detect OLD binary format from real CCLE file."""
        if not data_files_exist["ccle_muts"]:
            pytest.skip("CCLE_muts_binary.csv not found")

        detected = FormatDetector.detect(str(CCLE_MUTS_FILE))
        assert detected == "old", f"Expected 'old', got {detected}"

    def test_format_override_valid(self, data_files_exist):
        """Verify format override for valid formats."""
        if not data_files_exist["ccle_muts"]:
            pytest.skip("CCLE_muts_binary.csv not found")

        # Should return override value without inspecting file
        detected = FormatDetector.detect(str(CCLE_MUTS_FILE), format_override="26q1")
        assert detected == "26q1"

        detected = FormatDetector.detect(str(CCLE_MUTS_FILE), format_override="old")
        assert detected == "old"

    def test_format_override_invalid(self, data_files_exist):
        """Reject invalid format_override values."""
        if not data_files_exist["ccle_muts"]:
            pytest.skip("CCLE_muts_binary.csv not found")

        with pytest.raises(ValueError, match="Invalid format_override"):
            FormatDetector.detect(str(CCLE_MUTS_FILE), format_override="invalid")


# ============================================================================
# Test OldBinaryMatrixParser
# ============================================================================

class TestOldBinaryMatrixParser:
    """Test OLD binary format parsing (genes × SIDM)."""

    def test_parse_old_format_real_file(self, data_files_exist):
        """Parse OLD binary format from real CCLE file."""
        if not data_files_exist["ccle_muts"]:
            pytest.skip("CCLE_muts_binary.csv not found")

        parser = OldBinaryMatrixParser(verbose=False)
        df, metadata = parser.load(str(CCLE_MUTS_FILE))

        # Check structure
        assert isinstance(df, pd.DataFrame), "Output should be DataFrame"
        assert df.index.name == "gene_symbol", "Index should be named 'gene_symbol'"
        assert all(df.columns.str.startswith("SIDM")), "Columns should be SIDM identifiers"

        # Check metadata
        assert metadata["format"] == "old"
        assert "shape" in metadata
        n_genes, n_samples = metadata["shape"]
        assert n_genes == len(df), "Metadata shape should match DataFrame"
        assert n_samples == len(df.columns)

        # Check values are binary (0.0, 1.0) or NaN (preserved for missing data)
        unique_vals = pd.unique(df.values.flatten())
        non_nan = [v for v in unique_vals if not pd.isna(v)]
        assert all(v in (0, 0.0, 1, 1.0) for v in non_nan), "All non-NaN values should be 0.0 or 1.0"

        # Check gene symbols are uppercase
        assert all(gene.isupper() for gene in df.index), "Gene symbols should be uppercase"

    def test_parse_old_format_column_types(self, data_files_exist):
        """Verify correct parsing and NaN handling in OLD format."""
        if not data_files_exist["ccle_muts"]:
            pytest.skip("CCLE_muts_binary.csv not found")

        parser = OldBinaryMatrixParser(verbose=False)
        df, _ = parser.load(str(CCLE_MUTS_FILE))

        # Verify all columns are float64 (to accommodate NaN preservation)
        assert df.dtypes.apply(lambda x: x == float).all(), "All columns should be float64 (for NaN support)"
        # Verify NaN values are preserved in the matrix (not filled to 0)
        assert df.isna().sum().sum() > 0, "Matrix should contain NaN values"

    def test_nan_preservation_in_old_parser(self, data_files_exist):
        """Verify NaN values are preserved (not coerced to 0) in OLD parser."""
        if not data_files_exist["ccle_muts"]:
            pytest.skip("CCLE_muts_binary.csv not found")

        parser = OldBinaryMatrixParser(verbose=False)
        df, _ = parser.load(str(CCLE_MUTS_FILE))

        # Check that NaN count is > 0 (empty cells preserved as NaN, not filled to 0)
        nan_count = df.isna().sum().sum()
        assert nan_count > 0, "OLD parser should preserve NaN values from empty cells"


# ============================================================================
# Test ModelIDBinaryMatrixParser (requires Model.csv)
# ============================================================================

class TestModelIDBinaryMatrixParser:
    """Test 26Q1 binary format parsing (ModelID × genes, transposed)."""

    def test_parser_initialization(self):
        """Verify parser can be initialized."""
        parser = ModelIDBinaryMatrixParser(verbose=False)
        assert isinstance(parser, ModelIDBinaryMatrixParser)

    def test_parse_modelid_format_with_partial_registry(self, tmp_path):
        """Parse 26Q1 format and keep only ModelIDs that exist in the registry."""
        activity_file = tmp_path / "binary_26q1.csv"
        activity_file.write_text(
            'ModelID,gene_a,gene_b\n'
            'ACH-000001,1,0\n'
            'ACH-999999,0,1\n',
            encoding="utf-8",
        )

        model_registry = tmp_path / "Model.csv"
        model_registry.write_text(
            'ModelID,SangerModelID\n'
            'ACH-000001,SIDM00001\n',
            encoding="utf-8",
        )

        parser = ModelIDBinaryMatrixParser(verbose=False)
        df, metadata = parser.load(str(activity_file), model_registry=str(model_registry))

        assert df.index.name == "gene_symbol"
        assert list(df.columns) == ["SIDM00001"]
        assert df.loc["GENE_A", "SIDM00001"] == 1
        assert df.loc["GENE_B", "SIDM00001"] == 0
        assert metadata["model_ids_mapped"] == 1
        assert metadata["model_ids_not_found"] == 1


# ============================================================================
# Test get_binary_parser Factory
# ============================================================================

class TestGetBinaryParser:
    """Test factory function for binary matrix parsers."""

    def test_get_old_parser(self):
        """Factory returns OldBinaryMatrixParser for 'old' format."""
        parser = get_binary_parser("old")
        assert isinstance(parser, OldBinaryMatrixParser)

    def test_get_26q1_parser(self):
        """Factory returns ModelIDBinaryMatrixParser for '26q1' format."""
        parser = get_binary_parser("26q1")
        assert isinstance(parser, ModelIDBinaryMatrixParser)

    def test_get_parser_case_insensitive(self):
        """Factory should accept both uppercase and lowercase format names."""
        parser_old = get_binary_parser("OLD")
        assert isinstance(parser_old, OldBinaryMatrixParser)

        parser_26q1 = get_binary_parser("26Q1")
        assert isinstance(parser_26q1, ModelIDBinaryMatrixParser)

    def test_get_parser_invalid_format(self):
        """Factory rejects invalid format names."""
        with pytest.raises(ValueError, match="Unknown binary matrix format"):
            get_binary_parser("invalid_format")

    def test_get_parser_verbose(self):
        """Factory passes verbose flag to parser."""
        parser = get_binary_parser("old", verbose=True)
        assert parser.verbose is True


# ============================================================================
# Test End-to-End Integration
# ============================================================================

class TestIntegration:
    """Integration tests for full parsing pipeline."""

    def test_old_format_end_to_end(self, data_files_exist):
        """Full pipeline: detect -> parse -> validate OLD format."""
        if not data_files_exist["ccle_muts"]:
            pytest.skip("CCLE_muts_binary.csv not found")

        # Auto-detect
        format_name = FormatDetector.detect(str(CCLE_MUTS_FILE))
        assert format_name == "old"

        # Get parser
        parser = get_binary_parser(format_name)

        # Load
        df, metadata = parser.load(str(CCLE_MUTS_FILE))

        # Validate
        assert isinstance(df, pd.DataFrame)
        assert df.shape[0] > 0, "Should have at least one gene"
        assert df.shape[1] > 0, "Should have at least one sample"
        assert metadata["format"] == "old"
        assert metadata["shape"] == (df.shape[0], df.shape[1])

    def test_format_override_end_to_end(self, data_files_exist):
        """Full pipeline with explicit format override."""
        if not data_files_exist["ccle_muts"]:
            pytest.skip("CCLE_muts_binary.csv not found")

        # Force OLD format via override
        format_name = FormatDetector.detect(str(CCLE_MUTS_FILE), format_override="old")
        parser = get_binary_parser(format_name)
        df, metadata = parser.load(str(CCLE_MUTS_FILE))

        assert metadata["format"] == "old"
        assert df.shape[0] > 0


# ============================================================================
# Run Tests
# ============================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
