"""
test_activity_parser.py

Unit tests for activity file format detection and parsing.

Tests both old (rnaseq_tpm_20220624) and 26Q1 (rnaseq_tpm_coding_genes26Q1) formats.
"""

import pytest
import pandas as pd
import os
import tempfile
from pathlib import Path

from celios.features.activity_parser import (
    FormatDetector,
    OldFormatParser,
    Format26Q1Parser,
    get_parser,
)


class TestFormatDetector:
    """Tests for format auto-detection."""

    def test_detect_old_format_real_file(self):
        """Test detection of old format on real file."""
        old_file = Path(__file__).parent.parent.parent.parent / "data" / "activity_input" / "rnaseq_tpm_20220624.csv"
        if not old_file.exists():
            pytest.skip(f"Test file not found: {old_file}")
        
        detected = FormatDetector.detect(str(old_file))
        assert detected == "old", f"Expected 'old' format, got '{detected}'"

    def test_detect_26q1_format_real_file(self):
        """Test detection of 26Q1 format on real file."""
        q1_file = Path(__file__).parent.parent.parent.parent / "data" / "activity_input" / "rnaseq_tpm_coding_genes26Q1.csv"
        if not q1_file.exists():
            pytest.skip(f"Test file not found: {q1_file}")
        
        detected = FormatDetector.detect(str(q1_file))
        assert detected == "26Q1", f"Expected '26Q1' format, got '{detected}'"

    def test_format_override_valid(self):
        """Test format_override parameter with valid value."""
        detected = FormatDetector.detect("dummy.csv", format_override="old")
        assert detected == "old"

    def test_format_override_invalid(self):
        """Test format_override raises error for invalid format."""
        with pytest.raises(ValueError, match="Unknown format override"):
            FormatDetector.detect("dummy.csv", format_override="invalid_format")


class TestOldFormatParser:
    """Tests for old format parser."""

    def test_parse_old_format_real_file(self):
        """Test parsing of old format file."""
        old_file = Path(__file__).parent.parent.parent.parent / "data" / "activity_input" / "rnaseq_tpm_20220624.csv"
        if not old_file.exists():
            pytest.skip(f"Test file not found: {old_file}")
        
        parser = OldFormatParser(verbose=True)
        df, metadata = parser.load(str(old_file))
        
        # Check output structure
        assert isinstance(df, pd.DataFrame), "Parser should return DataFrame"
        assert isinstance(metadata, dict), "Parser should return metadata dict"
        
        # Check DataFrame properties
        assert df.index.name == "symbol", "Index should be named 'symbol'"
        assert len(df) > 0, "DataFrame should have rows (genes)"
        assert len(df.columns) > 0, "DataFrame should have columns (samples)"
        
        # Check metadata
        assert metadata["format"] == "old"
        assert "n_genes" in metadata
        assert "n_samples" in metadata
        assert metadata["n_genes"] == df.shape[0]
        assert metadata["n_samples"] == df.shape[1]
        
        # Check gene symbols are uppercase
        assert all(s == s.upper() for s in df.index.str.upper()), "Gene symbols should be normalized to uppercase"

    def test_parse_old_format_column_types(self):
        """Test that old format parser returns numeric data."""
        old_file = Path(__file__).parent.parent.parent.parent / "data" / "activity_input" / "rnaseq_tpm_20220624.csv"
        if not old_file.exists():
            pytest.skip(f"Test file not found: {old_file}")
        
        parser = OldFormatParser()
        df, _ = parser.load(str(old_file))
        
        # Check that all values are numeric
        assert df.select_dtypes(include=['number']).shape[1] > 0, "Should have numeric columns"


class TestFormat26Q1Parser:
    """Tests for 26Q1 format parser."""

    def test_parse_26q1_format_real_file(self):
        """Test parsing of 26Q1 format file."""
        q1_file = Path(__file__).parent.parent.parent.parent / "data" / "activity_input" / "rnaseq_tpm_coding_genes26Q1.csv"
        if not q1_file.exists():
            pytest.skip(f"Test file not found: {q1_file}")
        
        parser = Format26Q1Parser(verbose=True)
        df, metadata = parser.load(str(q1_file))
        
        # Check output structure
        assert isinstance(df, pd.DataFrame), "Parser should return DataFrame"
        assert isinstance(metadata, dict), "Parser should return metadata dict"
        
        # Check DataFrame properties
        assert df.index.name == "symbol", "Index should be named 'symbol' (gene symbols)"
        assert len(df) > 0, "DataFrame should have rows (genes)"
        assert len(df.columns) > 0, "DataFrame should have columns (samples/ModelIDs)"
        
        # Check metadata
        assert metadata["format"] == "26Q1"
        assert "n_genes" in metadata
        assert "n_samples" in metadata
        assert metadata["n_genes"] == df.shape[0]
        assert metadata["n_samples"] == df.shape[1]

    def test_parse_26q1_gene_symbol_extraction(self):
        """Test that gene symbols are extracted correctly from 'NAME (ID)' format."""
        q1_file = Path(__file__).parent.parent.parent.parent / "data" / "activity_input" / "rnaseq_tpm_coding_genes26Q1.csv"
        if not q1_file.exists():
            pytest.skip(f"Test file not found: {q1_file}")
        
        parser = Format26Q1Parser()
        df, _ = parser.load(str(q1_file))
        
        # Check that gene symbols don't contain parentheses
        for symbol in df.index:
            assert "(" not in str(symbol), f"Gene symbol '{symbol}' should not contain '('"
            assert ")" not in str(symbol), f"Gene symbol '{symbol}' should not contain ')'"
        
        # Check that symbols are uppercase
        assert all(s == s.upper() for s in df.index), "Gene symbols should be normalized to uppercase"

    def test_parse_26q1_synthetic_modelid_mapping(self):
        """Test 26Q1 parsing with synthetic data: ModelID should be used as column names."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a synthetic 26Q1 CSV with ModelID and SequencingID
            synthetic_file = os.path.join(tmpdir, "synthetic_26q1.csv")
            
            # Create synthetic data
            synthetic_data = """SequencingID,ModelConditionID,ModelID,IsDefaultEntryForMC,IsDefaultEntryForModel,BRCA1 (672),TP53 (7157),MYC (4609)
SQ00000001,MCC0001,ACH-000111,true,true,10.5,5.2,3.1
SQ00000002,MCC0002,ACH-000222,true,true,8.2,4.9,2.8
SQ00000003,MCC0003,ACH-000333,true,true,9.1,6.3,3.5"""
            
            with open(synthetic_file, 'w') as f:
                f.write(synthetic_data)
            
            # Parse with Format26Q1Parser
            parser = Format26Q1Parser()
            df, metadata = parser.load(synthetic_file)
            
            # Verify output
            assert df.shape == (3, 3), f"Expected shape (3, 3), got {df.shape}"
            assert df.index.name == "symbol", "Index should be named 'symbol'"
            assert list(df.index) == ["BRCA1", "TP53", "MYC"], "Gene symbols should be extracted"
            
            # CRITICAL: Columns should be ModelID values, NOT SequencingID
            expected_modelids = ["ACH-000111", "ACH-000222", "ACH-000333"]
            assert list(df.columns) == expected_modelids, \
                f"Columns should be ModelID values {expected_modelids}, got {list(df.columns)}"
            
            # Verify metadata has model_ids
            assert "model_ids" in metadata, "Metadata should contain model_ids"
            assert set(metadata["model_ids"]) == set(expected_modelids), \
                f"Metadata model_ids should contain {expected_modelids}"
            
            # Verify numeric values
            assert df.loc["BRCA1", "ACH-000111"] == 10.5, "Expression values should be parsed correctly"
            assert df.loc["TP53", "ACH-000222"] == 4.9, "Expression values should be parsed correctly"


class TestGetParser:
    """Tests for factory function."""

    def test_get_old_parser(self):
        """Test that get_parser returns OldFormatParser for 'old' format."""
        parser = get_parser("old")
        assert isinstance(parser, OldFormatParser)

    def test_get_26q1_parser(self):
        """Test that get_parser returns Format26Q1Parser for '26Q1' format."""
        parser = get_parser("26Q1")
        assert isinstance(parser, Format26Q1Parser)

    def test_get_parser_invalid_format(self):
        """Test that get_parser raises error for invalid format."""
        with pytest.raises(ValueError, match="Unknown format"):
            get_parser("invalid_format")

    def test_get_parser_verbose(self):
        """Test that get_parser passes verbose flag."""
        parser = get_parser("old", verbose=True)
        assert parser.verbose is True

        parser_quiet = get_parser("old", verbose=False)
        assert parser_quiet.verbose is False


class TestIntegration:
    """Integration tests using real files."""

    def test_old_format_end_to_end(self):
        """Test complete workflow with old format file."""
        old_file = Path(__file__).parent.parent.parent.parent / "data" / "activity_input" / "rnaseq_tpm_20220624.csv"
        if not old_file.exists():
            pytest.skip(f"Test file not found: {old_file}")
        
        # Detect format
        detected_format = FormatDetector.detect(str(old_file))
        assert detected_format == "old"
        
        # Get parser
        parser = get_parser(detected_format)
        
        # Load file
        df, metadata = parser.load(str(old_file))
        
        # Verify output contract
        assert isinstance(df, pd.DataFrame)
        assert isinstance(metadata, dict)
        assert metadata["format"] == "old"
        assert df.shape[0] > 0 and df.shape[1] > 0

    def test_26q1_format_end_to_end(self):
        """Test complete workflow with 26Q1 format file."""
        q1_file = Path(__file__).parent.parent.parent.parent / "data" / "activity_input" / "rnaseq_tpm_coding_genes26Q1.csv"
        if not q1_file.exists():
            pytest.skip(f"Test file not found: {q1_file}")
        
        # Detect format
        detected_format = FormatDetector.detect(str(q1_file))
        assert detected_format == "26Q1"
        
        # Get parser
        parser = get_parser(detected_format)
        
        # Load file
        df, metadata = parser.load(str(q1_file))
        
        # Verify output contract
        assert isinstance(df, pd.DataFrame)
        assert isinstance(metadata, dict)
        assert metadata["format"] == "26Q1"
        assert df.shape[0] > 0 and df.shape[1] > 0
        
        # Verify that columns are ModelID values (for downstream SIDM mapping)
        # This is critical: columns should be compatible with alias_map lookup
        assert "model_ids" in metadata, "Metadata should contain model_ids"
        assert len(metadata["model_ids"]) == df.shape[1], "model_ids count should match column count"

    def test_26q1_expression_builder_integration(self):
        """Test that 26Q1 parser output integrates with expression_builder for SIDM mapping."""
        from celios.features.expression_builder import prepare_expression_matrix
        
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a synthetic 26Q1 CSV
            synthetic_file = os.path.join(tmpdir, "synthetic_26q1.csv")
            
            synthetic_data = """SequencingID,ModelConditionID,ModelID,IsDefaultEntryForMC,IsDefaultEntryForModel,BRCA1 (672),TP53 (7157)
SQ00000001,MCC0001,ACH-000111,true,true,10.5,5.2
SQ00000002,MCC0002,ACH-000222,true,true,8.2,4.9"""
            
            with open(synthetic_file, 'w') as f:
                f.write(synthetic_data)
            
            # Create alias_map: ModelID → SIDM
            alias_map = {
                "ACH-000111": "SIDM00001",
                "ACH-000222": "SIDM00002",
            }
            sidm_list = ["SIDM00001", "SIDM00002"]
            
            # Prepare expression matrix with SIDM mapping
            df, metadata = prepare_expression_matrix(
                activity_file=synthetic_file,
                format_override="26Q1",
                alias_map=alias_map,
                sidm_list=sidm_list,
                verbose=True,
                deep_debug=False,
            )
            
            # Verify output
            assert df.shape[0] == 2, f"Expected 2 genes, got {df.shape[0]}"
            assert df.shape[1] == 2, f"Expected 2 SIDMs, got {df.shape[1]}"
            assert df.index.name == "symbol", "Index should be 'symbol'"
            
            # Verify columns are now SIDM (not SequencingID or ModelID)
            expected_sidms = {"SIDM00001", "SIDM00002"}
            actual_sidms = set(df.columns)
            assert actual_sidms == expected_sidms, \
                f"Expected columns to be {expected_sidms}, got {actual_sidms}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
