import pytest
import os
import sys
import pytest

# Make repository root importable
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from celios.core import run_celios

def test_mutations_only_no_activity_file():
    """Test that CELIOS can run without activity file when expression is not needed."""
    
    # Configuration for mutations/CNV only mode (no activity_file, no expression)
    config_mutations_only = {
        "paths": {
            "base": ".",
            "input": "data",
            "output": "results",
            "cellfiles_dir": "results/cell_lines_mutations_only",
        },
        "steps": {
            "Node": {
                "node_input": "node_dic_input/cell_fate_plus.sif",
                "hgnc_symbols_file": "node_dic_input/hgnc_complete_set.txt",
                "manual_symbols_file": "node_dic_input/manual_symbols.csv",
                "include_alias_previous_symbols": False,
                "directory_output": "results",
            },
            "Activity": {
                # NOTE: No activity_file specified - using only binary matrices
                "cell_line_file": "vis_2024/cell_line_list.csv",
                "mutations_file": "activity_input/CCLE_muts_binary.csv",
                "cnv_file": "activity_input/CCLE_CNV_binary.csv",
                "directory_output": "results",
                "data_sources": ["mutations", "cnv"],  # No 'expression' - activity_file not needed
            },
        },
    }
    
    # Check required data files exist
    required_files = {
        "cell_line_file": os.path.join(ROOT, "data", "vis_2024", "cell_line_list.csv"),
        "mutations_file": os.path.join(ROOT, "data", "activity_input", "CCLE_muts_binary.csv"),
        "cnv_file": os.path.join(ROOT, "data", "activity_input", "CCLE_CNV_binary.csv"),
        "node_input": os.path.join(ROOT, "data", "node_dic_input", "cell_fate_plus.sif"),
    }
    
    missing = [name for name, path in required_files.items() if not os.path.exists(path)]
    if missing:
        pytest.skip(f"Required data files not found: {', '.join(missing)}")
    
    # Run pipeline without activity_file
    artifacts = run_celios(config=config_mutations_only, verbose=True)
    
    # Verify results
    assert artifacts is not None, "Pipeline should return artifacts"
    assert isinstance(artifacts, dict), "Artifacts should be dict"
    
    # Check the activity matrix
    matrix = artifacts.get("activity_matrix")
    print(f"\nActivity matrix type: {type(matrix)}")
    if matrix is not None:
        print(f"Activity matrix shape: {matrix.shape}")
        print(f"Activity matrix columns: {list(matrix.columns)[:10]}")  # Show first 10 cols
        assert len(matrix) > 0, "Activity matrix should not be empty"
    else:
        print("Activity matrix is None (no results for matched cell lines)")

if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
