"""
Test script for CELIOS pipeline.

Tests both legacy mode (cellfiles_dir) and tissue-aware mode (tissue_dir).

Requirements:
- pandas
- All required data files in data/ directory
- For legacy mode: data/activity_input/cell_line_list.csv
- For tissue mode: data/vis_2024/all_cell_lines.csv
"""


import traceback
import os
import sys

# Make repository root importable when running the script directly
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

# Check for required dependencies
try:
    import pandas
except ImportError:
    print("ERROR: pandas is required to run this test. Please install it with: pip install pandas")
    sys.exit(1)

from celios.core import run_celios
config_legacy = {
    "paths": {
        "base": ".",
        "input": "data",
        "output": "results",
        "cellfiles_dir": "results/cell_lines",
    },
    "steps": {
        "Node": {
            # relative to `paths.input`
            "node_input": "node_dic_input/cell_fate_plus.sif",
            "hgnc_symbols_file": "node_dic_input/hgnc_complete_set.txt",
            "manual_symbols_file": "node_dic_input/manual_symbols.csv",
            "include_alias_previous_symbols": False,
            "directory_output": "results",
        },
        "Activity": {
            "activity_file": "activity_input/rnaseq_tpm_20220624.csv",
            "cell_line_file": "activity_input/cell_line_list.csv",
            "tf_activity_file": "activity_input/ccle_tf_activities.csv",
            "mutations_file": "activity_input/CCLE_muts_binary.csv",
            "cnv_file": "activity_input/CCLE_CNV_binary.csv",
            "directory_output": "results",
            "data_sources": ["mutations", "cnv", "TF",],
        },
    },
}

# Tissue-aware mode config (tissue_dir)
# Note: Requires data/vis_2024/all_cell_lines.csv with columns: Tissue, cell_line_name, SIDM
config_tissue = {
    "paths": {
        "base": ".",
        "input": "data",
        "output": "results",
        "tissue_dir": "results/tissue_folders",
    },
    "steps": {
        "Node": {
            # relative to `paths.input`
            "node_input": "node_dic_input/cell_fate_plus.sif",
            "hgnc_symbols_file": "node_dic_input/hgnc_complete_set.txt",
            "manual_symbols_file": "node_dic_input/manual_symbols.csv",
            "include_alias_previous_symbols": False,
            "directory_output": "results",
        },
        "Activity": {
            "activity_file": "activity_input/rnaseq_tpm_20220624.csv",
            "cell_line_file": "vis_2024/all_cell_lines.csv",  # Includes tissue column for tissue-aware mode
            "tf_activity_file": "activity_input/ccle_tf_activities.csv",
            "mutations_file": "activity_input/CCLE_muts_binary.csv",
            "cnv_file": "activity_input/CCLE_CNV_binary.csv",
            "directory_output": "results",
            "data_sources": ["mutations", "cnv", "TF",],
        },
    },
}


def test_celios_config(config, mode_name):
    """Test CELIOS with a given configuration."""
    print(f"\n{'='*50}")
    print(f"Testing {mode_name} mode")
    print(f"{'='*50}")
    
    try:
        # Run the full pipeline (plan=False). This will perform Node and
        # Activity steps and write outputs to the configured `paths.output`.
        artifacts = run_celios(config=config, plan=False, verbose=True)
        print(f'Success: {mode_name} mode returned artifacts of type', type(artifacts))
        try:
            # pretty-print keys and small summaries
            if isinstance(artifacts, dict):
                for k, v in artifacts.items():
                    if k == 'activity_matrix' and v is not None:
                        print(f"- {k}: DataFrame with shape {getattr(v, 'shape', 'unknown')}")
                    else:
                        print(f"- {k}: {type(v)}")
        except Exception:
            pass
        return True
    except Exception as e:
        print(f'Error running {mode_name} mode:')
        traceback.print_exc()
        return False


if __name__ == '__main__':
    # Test both legacy and tissue-aware modes
    results = []
    
    # Test legacy mode
    results.append(test_celios_config(config_legacy, "Legacy"))
    
    # Test tissue-aware mode (only if tissue cell line file exists)
    tissue_cell_file = os.path.join(ROOT, "data", "vis_2024", "all_cell_lines.csv")
    if os.path.exists(tissue_cell_file):
        results.append(test_celios_config(config_tissue, "Tissue-aware"))
    else:
        print(f"\n{'='*50}")
        print("Skipping Tissue-aware mode test")
        print(f"{'='*50}")
        print(f"Tissue cell line file not found: {tissue_cell_file}")
        print("Create this file with columns: Tissue, cell_line_name, SIDM")
        results.append(None)
    
    # Summary
    print(f"\n{'='*50}")
    print("TEST SUMMARY")
    print(f"{'='*50}")
    print(f"Legacy mode: {'PASSED' if results[0] else 'FAILED'}")
    if results[1] is not None:
        print(f"Tissue-aware mode: {'PASSED' if results[1] else 'FAILED'}")
    else:
        print("Tissue-aware mode: SKIPPED (no tissue cell line file)")
    
    all_passed = all(r for r in results if r is not None)
    print(f"\nOverall: {'ALL TESTS PASSED' if all_passed else 'SOME TESTS FAILED'}")
