import os
import sys
import json

# Make repo src importable
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
SRC = os.path.join(ROOT, 'src')
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from celios.core import run_celios


def main():
    """
    Minimal example demonstrating the CELIOS Python API.
    
    This example runs the full pipeline (Node extraction + Activity calculation)
    using toy data to demonstrate all key features.
    """
    # Define the CELIOS configuration using toy data
    config = {
        "paths": {
            "base": os.path.dirname(os.path.abspath(__file__)),
            "input": os.path.join(ROOT, "data", "omics_input"),
            "output": "examples/results",
            "cellfiles_dir": "examples/results/cell_lines",
        },
        "steps": {
            # STEP 1: Build the node dictionary from the network
            "Node": {
                "node_input": os.path.join(ROOT, "data", "node_dic_input", "toy_model.sif"),
                "hgnc_symbols_file": os.path.join(ROOT, "data", "node_dic_input", "hgnc_complete_set.txt"),
                "manual_symbols_file": os.path.join(ROOT, "data", "node_dic_input", "manual_symbols.csv"),
                "include_alias_previous_symbols": False,
                "directory_output": "examples/results",
            },
            # STEP 2: Extract and process omics data
            "Activity": {
                "cell_line_file": os.path.join(ROOT, "data", "omics_input", "toy_cells_list.csv"),
                "activity_file": os.path.join(ROOT, "data", "omics_input", "rnaseq_tpm_coding_genes26Q1.csv"),
                "mutations_file": os.path.join(ROOT, "data", "omics_input", "26Q1_mutations.csv"),
                "cnv_file": os.path.join(ROOT, "data", "omics_input", "CCLE_CNV_binary.csv"),
                "tf_activity_file": os.path.join(ROOT, "data", "omics_input", "ccle_tf_activities.csv"),
                "format_override": "26Q1",
                "directory_output": "examples/results",
                "data_sources": ["mutations", "cnv", "TF", "expression"],
            },
        },
    }

    print('Running CELIOS example with toy data...')
    print(f'Config base: {config["paths"]["base"]}')
    print(f'Results will be saved to: {config["paths"]["output"]}\n')
    
    # Execute the pipeline
    artifacts = run_celios(config=config, verbose=True)

    print('\n' + '='*60)
    print('Pipeline finished successfully!')
    print('='*60)
    
    if artifacts is None:
        print('ERROR: No artifacts returned')
        return 1

    print('\nGenerated artifacts:')
    for name, obj in artifacts.items():
        print(f'  ✓ {name}: {type(obj).__name__}')

    print(f'\nResults saved to: {config["paths"]["output"]}/')
    print('\nNext step: Open notebooks/celios_visualize_example.ipynb to visualize the results')
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
