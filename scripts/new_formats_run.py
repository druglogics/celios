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
    config = {
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
                "cell_line_file": "vis_2024/cell_line_list.csv",
                "mutations_file": "activity_input/CCLE_muts_binary.csv",
                "cnv_file": "activity_input/CCLE_CNV_binary.csv",
                "directory_output": "results",
                "data_sources": ["mutations", "cnv"],
            },
        },
    }

    print('Running CELIOS mutations/CNV-only pipeline...')
    artifacts = run_celios(config=config, verbose=True)

    print('\nPipeline finished. Artifacts:')
    if artifacts is None:
        print('No artifacts returned')
        return 1

    print(json.dumps({k: str(type(v)) for k, v in artifacts.items()}, indent=2))

    return 0


if __name__ == '__main__':
    sys.exit(main())
