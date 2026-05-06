CELIOS

CEll LIne OmicS processor for extracting omics data into integrated activity datasets, from which calibration files can be created for Boolean models used in the DrugLogics and the TRAFIKK pipelines.

> **New to CELIOS?** Start with [QUICKSTART.md](QUICKSTART.md) for a 5-minute introduction, or see [INSTALL.md](INSTALL.md) for detailed setup instructions.

## 📚 Documentation Map

| Document | Purpose |
|----------|---------|
| **[QUICKSTART.md](QUICKSTART.md)** | 5-minute quick reference with common commands |
| **[INSTALL.md](INSTALL.md)** | Installation guide, virtual environments, troubleshooting |
| **[PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md)** | Code organization, module descriptions, architecture |
| **[README.md](README.md)** | Full API reference, usage examples, configuration details (you are here) |
| **[notebooks/celios_visualize_bashi_final.ipynb](notebooks/celios_visualize_bashi_final.ipynb)** | Interactive example: tissue selection, CELIOS execution, visualization |

---

## Usage

### Running the Full Pipeline

#### Via CLI (recommended)

Run the pipeline from the repository root with a configuration file (JSON or YAML):

```cmd
python -m celios.cli run --config path\to\config.yaml --verbose
```

Options:
- `--config` (required): Path to pipeline config (JSON or YAML)
- `--verbose`: Print detailed output
- `--plan`: Only print execution plan without running
- `--stop-after`: Stop execution after a specific step

#### Via Python Script

Alternatively, define the configuration directly in a Python script and call `run_celios()`:

```python
from celios.core import run_celios

config = {
    "paths": {
        "base": ".",
        "input": "data",
        "output": "results",
        "cellfiles_dir": "results/cell_lines",
    },
    "steps": {
        "Node": {
            "node_input": "node_dic_input/DNAdamage.sif",
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
            "data_sources": ["mutations", "cnv", "TF"],
        },
    },
}

artifacts = run_celios(config=config, plan=False, verbose=True)
```

See `tests/test_run_celios.py` and `tests/celios_consensus.py` for additional examples.

### Tissue-Organized Output (Optional)

CELIOS supports organizing DrugLogics training files by tissue type. When `paths.tissue_dir` is specified, training files are written to `tissue_dir/<Tissue>/<cell_line_name>/` based on the tissue information in `cell_line_file`.

This is useful when:
- You have cell lines from multiple tissues and want organized output
- Existing tissue folders contain additional files (they will be preserved)
- You want to create new tissue/cell-line directories automatically

Example config:
```yaml
paths:
  tissue_dir: "results/tissue_folders"  # Enable tissue-organized output
steps:
  Activity:
    cell_line_file: "data/cell_line_list.csv"  # Must contain 'tissue', 'SIDM', and 'cell_line_name' columns
```

The CSV file must include columns for tissue, SIDM (unique identifier), and cell_line_name (folder names) when using tissue-organized output. For legacy mode (`cellfiles_dir`), no specific columns are required.

### Skipping the Node Step (Using Pre-built Node Dictionary)

If you already have a pre-built node dictionary file, you can skip **STEP 1 (Node dictionary generation)** by:

1. **Omitting** the `"Node"` section from `steps` in your config
2. **Adding** `"node_dic"` to the `Activity` section pointing to your CSV file

Example:

```python
config = {
    "paths": {
        "base": ".",
        "input": "consensus",
        "output": "consensus/hgsoc_results",
    },
    "steps": {
        "Activity": {
            "node_dic": "hgsoc_net/NODE_HGNC_equivalences.csv",  # Pre-built node dictionary
            "activity_file": "activity_input/rnaseq_tpm_20220624.csv",
            "cell_line_file": "activity_input/cell_line_list.csv",
            "tf_activity_file": "activity_input/ccle_tf_activities.csv",
            "mutations_file": "activity_input/CCLE_muts_binary.csv",
            "cnv_file": "activity_input/CCLE_CNV_binary.csv",
            "directory_output": "consensus/hgsoc_results",
            "data_sources": ["mutations", "cnv", "TF"],
        },
    },
}

artifacts = run_celios(config=config, verbose=True)
```

The pipeline will:
- **Skip STEP 1** (Node dictionary generation) since `"Node"` is not defined
- Load your CSV file as the node dictionary
- Proceed to **STEP 2** (Activity extraction)

Your CSV file should have at least two columns:
- First column: node names
- Second column: symbols (comma-separated list of gene symbols)

See `tests/celios_hgsoc.py` for a real example using a pre-built node dictionary.

### Pipeline Behavior

**Scenario 1: `"Node"` step is defined in config**
- STEP 1 will run (Node dictionary generation), even if `node_dic` is also in Activity
- The generated node dictionary will be used for Activity extraction

**Scenario 2: `"Node"` step is NOT defined, but `node_dic` is in Activity config**
- STEP 1 is skipped
- The provided CSV file is loaded and used directly

**Scenario 3: Neither `"Node"` step nor `node_dic` is provided**
- The pipeline will raise an error (no node dictionary source available)

### Feature Helpers

Call Node helpers directly with CLI subcommands:

```cmd
python -m celios.cli node-from-sif --sif examples\DNAdamage.sif --hgnc examples\hgnc_complete_set.txt --out results\node_dict.csv

python -m celios.cli node-from-object --input "TP53,BRCA1,EGFR" --hgnc examples\hgnc_complete_set.txt --out results\node_dict.csv --include_alias_prev
```

## Installation

### Install from PyPI

Once published, install CELIOS directly with:

```bash
pip install celios
```

### Development Install

If you are developing CELIOS from a local checkout:

```bash
pip install -e .
```

For detailed installation instructions, troubleshooting, and virtual environment setup, see [INSTALL.md](INSTALL.md).

### Run from Source

Alternatively, run CLI commands from the repository root without installing:

```bash
python -m celios.cli run --config path\to\config.yaml --verbose
```

---

## Architecture & Design

### Pipeline Overview

CELIOS is a two-step pipeline:

1. **Node Extraction** - Extract nodes from a biological network (SIF format) and map them to standardized gene symbols (HGNC)
2. **Activity Calculation** - Integrate multi-omics data (mutations, CNV, TF activity, gene expression) into activity matrices by cell line

Each step can be skipped or customized via configuration.

### Configuration-Driven Design

The pipeline is entirely controlled via JSON or YAML configuration files, making it easy to:
- Define data sources and output locations
- Skip pipeline steps
- Reproduce analyses
- Scale to multiple datasets

See [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) for detailed code organization and module descriptions.

## Notes

- **YAML support** - YAML config files require `pyyaml` in your environment (optional). JSON configs work without extra packages.
  ```bash
  pip install pyyaml
  ```
- **Example configurations** - See `src/tests/test_run_celios.py` and the notebooks for real-world examples
- **Interactive tutorials** - See `notebooks/1_select_visualize.ipynb` for a step-by-step walkthrough
- **Project structure** - See [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) for code organization
This works as long as you're in the repository directory.

## Release Notes

Maintainers should publish releases through GitHub Actions using the workflow in [.github/workflows/publish.yml](.github/workflows/publish.yml). The recommended path is to publish to TestPyPI first, validate `pip install celios` from that index, and then publish the tagged GitHub release to PyPI.

