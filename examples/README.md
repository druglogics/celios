# CELIOS Example: Getting Started

This folder contains a complete public example demonstrating the CELIOS pipeline using toy data. It serves as the primary onboarding experience for new CELIOS users.

## What This Example Demonstrates

CELIOS is a two-step pipeline for integrating omics data with biological networks:

1. **Node Extraction (Step 1)**: Build a mapping from gene symbols to network nodes using a biological interaction network (SIF format) and the HGNC gene database.

2. **Activity Calculation (Step 2)**: Extract and integrate omics data (gene expression, mutations, CNVs, TF activities) to compute activity scores for each network node across your samples.

This example uses a small toy network and toy omics data to demonstrate the full workflow. The output includes activity scores for nodes across cell lines, suitable for downstream analysis or machine learning.

## Folder Structure

```
examples/
├── README.md                  # This file
├── config.yaml               # Configuration file (source of truth)
├── celios_run_example.py     # Minimal Python API example
└── results/                  # Output directory (created by pipeline)
    ├── node_HGNC_dict.csv           # Mapping from network nodes to genes
    ├── identifiers.csv              # Cell line identifier resolution
    ├── activity_master_matrix.csv   # Full activity matrix (all nodes/cell lines)
    ├── activity_from_master.csv     # Selected activity values used in training
    ├── cell_lines/                  # Per-cell-line training files
    │   ├── COLO205/
    │   ├── DLD1/
    │   ├── HT29/
    │   └── SW620/
    └── run_log.txt            # Pipeline execution log
```

## Understanding config.yaml

The `config.yaml` file is the source of truth for all pipeline parameters. It has three main sections:

### Paths Section

Defines where input data and outputs are located:

```yaml
paths:
  base: "."                          # Base directory (where you run from)
  input: "../data/omics_input"      # Input data location
  output: "results"                 # Output directory
  cellfiles_dir: "results/cell_lines" # Per-cell-line files location
```

**Note**: Paths can be absolute or relative. When running from the examples/ directory, use `../data/` to reference the shared data folder.

### Steps Section: Node Dictionary Extraction

Step 1 builds a mapping from network nodes to gene symbols:

| Parameter | Purpose |
|-----------|---------|
| `node_input` | SIF network file (tab-separated triplets: source → interaction → target) |
| `hgnc_symbols_file` | HGNC gene database for standardizing gene symbols |
| `manual_symbols_file` | Optional manual overrides for genes not in HGNC |
| `include_alias_previous_symbols` | Whether to include HGNC aliases and previous symbols |
| `directory_output` | Where to save the node dictionary |

### Steps Section: Activity Calculation

Step 2 extracts and integrates omics data:

| Parameter | Purpose |
|-----------|---------|
| `cell_line_file` | CSV list of cell lines to analyze (with cell line names/IDs) |
| `activity_file` | Gene expression data (genes × samples, TPM or log2 format) |
| `mutations_file` | Binary mutation matrix (0=wild-type, 1=mutated) |
| `cnv_file` | Copy number variation matrix (0=loss, 1=neutral, 2=gain) |
| `tf_activity_file` | Optional transcription factor activity scores |
| `data_sources` | Which omics types to integrate: `["mutations", "cnv", "expression", "TF"]` |
| `directory_output` | Where to save activity results |

## Running CELIOS from the Command Line

### 1. Navigate to the examples directory

```bash
cd examples/
```

### 2. Run the Python example script

```bash
python celios_run_example.py
```

The script will:
- Read the configuration from `config.yaml`
- Extract the network node dictionary
- Process cell line identifiers
- Integrate omics data
- Generate activity scores
- Save results to `results/`

Expected output:
```
Running CELIOS example with toy data...
Config base: .../examples
Results will be saved to: results

STEP 1: Node dictionary
...
STEP 2: Resolving cell line identifiers
...
STEP 3: Extracting omics - activity matrix
...
STEP 3: Writing DrugLogic pipeline training files

Pipeline finished successfully!
```

## Running CELIOS from Python

For interactive workflows or integration into larger pipelines, use the Python API:

```python
import os
from celios.core import run_celios

# Define configuration
config = {
    "paths": {
        "base": ".",
        "input": "../data/omics_input",
        "output": "results",
        "cellfiles_dir": "results/cell_lines",
    },
    "steps": {
        "Node": {
            "node_input": "../data/node_dic_input/toy_model.sif",
            "hgnc_symbols_file": "../data/node_dic_input/hgnc_complete_set.txt",
            "manual_symbols_file": "../data/node_dic_input/manual_symbols.csv",
            "include_alias_previous_symbols": False,
            "directory_output": "results",
        },
        "Activity": {
            "cell_line_file": "../data/omics_input/toy_cells_list.csv",
            "activity_file": "../data/omics_input/rnaseq_tpm_coding_genes26Q1.csv",
            "mutations_file": "../data/omics_input/26Q1_mutations.csv",
            "cnv_file": "../data/omics_input/CCLE_CNV_binary.csv",
            "directory_output": "results",
            "data_sources": ["mutations", "cnv", "expression"],
        },
    },
}

# Run the pipeline
artifacts = run_celios(config=config, verbose=True)

# Access outputs
activity_matrix = artifacts["activity_matrix"]
node_dict = artifacts["node_dict"]
```

The `run_celios()` function returns a dictionary of artifacts:
- `celios_configuration`: The merged configuration used
- `node_dict`: Dictionary mapping nodes to genes
- `activity_matrix`: DataFrame with activity scores (nodes × cell lines)
- `pipeline_dir`: Path to results directory

### Advanced Options

```python
# Run with additional options
artifacts = run_celios(
    config=config,
    verbose=True,           # Print detailed progress
    plan=False,             # Set to True to just validate config (don't run)
    get_DLPfiles=True,      # Generate DrugLogics training files
    proliferation_state=True, # Include proliferation state in outputs
    dna_damage_state=False,  # Include DNA damage state
)
```

## Expected Outputs

After running the pipeline, you'll find the following files in the `results/` directory:

### Core Output Files

**`node_HGNC_dict.csv`**
- Maps network nodes to HGNC gene symbols
- Columns: `node_id`, `symbol`, `Entrez_ID` (optional)
- Use this to interpret node IDs in activity matrices

**`identifiers.csv`**
- Cell line identifier resolution mapping
- Shows how input cell line names map to SIDM identifiers
- Useful for tracking sample provenance

**`activity_master_matrix.csv`**
- Complete activity scores for all nodes across all cell lines
- Rows: network nodes
- Columns: cell line + data source combinations
- Format: `{cell_line}_{data_source}` for each omics type
- Includes: `{cell_line}_expression`, `{cell_line}_mutations`, `{cell_line}_cnv`, etc.
- Use this for comprehensive analysis and visualization

**`activity_from_master.csv`**
- Subset of activity matrix selected for downstream use
- Contains only the activity scores configured in `data_sources`
- More compact than the full master matrix
- Default rows/columns depend on pipeline configuration

### Cell Line Directories

**`cell_lines/{cell_line}_{data_type}_training/`**
- Per-cell-line training files for DrugLogics pipeline integration
- Contains preprocessed omics and activity data for machine learning
- One directory per cell line and data type combination
- Files are ready for downstream computational pipeline

### Pipeline Documentation

**`run_log.txt`**
- Detailed execution log
- Useful for debugging and verifying pipeline run parameters

## Visualizing Results

The companion notebook `notebooks/celios_visualize_example.ipynb` demonstrates how to visualize and explore the pipeline outputs:

### Quick Start

1. Navigate to the notebooks folder
2. Open `celios_visualize_example.ipynb` in Jupyter
3. Run all cells to generate visualizations

### What the Notebook Does

- Loads the activity matrix from `results/activity_master_matrix.csv`
- Generates a heatmap of activity scores across all cell lines and nodes
- Inspects individual cell line activity patterns
- Demonstrates how to select and visualize specific nodes or cell lines
- Shows how to extract data for downstream analysis

### Example Visualization: Node Activity Across Cell Lines

```python
import pandas as pd
import plotly.graph_objects as go

# Load the activity matrix
activity_matrix = pd.read_csv("results/activity_master_matrix.csv", index_col=0)

# Create a heatmap
fig = go.Figure(data=go.Heatmap(z=activity_matrix.values, x=activity_matrix.columns, y=activity_matrix.index))
fig.update_layout(title="Node Activity Across Cell Lines")
fig.show()
```

### Example: Inspect a Single Cell Line

```python
# Select a specific cell line
cell_line = "COLO205"
cell_line_data = activity_matrix[[col for col in activity_matrix.columns if cell_line in col]]
print(cell_line_data)
```

## Next Steps

1. **Understand the outputs**: Browse `results/activity_master_matrix.csv` and `results/node_HGNC_dict.csv`
2. **Visualize the results**: Run `notebooks/celios_visualize_example.ipynb`
3. **Adapt to your data**: Modify `config.yaml` with your own network and omics files
4. **Integrate downstream**: Use the activity matrices in your own analysis or machine learning pipelines

## Troubleshooting

### Error: "File not found"
- Check that paths in `config.yaml` are correct relative to where you're running the script
- From examples/ directory, data files should be accessible via `../data/`

### Error: "No matching cell lines"
- Verify cell line names in your activity file match those in `cell_line_file`
- Check the identifier resolution output: `results/identifiers.csv`

### Missing columns in output
- Ensure all omics files (activity, mutations, CNV, TF) have overlapping cell line columns
- Check `data_sources` in config.yaml matches available input files

### Empty activity matrix
- Verify cell lines in input files are aligned
- Check expression and other omics file formats are correct
- Review the pipeline output for coverage statistics

## Additional Resources

- Main README: See `../README.md` for CELIOS project overview
- Documentation: See `../documentation/` for detailed technical information
- Tests: See `../src/tests/` for unit tests and additional usage examples
