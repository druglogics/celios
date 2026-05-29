# CELIOS Usage Guide

Comprehensive guide to running the CELIOS pipeline, including all features and advanced options.

---

## Table of Contents

1. [Running the Pipeline](#running-the-pipeline)
2. [Configuration Guide](#configuration-guide)
3. [Input File Formats](#input-file-formats)
4. [Output Files](#output-files)
5. [Advanced Features](#advanced-features)
6. [Python API](#python-api)
7. [Troubleshooting](#troubleshooting)

---

## Running the Pipeline

### Via Command Line (Recommended)

```bash
celios run --config config.yaml --verbose
```

**Options:**
- `--config` (required): Path to JSON or YAML config file
- `--verbose`: Print detailed execution output
- `--plan`: Show execution plan without running pipeline
- `--stop-after`: Stop execution after a specific step (`node`, `activity`)

**Examples:**

```bash
# Run full pipeline with details
celios run --config analysis/config.yaml --verbose

# Preview what will execute
celios run --config analysis/config.yaml --plan

# Stop after node dictionary generation
celios run --config analysis/config.yaml --stop-after node

# Stop after cell-line resolution (before activity extraction)
celios run --config analysis/config.yaml --stop-after activity
```

### Via Python Script

```python
from celios.core import run_celios

config = {
    "paths": {
        "base": ".",
        "input": "data",
        "output": "results",
    },
    "steps": {
        "Node": {...},
        "Activity": {...},
    }
}

artifacts = run_celios(config=config, verbose=True, plan=False)

# Access results
print(f"Nodes processed: {len(artifacts['node_dictionary'])}")
print(f"Activity matrix shape: {artifacts['activity_matrix'].shape}")
```

---

## Configuration Guide

CELIOS is configured via JSON or YAML files. Below are complete examples for different scenarios.

### Minimal Configuration (Skip Node Step)

Use this if you have a pre-built node dictionary:

```yaml
paths:
  base: "."
  input: "data"
  output: "results"

steps:
  Activity:
    node_dic: "data/node_HGNC_dict.csv"
    activity_file: "data/rnaseq_tpm_20220624.csv"
    cell_line_file: "data/cell_line_list.csv"
    directory_output: "results"
```

### Full Configuration (Include All Steps)

```yaml
paths:
  base: "."
  input: "data"
  output: "results"
  cellfiles_dir: "results/cell_lines"           # Optional: per-cell-line training files
  tissue_dir: "results/tissue_folders"          # Optional: tissue-organized output

steps:
  Node:
    node_input: "data/network.sif"
    hgnc_symbols_file: "data/hgnc_complete_set.txt"
    manual_symbols_file: "data/manual_symbols.csv"     # Optional
    include_alias_previous_symbols: false
    directory_output: "results"

  Activity:
    activity_file: "data/rnaseq_tpm_20220624.csv"
    cell_line_file: "data/cell_line_list.csv"
    tf_activity_file: "data/tf_activities.csv"         # Optional
    mutations_file: "data/mutations.csv"               # Optional
    cnv_file: "data/cnv.csv"                           # Optional
    directory_output: "results"
    data_sources: ["mutations", "cnv", "TF", "expression"]
```

### Configuration with Format Overrides

CELIOS auto-detects file formats, but you can explicitly specify:

```yaml
steps:
  Activity:
    activity_file: "data/rnaseq_26Q1.csv"
    format_override: "26q1"                    # Force 26Q1 format for expression
    
    mutations_file: "data/mutations_modelid.csv"
    mutations_format_override: "26q1"          # Force 26Q1 format for mutations
    
    cnv_file: "data/cnv_sidm.csv"
    cnv_format_override: "old"                 # Force old format for CNV
```

### Configuration with Tissue Organization

```yaml
paths:
  base: "."
  tissue_dir: "results/tissue_folders"

steps:
  Activity:
    activity_file: "data/activity.csv"
    cell_line_file: "data/cell_line_list.csv"  # MUST include: tissue, SIDM, cell_line_name
    directory_output: "results"
```

Cell line file columns required for tissue organization:
- `tissue`: Tissue type (e.g., "breast", "lung")
- `SIDM`: Sample identifier (unique)
- `cell_line_name`: Folder-friendly name (e.g., "MCF7")

---

## Input File Formats

### 1. Network File (SIF Format)

**File:** Network definition in SIF (Simple Interaction Format)

**Format:** Tab-separated, 3 columns:
```
source   interaction_type   target
TP53     p>                 TP53_mutant
BRCA1    i>                 PARP
KLF4     i>                 KLF4_activity
```

**Notes:**
- One interaction per line
- Duplicate nodes are automatically deduplicated
- Comments (lines starting with #) are ignored
- Nodes extracted: all unique values in source and target columns

### 2. HGNC Symbols File

**File:** HGNC gene symbol reference

**Format:** Tab-separated, standard HGNC format:
```
HGNC ID  Symbol  Previous Symbols  Alias Symbols
HGNC:11857   TP53   P53             ...
```

**Notes:**
- Used to map network nodes to standardized gene symbols
- Column order: HGNC ID, Symbol, Previous Symbols, Alias Symbols (typical)
- `include_alias_previous_symbols` flag determines whether to use aliases/previous names

### 3. Manual Symbols File (Optional)

**File:** Corrections for nodes that don't map correctly

**Format:** CSV with 2 columns:
```
node,HGNC_symbol
custom_node_1,ENSG00000141510
TP53_mutant,TP53
```

**Notes:**
- Applied after HGNC mapping
- Overrides automatic mappings

### 4. Cell Line File

**Format:** CSV with columns:

For **legacy mode** (cellfiles_dir):
- No specific columns required
- Any columns used for reference

For **tissue-organized mode** (tissue_dir):
- `tissue` (required): Tissue type
- `SIDM` (required): Unique sample identifier
- `cell_line_name` (required): Folder name

**Example:**
```
tissue,SIDM,cell_line_name,source
breast,SIDM00456,MCF7,CCLE
lung,SIDM00123,H1299,CCLE
```

### 5. Activity File (Gene Expression)

**Format:** Two supported formats, auto-detected.

#### Legacy Format (Old)
- **Index:** Gene symbols
- **Columns:** SIDM identifiers
- **Values:** Gene expression (TPM, RPKM, or normalized)
- **File:** `rnaseq_tpm_20220624.csv`

Example:
```
,SIDM00001,SIDM00002,SIDM00003
TP53,5.2,3.1,0.0
BRCA1,1.5,2.3,4.1
```

#### 26Q1 Format (New)
- **Index:** Gene symbols
- **Columns:** DepMap ModelID (ACH-*)
- **Values:** Gene expression
- **File:** `rnaseq_tpm_coding_genes26Q1.csv`
- **ModelID mapping:** Automatically mapped to SIDM via `src/celios/features/Model.csv`

Example:
```
,ACH-000001,ACH-000002,ACH-000003
TP53,5.2,3.1,0.0
BRCA1,1.5,2.3,4.1
```

### 6. Binary Matrices (Mutations, CNV)

**Supported Formats:**

#### Old Format (Genes × SIDM)
- **Index:** Gene symbols
- **Columns:** SIDM identifiers
- **Values:** 0/1 (binary)

#### 26Q1 Format (ModelID × Genes)
- **Index:** DepMap ModelID (ACH-*)
- **Columns:** Gene symbols
- **Values:** 0/1 (binary)
- **ModelID mapping:** Automatically converted to SIDM

**Example (Old):**
```
,SIDM00001,SIDM00002
TP53,1,0
BRCA1,0,1
```

**Example (26Q1):**
```
,TP53,BRCA1
ACH-000001,1,0
ACH-000002,0,1
```

---

## Output Files

### Outputs from Step 1: Node Dictionary Generation

**File:** `node_HGNC_dict.csv`

**Format:** CSV with 2 columns
```
node_name,symbols
TP53,"TP53, P53"
BRCA1,BRCA1
```

**Use:** Input to Step 3 (Activity extraction). Defines nodes that will be in final activity matrix.

### Outputs from Step 2: Cell-Line Mapping

**File:** `identifiers.csv` (internal reference)

**Format:** CSV with SIDM → cell_line_name mappings

### Outputs from Step 3: Omics Integration

#### Master Activity Matrix

**File:** `activity_master_matrix.csv`

**Format:** CSV, nodes × multi-source activities
- **Rows:** Node names (from node dictionary)
- **Columns:** `{cell_line}__{data_source}` (e.g., `MCF7__mutations`, `MCF7__expression`)
- **Values:** Activity scores (0/1 for binary, continuous for TF/expression)

**Use:** View all omics data sources, debug, generate custom priority rankings

#### Priority-Selected Activity

**File:** `activity_from_master.csv`

**Format:** CSV, nodes × cell_lines
- **Rows:** Node names
- **Columns:** Cell line names
- **Values:** Single activity value per node-cell line (highest priority source)

**Priority order (default):** mutations → CNV → TF → expression
- Each node-cell combination uses the first available data source
- Can customize with `data_sources` configuration

**Use:** Ready for Boolean model calibration, DrugLogics compatibility

#### Per-Cell-Line Training Files (Optional)

**Location:** `cellfiles_dir/<cell_line_name>.csv` (legacy mode)
**Or:** `tissue_dir/<Tissue>/<CellLine>/training.csv` (tissue-organized mode)

**Format:** CSV, nodes × 1 (activity values for single cell line)

**Use:** One file per cell line, ready for DrugLogics pipeline

### Pipeline Log

**File:** `run_log.txt`

**Contents:**
- Execution timestamp
- Configuration summary
- Step execution status
- File statistics (nodes processed, cell lines, dimensions)
- Data source diagnostics
- Warnings and errors

---

## Advanced Features

### 1. Skipping Node Generation

If you have a pre-built node dictionary:

```yaml
paths:
  base: "."
  input: "data"
  output: "results"

steps:
  Activity:
    node_dic: "data/my_nodes.csv"     # Pre-built node dictionary
    activity_file: "data/activity.csv"
    cell_line_file: "data/cells.csv"
    directory_output: "results"
```

**Note:** Omit the `Node` step entirely. CELIOS will load your CSV file as-is.

### 2. Tissue-Organized Output

Organize outputs by tissue/cell-line hierarchy:

```yaml
paths:
  base: "."
  tissue_dir: "results/tissue_folders"

steps:
  Activity:
    activity_file: "data/activity.csv"
    cell_line_file: "data/cell_line_list.csv"  # Must have: tissue, SIDM, cell_line_name
    directory_output: "results"
```

**Result:**
```
tissue_folders/
├── breast/
│   ├── MCF7/
│   │   ├── training.csv
│   │   └── metadata.txt
│   └── T47D/
│       └── training.csv
├── lung/
│   └── H1299/
│       └── training.csv
```

### 3. Custom Data Source Priority

Change the priority order for activity selection:

```yaml
steps:
  Activity:
    activity_file: "data/activity.csv"
    cell_line_file: "data/cells.csv"
    mutations_file: "data/mutations.csv"
    cnv_file: "data/cnv.csv"
    tf_activity_file: "data/tf.csv"
    data_sources: ["TF", "mutations", "expression"]  # Custom priority
    directory_output: "results"
```

**Result:** For each node-cell pair, use TF activity first, then mutations, then expression.

### 4. Format Auto-Detection with Overrides

CELIOS detects formats automatically, but you can force specific formats:

```yaml
steps:
  Activity:
    # Force all formats explicitly
    activity_file: "data/activity.csv"
    format_override: "26q1"
    
    mutations_file: "data/muts.csv"
    mutations_format_override: "old"
    
    cnv_file: "data/cnv.csv"
    cnv_format_override: "26q1"
```

**Formats:** `"old"` or `"26q1"` or `null` (auto-detect)

### 5. Manual Node Corrections

Override automatic HGNC mappings:

```yaml
steps:
  Node:
    node_input: "data/network.sif"
    hgnc_symbols_file: "data/hgnc_complete_set.txt"
    manual_symbols_file: "data/corrections.csv"  # Your overrides
    directory_output: "results"
```

**Format of corrections.csv:**
```
node,HGNC_symbol
custom_kinase_1,MAP2K1
TP53_mutant,TP53
```

---

## Python API

### Main Entry Point

```python
from celios.core import run_celios

artifacts = run_celios(
    config=config_dict,
    plan=False,                 # Set True to preview only
    verbose=True,               # Detailed output
    stop_after=None             # "node" or "activity" to stop early
)
```

### Return Value

`artifacts` dictionary contains:

```python
{
    'node_dictionary': dict,           # Node → symbols mapping
    'activity_matrix': pd.DataFrame,   # Master activity matrix
    'activity_final': pd.DataFrame,    # Priority-selected activities
    'cell_lines': list,                # Cell line identifiers
    'node_count': int,                 # Number of nodes
    'run_log': str,                    # Execution log
}
```

### Example: Custom Analysis

```python
from celios.core import run_celios

config = {...}  # Your configuration

# Run pipeline
results = run_celios(config=config, verbose=True)

# Access outputs
nodes_dict = results['node_dictionary']
activity_matrix = results['activity_matrix']
final_activity = results['activity_final']

# Analyze
print(f"Total nodes: {results['node_count']}")
print(f"Cell lines processed: {len(results['cell_lines'])}")
print(f"Activity matrix shape: {activity_matrix.shape}")

# Save custom analysis
final_activity.to_csv("my_analysis.csv")
```

### Helper Functions

```python
from celios.features.node import Node

# Generate node dictionary from SIF file
node_dict, mapped = Node.from_sif(
    sif_path="data/network.sif",
    hgnc_symbols_file="data/hgnc.txt",
    directory_output="results",
    include_alias_prev=False
)

# Generate from gene list
node_dict2, mapped2 = Node.from_object(
    input_obj=["TP53", "BRCA1", "EGFR"],
    hgnc_symbols_file="data/hgnc.txt",
    include_alias_prev=False
)
```

---

## Troubleshooting

### Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `FileNotFoundError: config.yaml` | Config file not found | Use absolute path or verify file exists |
| `ModuleNotFoundError: No module named 'celios'` | Package not installed or not in PATH | Run: `pip install -e .` from repo root |
| `celios: command not found` | Console script not registered | Reinstall: `pip install -e .` |
| `YAML config not recognized` | PyYAML not installed | Run: `pip install pyyaml` |
| `Cannot auto-detect binary matrix format` | Unrecognized file headers | Use explicit format_override |
| `No ModelID values found in Model.csv` | Invalid ModelID format | Ensure IDs start with "ACH-" for 26Q1 format |
| `ValueError: No node dictionary source` | Missing both Node step and node_dic | Either add Node step or provide node_dic path |
| `KeyError: 'tissue'` | Missing required column in cell_line_file | Add columns: tissue, SIDM, cell_line_name |

### Debug Mode

Run with verbose output and execution plan:

```bash
# See what will run
celios run --config config.yaml --plan

# Run with detailed output
celios run --config config.yaml --verbose

# Stop after specific step for inspection
celios run --config config.yaml --stop-after node
ls -la results/
```

### Check Installation

```bash
# Verify CELIOS is installed
celios --help

# Test basic functionality
celios node-from-object --input "TP53" --hgnc data/hgnc.txt --out /tmp/test_node.csv

# Verify dependencies
python -c "import pandas; import yaml; print('OK')"
```

---

## FAQ

**Q: How do I use pre-built node dictionaries?**
A: Omit the Node step and provide `steps.Activity.node_dic` path to your CSV file.

**Q: Can I run just the node generation step?**
A: Yes: `celios run --config config.yaml --stop-after node`

**Q: What's the difference between activity_master_matrix and activity_from_master?**
A: Master includes all sources (columns per source). From_master is priority-selected (one value per node-cell).

**Q: How do I add more data sources?**
A: Pass additional files (mutations, cnv, tf_activity) in config and add to data_sources list.

**Q: Can I customize the priority order?**
A: Yes, reorder `data_sources` list in Activity config.

---

For more help, see [PIPELINE.md](PIPELINE.md) for architecture details or [INSTALL.md](INSTALL.md) for installation issues.
