# 📂 Project Structure

Understanding the CELIOS codebase organization.

---

## 📁 Repository Layout

```
CELIOS/
├── 📄 README.md                    # Main documentation & usage guide
├── 📄 INSTALL.md                   # Installation instructions & troubleshooting
├── 📄 QUICKSTART.md                # 5-minute quick start guide
├── 📄 PROJECT_STRUCTURE.md         # This file - code organization
├── 📄 LICENSE                      # License information
├── 📄 setup.py                     # Package configuration & dependencies
├── 📄 requirements.txt             # Python dependencies list
│
├── 📁 celios/                      # Main package directory
│   ├── __init__.py                 # Package initialization
│   ├── base_defaults.py            # Default configuration constants
│   ├── cli.py                      # Command-line interface entry point
│   ├── core.py                     # Core pipeline orchestration logic
│   │
│   ├── 📁 features/                # Feature modules for pipeline steps
│   │   ├── __init__.py
│   │   ├── activity.py             # Activity extraction step
│   │   ├── files.py                # File I/O operations
│   │   ├── node.py                 # Node dictionary generation step
│   │   ├── sifbase.py              # SIF network parsing utilities
│   │   ├── tissue.py               # Tissue-aware file organization
│   │   └── training.py             # Training file generation (DrugLogics format)
│   │
│   ├── 📁 utils/                   # Utility modules
│   │   ├── __init__.py
│   │   ├── io.py                   # Input/output helper functions
│   │   ├── report.py               # Report generation and logging
│   │   └── validate.py             # Configuration & data validation
│   │
│   └── 📁 __pycache__/             # Python cache (auto-generated)
│
├── 📁 data/                        # Sample and reference data
│   ├── 📁 activity_input/          # Example activity matrices
│   │   ├── CCLE_CNV_binary.csv
│   │   ├── CCLE_muts_binary.csv
│   │   ├── ccle_tf_activities.csv
│   │   └── rnaseq_tpm_20220624.csv
│   │
│   ├── 📁 node_dic_input/          # Network & annotation files
│   │   ├── cell_fate_plus.sif      # Network file (SIF format)
│   │   ├── hgnc_complete_set.txt   # HGNC gene symbol reference
│   │   └── manual_symbols.csv      # Manual gene symbol overrides
│   │
│   ├── 📁 oneil_2016/              # O'Neill 2016 dataset
│   │   ├── oneil_cells.csv
│   │   ├── oneil_cells_sidm.csv
│   │   └── 📁 results/             # Example outputs from O'Neill analysis
│   │
│   └── 📁 vis_2024/                # Visualization 2024 dataset
│       ├── cell_line_list.csv
│       ├── S1_cell_line_drug_info.xlsx
│       └── 📁 results/             # Example outputs from Vis 2024 analysis
│
├── 📁 notebooks/                   # Jupyter notebooks
│   └── 1_select_visualize.ipynb    # Interactive example: tissue selection & visualization
│
├── 📁 tests/                       # Test suite
│   ├── test_run_celios.py          # Main pipeline tests
│   └── __pycache__/                # Test cache (auto-generated)
│
├── 📁 celios.egg-info/             # Package metadata (auto-generated)
│   ├── entry_points.txt
│   ├── requires.txt
│   └── ...
│
└── 📁 siflex_example/              # Reference documentation style (for maintainers)
    ├── README.md
    ├── INSTALL.md
    ├── QUICKSTART.md
    └── PROJECT_STRUCTURE.md
```

---

## 🔧 Core Modules

### **celios/core.py** (~200 lines)
**Purpose:** Pipeline orchestration and configuration management

**Key Functions:**
- `run_celios(config, **kwargs)` - Main entry point; merges config with defaults, executes Node and Activity steps
- `load_config(path)` - Loads JSON/YAML configuration files
- `merge_config(user_config, defaults)` - Merges user config with defaults

**Key Features:**
- Handles both JSON and YAML configurations
- Validates configuration before execution
- Orchestrates Node → Activity pipeline
- Returns unified artifact dictionary

---

### **celios/cli.py** (~150 lines)
**Purpose:** Command-line interface for CELIOS

**Key Functions:**
- `main()` - CLI entry point (called by `celios` console script)
- `run()` - Subcommand for full pipeline execution
- `node_from_sif()` - Subcommand to generate nodes from SIF file
- `node_from_object()` - Subcommand to generate nodes from gene list

**Key Features:**
- Click-based CLI framework
- Arguments: `--config`, `--verbose`, `--plan`, `--stop-after`
- Direct node generation helpers (bypasses full pipeline)
- Integrates with core.py

---

### **celios/features/node.py** (~300 lines)
**Purpose:** Node dictionary generation from network SIF file

**Key Functions:**
- `extract_nodes_from_sif(sif_path)` - Parse SIF network; extract unique nodes
- `map_nodes_to_hgnc(nodes, hgnc_file)` - Map nodes to HGNC gene symbols
- `apply_manual_overrides(node_dict, manual_file)` - Apply manual symbol corrections
- `validate_nodes(node_dict)` - Validate node-to-symbol mappings

**Key Features:**
- Handles SIF network format (tab-separated: source → target → interaction_type)
- HGNC symbol resolution with fuzzy matching
- Support for alias and previous symbols (toggle via `include_alias_previous_symbols`)
- Manual overrides via CSV file (columns: `node`, `HGNC_symbol`)
- Outputs: `node_HGNC_dict.csv`

---

### **celios/features/activity.py** (~400 lines)
**Purpose:** Activity matrix extraction from multi-omics data

**Key Functions:**
- `build_activity_matrix(cell_line_file, activity_file, tf_file, mutations_file, cnv_file, data_sources)` - Main orchestrator
- `load_activity_data(activity_file, cell_lines)` - Load and filter gene expression data
- `load_omics_data(mutations_file, cnv_file, tf_file, cell_lines)` - Load binary/continuous omics matrices
- `merge_data_sources(data_dict, priority)` - Merge multi-source data with priority ordering
- `generate_training_files(activity_matrix, cell_lines, output_dir)` - Create per-cell-line DrugLogics files

**Key Features:**
- Supports multiple data sources: mutations (binary), CNV (binary), TF activity (continuous), gene expression (continuous)
- Priority-based data selection (mutations → CNV → TF → expression by default)
- Handles multi-index cell lines (SIDM + cell_line_name)
- Outputs: `activity_master_matrix.csv` (all sources), `activity_from_master.csv` (priority-selected)
- Optional: per-cell-line training files in DrugLogics format

---

### **celios/features/sifbase.py** (~100 lines)
**Purpose:** SIF network parsing utilities

**Key Functions:**
- `parse_sif(sif_path)` - Parse SIF file into edge list
- `extract_unique_nodes(edges)` - Get unique nodes from edge list
- `validate_sif_format(sif_path)` - Verify SIF file structure

**Key Features:**
- Handles tab-separated SIF format
- Supports comments and blank lines
- Error handling for malformed files

---

### **celios/utils/io.py** (~150 lines)
**Purpose:** File I/O operations and data loading

**Key Functions:**
- `read_csv(path, **kwargs)` - Wrapper around pandas.read_csv with error handling
- `write_csv(df, path, **kwargs)` - Write DataFrame to CSV with metadata
- `read_hgnc_symbols(hgnc_file)` - Load HGNC gene symbol reference
- `load_cell_lines(cell_line_file)` - Load cell line metadata (SIDM, cell_line_name)

**Key Features:**
- Consistent error messages for missing files
- Automatic path resolution (absolute/relative)
- Support for multi-index and wide-format DataFrames

---

### **celios/utils/validate.py** (~100 lines)
**Purpose:** Configuration and data validation

**Key Functions:**
- `validate_config(config)` - Validate config structure and required fields
- `validate_paths(config)` - Check that all file paths exist
- `validate_activity_data(activity_matrix)` - Validate activity matrix format

**Key Features:**
- Early error detection before pipeline execution
- Clear error messages for missing or invalid configuration

---

### **celios/utils/report.py** (~80 lines)
**Purpose:** Report generation and logging

**Key Functions:**
- `log_step(step_name, status, message)` - Log pipeline step execution
- `save_run_report(artifacts, output_dir)` - Save execution report with statistics
- `print_summary(artifacts)` - Print pipeline summary to console

**Key Features:**
- Pipeline execution logging
- Artifact statistics (nodes processed, cell lines, data sources)
- Run timestamps and durations

---

### **celios/features/files.py** (~100 lines)
**Purpose:** File format handling and conversion

**Key Functions:**
- `ensure_output_dir(dir_path)` - Create output directory if not exists
- `save_artifacts(artifacts, output_dir)` - Save all pipeline outputs to disk
- `cleanup_temp_files()` - Remove temporary files after execution

**Key Features:**
- Atomic file operations (prevent partial writes)
- Backup of existing outputs
- Clear file naming conventions

---

### **celios/features/training.py** (~150 lines)
**Purpose:** Generate DrugLogics-compatible training files

**Key Functions:**
- `generate_training_files(activity_matrix, cell_lines, output_dir)` - Create per-cell-line files
- `format_training_data(activity_vector)` - Convert activity vector to DrugLogics format
- `validate_training_format(training_file)` - Verify DrugLogics format compliance

**Key Features:**
- One training file per cell line
- DrugLogics-compatible format (nodes, activities, metadata)
- Metadata embedding (data sources, date, parameters)

---

### **celios/features/tissue.py** (~80 lines)
**Purpose:** Tissue-aware organization of DrugLogics training files

**Key Functions:**
- `write_tissue_files(activity_df, tissue_dir, cell_line_file, ...)` - Organize files by tissue/cell line

**Key Features:**
- Creates directory structure: `tissue_dir/Tissue/CellLine/`
- Maps cell lines to tissues using SIDM identifiers
- Applies cell line name cleaning for folder consistency
- Integrates with existing DL training file generation

---

## 🎯 Quick File Lookup

**I need to...**

| Task | File(s) |
|------|---------|
| **Run the pipeline** | `celios/cli.py`, `celios/core.py` |
| **Modify CLI commands** | `celios/cli.py` |
| **Change Node step logic** | `celios/features/node.py`, `celios/features/sifbase.py` |
| **Change Activity step logic** | `celios/features/activity.py` |
| **Modify output format** | `celios/features/files.py`, `celios/utils/io.py` |
| **Update configuration validation** | `celios/utils/validate.py` |
| **Change default parameters** | `celios/base_defaults.py` |
| **Add new tests** | `tests/test_run_celios.py` |
| **Improve documentation** | `README.md`, `INSTALL.md`, `QUICKSTART.md` |

---

## 📊 Project Statistics

| Metric | Value |
|--------|-------|
| **Python Files** | 13 |
| **Total Lines of Code** | ~1,800 |
| **Main Modules** | 8 |
| **Test Files** | 1 |
| **Primary Dependency** | pandas ≥ 1.0.0 |
| **Optional Dependencies** | pyyaml (for YAML configs) |
| **Python Version** | 3.8+ |
| **License** | See LICENSE |

---

## 🔄 Data Flow Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│ INPUT FILES                                                     │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  Network:          Gene Symbols:      Gene Data:               │
│  cell_fate.sif  +  hgnc_symbols.txt + rnaseq.csv      +        │
│                    manual_symbols.csv   mutations.csv          │
│                                        cnv.csv                 │
│                                        tf_activities.csv       │
│                                        cell_lines.csv          │
└─────────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────────┐
│ STEP 1: NODE EXTRACTION (celios/features/node.py)              │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  1. Parse SIF network                                           │
│  2. Extract unique nodes from network topology                 │
│  3. Map nodes → HGNC gene symbols (fuzzy matching)             │
│  4. Apply manual overrides                                      │
│  5. Validate mappings                                           │
│                                                                 │
│  Output: node_HGNC_dict.csv (node ↔ gene symbol mapping)       │
└─────────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────────┐
│ STEP 2: ACTIVITY EXTRACTION (celios/features/activity.py)      │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  1. Load cell lines from metadata (SIDM, cell_line_name)       │
│  2. Load gene expression data, filter to nodes                 │
│  3. Load binary mutation/CNV data, binarize                    │
│  4. Load TF activity scores                                    │
│  5. Merge all data sources by priority:                        │
│     mutations → CNV → TF → expression                          │
│  6. Create activity matrix (nodes × cell_lines__source)        │
│  7. Generate priority-selected matrix (1 source per cell)      │
│  8. Generate per-cell-line training files (optional)           │
│                                                                 │
│  Outputs:                                                       │
│    - activity_master_matrix.csv (all sources)                  │
│    - activity_from_master.csv (priority-selected)              │
│    - tissue_folders/ (optional, tissue-organized)              │
│      ├── Tissue_A/CellLine1/CellLine1_training                   │
│      └── Tissue_B/CellLine2/CellLine2_training                   │
│    - cell_lines/SIDM001/ (legacy, optional)                    │
│      └── SIDM001_training                                      │
│    - run_log.txt (execution metadata)                          │
└─────────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────────┐
│ OUTPUT FILES                                                    │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  node_HGNC_dict.csv              → Node mapping reference      │
│  activity_master_matrix.csv      → All data sources combined   │
│  activity_from_master.csv        → Single priority value/cell  │
│  run_log.txt                     → Execution summary           │
│  tissue_folders/                 → Tissue-organized structure   │
│    ├── Tissue_A/                 → e.g., Breast/               │
│    │   ├── CellLine1/             → e.g., MCF7/                 │
│    │   │   ├── CellLine1_training → DrugLogics training file    │
│    │   │   └── ...               → Other files preserved       │
│    │   └── CellLine2/                                             │
│    └── Tissue_B/                                                 │
│        └── ...                                                   │
│                                                                 │
│  cell_lines/ (legacy)           → Flat cell-line structure     │
│    ├── SIDM001/                                                 │
│    ├── SIDM002/                                                 │
│    └── ...                                                      │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

---

## 🎨 Architecture Patterns

### **Configuration Pattern**
- Configuration-driven pipeline (JSON/YAML)
- Defaults merged with user config at runtime
- Validation layer before execution

### **Modular Step Design**
- Node and Activity steps are independent modules
- Each step defines input/output contracts
- Steps can be skipped via configuration

### **Data Flow**
- Files → DataFrames (pandas) → Processing → Output files
- Consistent column naming (SIDM, symbol, source)
- Multi-index structures for cell line × data source

### **Error Handling**
- Early validation of configuration and paths
- Clear error messages for troubleshooting
- Graceful fallbacks for optional parameters

---

## 📝 Documentation Files

| File | Purpose |
|------|---------|
| **README.md** | Main usage guide, API reference, examples |
| **INSTALL.md** | Installation, setup, troubleshooting |
| **QUICKSTART.md** | 5-minute quick reference, common commands |
| **PROJECT_STRUCTURE.md** | This file - code organization (you are here) |
| **notebooks/celios_visualize_bashi_final.ipynb** | Interactive example with tissue filtering & visualization |

---

**Last Updated:** February 2026 | **Maintained by:** CELIOS Development Team
