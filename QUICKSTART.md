# ⚡ QUICKSTART

Get CELIOS running in **5 minutes**.

---

## 🚀 Get Started in 3 Steps

### **Step 1: Install**
```bash
git clone https://github.com/yourusername/celios.git
cd CELIOS
pip install -e .
```

### **Step 2: Create a Config File**
Save as `celios_config.yaml`:
```yaml
paths:
  base: "."
  input: "data"
  output: "results"
  cellfiles_dir: "results/cell_lines"
  tissue_dir: "results/tissue_folders"  # Optional: organize DL files by tissue/cell_line

steps:
  Node:
    node_input: "data/node_dic_input/cell_fate_plus.sif"
    hgnc_symbols_file: "data/node_dic_input/hgnc_complete_set.txt"
    manual_symbols_file: "data/node_dic_input/manual_symbols.csv"
    include_alias_previous_symbols: false
    directory_output: "results"
  
  Activity:
    activity_file: "data/activity_input/rnaseq_tpm_20220624.csv"
    cell_line_file: "data/activity_input/cell_line_list.csv"
    tf_activity_file: "data/activity_input/ccle_tf_activities.csv"
    mutations_file: "data/activity_input/CCLE_muts_binary.csv"
    cnv_file: "data/activity_input/CCLE_CNV_binary.csv"
    directory_output: "results"
    data_sources: ["mutations", "cnv", "TF"]
```

**Note:** For tissue-organized output, `cell_line_list.csv` must contain `tissue`, `SIDM`, and `cell_line_name` columns. For legacy mode, no specific columns are required.

**Format Overrides (optional):** CELIOS auto-detects file formats. To override:
```yaml
Activity:
  format_override: "26q1"  # Force activity file format
  mutations_format_override: "26q1"  # Force mutations format (old | 26q1)
  cnv_format_override: "old"  # Force CNV format (old | 26q1)
```

### **Step 3: Run**
```bash
celios run --config celios_config.yaml --verbose
```

✅ Done! Check `results/` for outputs.

---

## 📚 Documentation Map

| Need | File | Purpose |
|------|------|---------|
| **Installation help** | [INSTALL.md](INSTALL.md) | Detailed setup, troubleshooting, virtual environments |
| **Project organization** | [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) | Code architecture, module descriptions, data flow |
| **Full documentation** | [README.md](README.md) | Usage examples, API reference, input/output formats |
| **Interactive tutorial** | `notebooks/1_select_visualize.ipynb` | Step-by-step notebook with visualization |

---

## ⚡ Common Commands

| Task | Command |
|------|---------|
| **Run full pipeline** | `celios run --config config.yaml` |
| **Run with verbose output** | `celios run --config config.yaml --verbose` |
| **See execution plan (no run)** | `celios run --config config.yaml --plan` |
| **Stop after Node step** | `celios run --config config.yaml --stop-after node` |
| **Stop after Activity step** | `celios run --config config.yaml --stop-after activity` |
| **Get help** | `celios --help` |
| **Generate nodes from SIF** | `celios node-from-sif --sif network.sif --hgnc hgnc.txt --out nodes.csv` |
| **Generate nodes from list** | `celios node-from-object --input "TP53,BRCA1" --hgnc hgnc.txt --out nodes.csv` |

---

## 🐍 Python API Usage

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
        "Node": {...},  # See config example above
        "Activity": {...},
    },
}

# Run the pipeline
artifacts = run_celios(config=config, verbose=True)

# Access outputs
activity_matrix = artifacts['activity_matrix']
print(activity_matrix.head())
```

---

## 🧪 Testing

```bash
# Run all tests
pytest tests/ -v

# Run specific test file
pytest tests/test_run_celios.py -v

# Run with coverage (if pytest-cov installed)
pytest tests/ --cov=celios
```

---

## 🔧 Configuration Options

| Option | Type | Required | Default | Description |
|--------|------|----------|---------|-------------|
| `paths.base` | str | ✓ | — | Base directory for relative paths |
| `paths.input` | str | ✓ | — | Input data directory |
| `paths.output` | str | ✓ | — | Output directory |
| `paths.cellfiles_dir` | str | ✗ | — | Directory for per-cell-line training files (legacy mode) |
| `paths.tissue_dir` | str | ✗ | — | Root directory for tissue-organized cell line folders |
| `paths.tissue_dir` | str | ✗ | — | Root directory for tissue-organized cell line folders |
| `steps.Node` | dict | ✗ | — | Network analysis step config (skip if using pre-built node_dic) |
| `steps.Activity` | dict | ✓ | — | Activity extraction step config |
| `steps.Activity.data_sources` | list | ✗ | ["mutations","cnv","TF","expression"] | Priority order for activity data |

---

## 🚨 Quick Troubleshooting

| Issue | Quick Fix |
|-------|-----------|
| `celios: command not found` | Reinstall: `pip install -e .` |
| `ModuleNotFoundError: No module named 'celios'` | Run from repo root: `python -m celios.cli run --config config.yaml` |
| `FileNotFoundError: config.yaml` | Use absolute path or verify file exists: `ls config.yaml` |
| `YAML config not recognized` | Install PyYAML: `pip install pyyaml` (JSON configs work by default) |
| `Cannot auto-detect binary matrix format` | Check file headers or use `mutations_format_override` or `cnv_format_override` |
| `ValueError: No ModelID values found in Model.csv` | Verify ModelID values exist in `src/celios/features/Model.csv` |

See [INSTALL.md](INSTALL.md) for more troubleshooting.

---

## 📊 What You Get

✅ **Node Dictionary** (`node_HGNC_dict.csv`)
- Network nodes mapped to HGNC gene symbols
- Enables consistent gene annotation

✅ **Activity Master Matrix** (`activity_master_matrix.csv`)
- All data sources combined (mutations, CNV, TF, expression)
- Format: nodes × cell_lines__data_source

✅ **Priority-Filtered Activity** (`activity_from_master.csv`)
- One activity value per node-cell line (highest priority source)
- Ready for downstream analysis

✅ **Per-Cell-Line Training Files** (optional)
- DrugLogics-compatible format
- Saved to `cellfiles_dir/` if specified

---

## 🎯 Next Steps

1. ✅ **Installed?** → Run your first analysis with the config above
2. 📖 **Need details?** → Read [README.md](README.md) for API reference and examples
3. 🔍 **Explore code?** → See [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md)
4. 📓 **Want visuals?** → Check `notebooks/1_select_visualize.ipynb`

---

**Version:** 0.0.1 | **Status:** Active | **License:** See [LICENSE](LICENSE)
