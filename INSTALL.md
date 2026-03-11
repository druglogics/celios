# Installation Guide

## Prerequisites

Before installing CELIOS, ensure you have the following installed on your system:

- **Python 3.8 or higher** (check with `python --version`)
- **pip** (Python package manager, usually comes with Python)
- **Git** (optional, for cloning the repository)
- **Virtual environment tool** (venv or conda, recommended but not required)

## Installation Steps

### Option A: Development Installation (Recommended)

If you plan to modify the code or develop features, install in editable mode:

```bash
# Clone the repository (if not already done)
git clone https://github.com/yourusername/celios.git
cd CELIOS

# Create a virtual environment (optional but recommended)
python -m venv celios_env

# Activate the virtual environment
# On Windows:
celios_env\Scripts\activate
# On macOS/Linux:
source celios_env/bin/activate

# Install in editable mode
pip install -e .
```

### Option B: Standard Installation

For users who only need to run CELIOS without modifying code:

```bash
cd CELIOS
pip install .
```

### Option C: Installation with Development Dependencies

If you want to run tests and use development tools:

```bash
cd CELIOS
pip install -e ".[dev]"
```

**Note:** Development dependencies (like `pytest`) should be added to `setup.py` as `extras_require`.

### Option D: Using Conda

If you prefer using Conda instead of pip:

```bash
# Create a new conda environment
conda create -n celios python=3.10

# Activate the environment
conda activate celios

# Install CELIOS from the repository
cd CELIOS
pip install -e .
```

## Post-Installation Setup

### 1. Verify Installation

Test that CELIOS installed correctly:

```bash
# Check CLI is available
celios --help

# Run the test suite (requires pytest)
pytest tests/ -v
```

### 2. Prepare Your Workspace

Create a directory structure for your CELIOS analysis:

```
my_celios_project/
├── data/
│   ├── activity_input/          # Gene expression, mutations, CNV, TF activities
│   ├── node_dic_input/          # Network SIF file, HGNC symbols, manual overrides
│   └── <your_data>/
├── configs/
│   └── celios_config.yaml        # CELIOS configuration file
└── results/                       # Output directory
```

### 3. Create a Configuration File

Create a YAML or JSON config file to specify your inputs and outputs. See [QUICKSTART.md](QUICKSTART.md) for a configuration template.

### 4. Run Your First Analysis

```bash
celios run --config configs/celios_config.yaml --verbose
```

## Troubleshooting

| Problem | Cause | Solution |
|---------|-------|----------|
| `ModuleNotFoundError: No module named 'celios'` | CELIOS not installed or not on Python path | Install with `pip install -e .` or run from repo root with `python -m celios.cli` |
| `celios: command not found` | Console script not in PATH | Reinstall with `pip install -e .` and ensure virtual environment is activated |
| `ImportError: No module named 'pandas'` | Missing pandas dependency | Run `pip install pandas` or reinstall CELIOS with `pip install -e .` |
| `No such file or directory: config.yaml` | Configuration file path is incorrect | Verify config path is correct (absolute or relative to current directory) |
| `Tests fail after installation` | Testing dependencies missing | Install with `pip install -e ".[dev]"` to include pytest and other dev tools |
| `YAML config not recognized` | PyYAML not installed | Install with `pip install pyyaml` for YAML support (optional; JSON configs work by default) |
| `Permission denied` on Linux/macOS | Virtual environment permissions issue | Ensure you've activated the virtual environment before installing |

## Uninstallation

To remove CELIOS from your environment:

```bash
pip uninstall celios
```

## Upgrading CELIOS

To upgrade to the latest version:

```bash
pip install --upgrade -e .
```

## Next Steps

- **Quick Start:** See [QUICKSTART.md](QUICKSTART.md) for a 5-minute guide to running your first analysis
- **Project Structure:** See [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) to understand the codebase organization
- **Full Documentation:** See [README.md](README.md) for detailed usage, API documentation, and examples
- **Example Notebook:** See `notebooks/1_select_visualize.ipynb` for an interactive walkthrough with visualization

## Getting Help

If you encounter issues not covered in this guide:

1. Check the [Troubleshooting](#troubleshooting) section above
2. Review the example configurations in `tests/test_run_celios.py`
3. Check existing GitHub issues or open a new one with error details
4. Contact the development team at viviamsb@ntnu.no
