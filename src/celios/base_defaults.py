# celios/pipeline_steps_arg.py

# PIPELINE STEPS and ARGUMENTS constants for Celios

"""Pipeline configuration template.

This module exposes a single `CONFIG_TEMPLATE` dict intended as the
authoritative example/default configuration for the pipeline. It replaces the
older `BASE_DEFAULTS` concept: callers should merge user-provided `steps`
into `CONFIG_TEMPLATE['steps']` (see `core._build_config`).
"""

CONFIG_TEMPLATE = {
    # Top-level path settings used to resolve relative filenames
    "paths": {
        "base": ".",
        "input": "data/activity_input",
        "output": None,
        "cellfiles_dir": "results/cell_lines",
        "tissue_dir": None,  # Root directory for tissue-organized cell line folders
    },

    # General run-level options
    "general": {
        "verbose": False,
        "log_level": "INFO",
        "force": False,
        "seed": 42,
    },

    # Per-step defaults. Tests and code expect keys like 'Node' and 'Activity'.
    "steps": {
        "run": {
            "work_dir": None,
            "log_level": "INFO",
            "force": False,
            "seed": 42,
        },

        "Node": {
            "node_input": None,
            "hgnc_symbols_file": None,
            "directory_output": None,
            "manual_symbols_file": None,
            "include_alias_previous_symbols": False,
            "verbose": False,
        },

        "Activity": {
            "activity_file": None,
            "cell_line_file": None,
            "format_override": None,
            "node_dict": None,
            "tf_activity_file": None,
            "mutations_file": None,
            "mutations_format_override": None,
            "cnv_file": None,
            "cnv_format_override": None,
            "directory_output": None,
            "verbose": False,
            "data_sources": ["mutations", "cnv", "TF", "expression"],
        },
    },
}