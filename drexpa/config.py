"""
DREXPA Pipeline Configuration

This module contains the configuration for the DREXPA drug synergy analysis pipeline.

For released module usage:
1. Copy get_default_config() and modify the paths to point to your data files
2. Or create a custom config JSON file and use --config option
3. Or modify the default config directly for your use case

Required data files:
- drug_names_file: Text file with one identifier per line, or CSV with identifier columns
- synergy_data_file: CSV file with drug combination screening data
- node_dict_file: CSV file mapping gene symbols to model nodes
- db_file: SQLite database with drug-target interactions

Example config JSON:
{
  "global": {"output_dir": "results", "base_data_dir": "my_data"},
  "paths": {
    "drug_names_file": "drugs.txt",
    "synergy_data_file": "screen_data.csv",
    "node_dict_file": "gene_nodes.csv"
  }
}
"""
from copy import deepcopy
import os

from .resources.database import get_internal_database_path


def get_default_config():
    """
    Return a default configuration dictionary for the DREXPA pipeline.
    Adjust values to match your dataset and file paths.
    
    For released module, users should modify these paths to point to their data files.
    
    Note: The drug-target interaction database (db_file) is managed internally by DREXPA
    and is shipped with the package. It should not be changed by users under normal circumstances.
    For testing or advanced scenarios, it can be overridden via custom config JSON.
    """
    return {
        "global": {
            "output_dir": "output",
            "verbose": False,
            "save": True,
            "base_data_dir": "data"
        },
        "paths": {
            # Modify these paths to point to your specific data files
            "drug_names_file": "drug_names.txt",  # Path to drug names file (one drug per line)
            "synergy_data_file": "synergy_data.csv",  # Path to synergy screening data (CSV format)
            "node_dict_file": "node_dict.csv",  # Path to node dictionary (gene-to-node mapping)
            "tissue_cline_file": "tissue_cline.csv",  # Path to tissue-cell line mapping (CSV format)
            "db_file": None,  # Auto-resolved to internal database; override only for advanced use
            "manual_chembl_csv": "manual_chembl.csv"  # Path to manual ChEMBL IDs CSV (optional)
        },
        "columns": {
            "drug_name": "drug_name",
            "chembl_id": "ChEMBL_ID",
            "concentration": "concentration",
            "cell_line": "cell_line",
            "synergy": "synergy",
            "drug_name_A": "drug_name_A",
            "drug_name_B": "drug_name_B",
            "conc_A": "conc_A",
            "conc_B": "conc_B",
            "targets_A": "targets_A",
            "targets_B": "targets_B",
            "targets": "targets"
        },
        "options": {
            "synergy_threshold": 0.0,
            "ic50_value": None,
            "double_drug_screen": True,
            "original_target_merge": "fill_missing",
            "profile_values": "drug_name"
        }
    }


def merge_config(base_config, override_config):
    """
    Deep-merge configuration dictionaries.

    - Nested dictionaries are merged recursively.
    - Scalar values and non-dict objects in override_config replace base values.
    - Inputs are not mutated.
    """
    if override_config is None:
        return deepcopy(base_config)

    if not isinstance(base_config, dict) or not isinstance(override_config, dict):
        return deepcopy(override_config)

    merged = deepcopy(base_config)
    for key, value in override_config.items():
        if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
            merged[key] = merge_config(merged[key], value)
        else:
            merged[key] = deepcopy(value)
    return merged


def validate_config(config):
    """
    Validate configuration dictionary structure.
    
    Args:
        config (dict): Configuration dictionary to validate
        
    Raises:
        ValueError: If configuration is invalid
    """
    if not isinstance(config, dict):
        raise ValueError("Configuration must be a dictionary")
    
    required_sections = ["global", "paths"]
    for section in required_sections:
        if section not in config:
            raise ValueError(f"Missing required configuration section: {section}")
    
    # Check global section has required keys
    required_global = ["output_dir", "verbose", "save"]
    for key in required_global:
        if key not in config["global"]:
            raise ValueError(f"Missing required global config key: {key}")
    
    # Check paths section has required keys
    required_paths = ["drug_names_file", "node_dict_file"]
    for key in required_paths:
        if key not in config["paths"]:
            raise ValueError(f"Missing required paths config key: {key}")


class Config:
    """
    Configuration class for DREXPA pipeline.
    Provides methods to get parameter dictionaries for each module.
    """
    def __init__(self, config_dict=None):
        self.config = config_dict or get_default_config()
        self.global_config = self.config.get('global', {})
        self.paths = self.config.get('paths', {})
        self.columns = self.config.get('columns', {})
        self.options = self.config.get('options', {})

    def get_chembl_config(self):
        """Get config for ChEMBLIDResolver."""
        drug_names_file_path = os.path.join(self.global_config.get('base_data_dir', ''), self.paths.get('drug_names_file', ''))
        manual_chembl_file_path = os.path.join(self.global_config.get('base_data_dir', ''), self.paths.get('manual_chembl_csv', ''))
        
        return {
            'directory_output': self.global_config.get('output_dir', ''),
            'verbose': self.global_config.get('verbose', False),
            'save': self.global_config.get('save', False),
            'drugnames_file': drug_names_file_path,
            'column_drugname': self.columns.get('drug_name', 'drug_name'),
            'column_chembl': self.columns.get('chembl_id', 'ChEMBL_ID'),
            'column_pubchem': 'PubChemID',
            'manual_chembl_file': manual_chembl_file_path if os.path.exists(manual_chembl_file_path) else None
        }

    def get_doses_config(self):
        """Get config for DoseProcessor."""
        return {
            'drugnames_file': os.path.join(self.global_config.get('base_data_dir', ''), self.paths.get('drug_names_file', '')),
            'dual_drug_mode': self.options.get('double_drug_screen', False),
            'column_drugname_b': self.columns.get('drug_name_B'),
            'column_concentration_b': self.columns.get('conc_B'),
            'directory_output': self.global_config.get('output_dir', ''),
            'verbose': self.global_config.get('verbose', False),
            'save': self.global_config.get('save', False)
        }

    def get_targets_config(self):
        """Get config for TargetProcessor."""
        db_file = self.paths.get('db_file')
        if db_file is None:
            db_file = get_internal_database_path()
        else:
            db_file = os.path.join(self.global_config.get('base_data_dir', ''), db_file)
        
        return {
            'db_file': db_file,
            'ic50_value': self.options.get('ic50_value'),
            'merge_strategy': self.options.get('original_target_merge', 'fill_missing'),
            'directory_output': self.global_config.get('output_dir', ''),
            'verbose': self.global_config.get('verbose', False),
            'save': self.global_config.get('save', False)
        }

    def get_profiles_config(self):
        """Get config for ProfileBuilder."""
        return {
            'profile_values': self.options.get('profile_values', 'drug_name'),
            'directory_output': self.global_config.get('output_dir', ''),
            'verbose': self.global_config.get('verbose', False),
            'save': self.global_config.get('save', False)
        }

    def get_combinations_config(self):
        """Get config for CombinationProcessor."""
        return {
            'column_cell_line_name': self.columns.get('cell_line', 'cell_line'),
            'column_synergy': self.columns.get('synergy', 'synergy')
        }

    def get_node_targets_config(self):
        """Get config for NodeTargetMapper."""
        return {
            'node_dict_file': os.path.join(self.global_config.get('base_data_dir', ''), self.paths.get('node_dict_file', '')),
            'column_drugname': self.columns.get('drug_name', 'drug_name'),
            'column_target': self.columns.get('targets', 'targets'),
            'directory_output': self.global_config.get('output_dir', ''),
            'verbose': self.global_config.get('verbose', False),
            'save': self.global_config.get('save', False)
        }

    def get_panel_config(self):
        """Get config for PanelMaker."""
        return {
            'directory_main_output': self.global_config.get('output_dir', ''),
            'directory_suppl_output': self.global_config.get('output_dir', ''),
            'save': self.global_config.get('save', False)
        }

    def get_perturbations_config(self):
        """Get config for PerturbationPanelBuilder."""
        return {
            'output_directory': self.global_config.get('output_dir', ''),
            'tissue_cline_file': os.path.join(self.global_config.get('base_data_dir', ''), self.paths.get('tissue_cline_file', '')),
            'save': self.global_config.get('save', False)
        }

    def get_synergies_config(self):
        """Get config for SynergyProcessor."""
        return {
            'output_directory': self.global_config.get('output_dir', ''),
            'tissue_cline_file': os.path.join(self.global_config.get('base_data_dir', ''), self.paths.get('tissue_cline_file', '')),
            'save': self.global_config.get('save', False),
            'threshold': self.options.get('synergy_threshold', 0.0)
        }
