"""
DREXPA Pipeline Main Module

This module provides the main pipeline execution logic using centralized configuration.
"""
import logging
import os
import time
import pandas as pd
from .config import Config
from .features.chembl_ids import ChEMBLIDResolver
from .features.doses import DoseProcessor
from .features.targets import TargetProcessor
from .features.node_targets import NodeTargetMapper
from .features.profiles import ProfileBuilder
from .features.panel import PanelMaker
from .features.combinations import CombinationProcessor
from .features.synergies import SynergyProcessor
from .features.perturbations import PerturbationPanelBuilder
from .step_registry import STEP_REGISTRY, resolve_steps
from .resources.database import get_internal_database_path


logger = logging.getLogger(__name__)


RUNTIME_MODE_NO_SYNERGY = "no_synergy_data"
RUNTIME_MODE_WITH_CONCENTRATIONS = "with_concentrations"
RUNTIME_MODE_WITHOUT_CONCENTRATIONS = "without_concentrations"


class DrexpaPipeline:
    """
    Main DREXPA pipeline class that orchestrates all processing steps.
    
    Args:
        config (Config): Configuration object containing all parameters.
        synergy_data_file (str): Path to synergy data file. If None, uses config.
    """
    
    def __init__(self, config: Config, synergy_data_file: str = None):
        self.config = config
        if synergy_data_file:
            self.synergy_data_file = synergy_data_file
        elif 'synergy_data_file' in config.paths:
            self.synergy_data_file = os.path.join(
                config.global_config['base_data_dir'], 
                config.paths['synergy_data_file']
            )
        else:
            self.synergy_data_file = None
        
        # Initialize data storage
        self.synergy_df = None
        self.drug_ids_df = None
        self.drugdoses_df = None
        self.chembl_targets_df = None
        self.node_targets_df = None
        self.drugprofiles_df = None
        self.panel_df = None
        self.combinations_short_df = None
        self.synergy_obs_df = None
        self.step_durations = {}
        self.runtime_mode = RUNTIME_MODE_NO_SYNERGY if self.synergy_data_file is None else None
        self._has_concentration_columns = None
    
    def run_pipeline(self, steps_to_run=None):
        """
        Run the DREXPA pipeline.
        
        Args:
            steps_to_run (str or list): Steps to run. Can be:
                - None: run all steps
                - "until_<step>": run until specified step (inclusive)
                - list: list of specific steps
        """
        # Determine which steps to execute
        steps = self._resolve_steps(steps_to_run)
        
        # Skip load_data step if no synergy data file is provided
        if self.synergy_data_file is None and 'load_data' in steps:
            steps.remove('load_data')
            # Also skip steps that depend on load_data but aren't needed for basic profiles
            steps_to_skip = ['doses', 'combinations', 'perturbations', 'synergies']
            steps = [s for s in steps if s not in steps_to_skip]

        self._preflight_validate(steps)
        self._set_runtime_mode(steps)
        
        print("=" * 80)
        print("RUNNING DREXPA PIPELINE")
        print(f"Runtime mode: {self.runtime_mode}")
        if steps_to_run:
            print(f"Steps to run: {', '.join(steps)}")
        else:
            print("Running all steps")
        print("=" * 80)
        
        # Execute the steps
        for step_name in steps:
            step_info = STEP_REGISTRY[step_name]
            print(f"\n[{step_name.upper()}] {step_info['description']}...")
            logger.info("event=step_start step=%s description=%s", step_name, step_info['description'])
            start_time = time.perf_counter()
            method = getattr(self, step_info['method'])
            method()
            duration = time.perf_counter() - start_time
            self.step_durations[step_name] = duration
            logger.info("event=step_end step=%s duration_seconds=%.3f", step_name, duration)
        
        print("\n" + "=" * 80)
        print("PIPELINE COMPLETED SUCCESSFULLY!")
        print("=" * 80)
        print(f"\nOutput files saved to: {self.config.global_config['output_dir']}")
        self._log_timing_summary()
    
    def _resolve_steps(self, steps_to_run):
        """
        Resolve which steps to run based on the input.
        
        Args:
            steps_to_run (str or list): Steps specification
            
        Returns:
            list: Ordered list of step names to execute
        """
        return resolve_steps(steps_to_run)

    def _required_path_for(self, path_key):
        """Resolve a configured path key relative to base_data_dir."""
        relative_path = self.config.paths.get(path_key)
        if not relative_path:
            return None
        return os.path.join(self.config.global_config.get('base_data_dir', ''), relative_path)

    def _preflight_validate(self, steps):
        """Validate required files and basic tabular schema before execution."""
        self._validate_required_files(steps)
        self._validate_required_columns(steps)

    def _set_runtime_mode(self, steps):
        """Determine explicit runtime mode for this execution."""
        if self.synergy_data_file is None or 'load_data' not in steps:
            self.runtime_mode = RUNTIME_MODE_NO_SYNERGY
        elif self._has_concentration_columns:
            self.runtime_mode = RUNTIME_MODE_WITH_CONCENTRATIONS
        else:
            self.runtime_mode = RUNTIME_MODE_WITHOUT_CONCENTRATIONS

        logger.info("event=runtime_mode mode=%s", self.runtime_mode)

    def _validate_required_files(self, steps):
        """Ensure required files exist for selected steps."""
        required_files = []

        if any(step in steps for step in ['chembl_ids', 'doses', 'targets', 'node_targets', 'profiles', 'panel']):
            required_files.append(('drug_names_file', self._required_path_for('drug_names_file')))

        if any(step in steps for step in ['node_targets', 'profiles', 'panel']):
            required_files.append(('node_dict_file', self._required_path_for('node_dict_file')))

        if any(step in steps for step in ['targets', 'combinations']):
            # DB is internal and managed by DREXPA; validate its availability
            try:
                db_path = get_internal_database_path()
                required_files.append(('db_file (internal)', db_path))
            except FileNotFoundError as e:
                raise FileNotFoundError(
                    f"Internal database validation failed: {e}. "
                    "This suggests a broken DREXPA installation."
                ) from e

        if any(step in steps for step in ['perturbations', 'synergies']):
            required_files.append(('tissue_cline_file', self._required_path_for('tissue_cline_file')))

        if 'load_data' in steps:
            required_files.append(('synergy_data_file', self.synergy_data_file))

        missing = []
        for label, file_path in required_files:
            if not file_path:
                missing.append(f"{label}: not configured")
                continue
            if not os.path.exists(file_path):
                missing.append(f"{label}: {file_path}")

        if missing:
            raise FileNotFoundError(
                "Preflight validation failed. Missing required files:\n- " + "\n- ".join(missing)
            )

    def _validate_required_columns(self, steps):
        """Validate required synergy-data columns when load_data is part of selected steps."""
        if 'load_data' not in steps or self.synergy_data_file is None:
            self._has_concentration_columns = False
            return

        columns_df = pd.read_csv(self.synergy_data_file, nrows=0)
        available_columns = set(columns_df.columns)

        required = set()
        if 'doses' in steps:
            required.add(self.config.columns.get('drug_name_A', 'drug_name_A'))
            if self.config.options.get('double_drug_screen', False):
                required.add(self.config.columns.get('drug_name_B', 'drug_name_B'))

        if 'combinations' in steps:
            required.add(self.config.columns.get('drug_name_A', 'drug_name_A'))
            required.add(self.config.columns.get('drug_name_B', 'drug_name_B'))
            required.add(self.config.columns.get('cell_line', 'cell_line'))
            required.add(self.config.columns.get('synergy', 'synergy'))

        if 'synergies' in steps:
            required.add(self.config.columns.get('cell_line', 'cell_line'))

        missing_columns = sorted([column for column in required if column not in available_columns])
        if missing_columns:
            raise ValueError(
                "Preflight validation failed. Missing required columns in synergy data: "
                + ", ".join(missing_columns)
            )

        concentration_columns = [
            self.config.columns.get('conc_A', 'conc_A'),
            self.config.columns.get('conc_B', 'conc_B')
        ]
        missing_concentrations = [col for col in concentration_columns if col not in available_columns]
        self._has_concentration_columns = len(missing_concentrations) == 0
        if 'doses' in steps and missing_concentrations:
            logger.warning(
                "event=preflight_warning message='Concentration columns missing; dose processing may be skipped' missing_columns=%s",
                ",".join(missing_concentrations)
            )

    def _log_timing_summary(self):
        """Log per-step timing summary after successful run."""
        if not self.step_durations:
            return

        total_seconds = sum(self.step_durations.values())
        logger.info("event=pipeline_summary total_seconds=%.3f", total_seconds)
        for step_name, duration in self.step_durations.items():
            logger.info("event=pipeline_summary_step step=%s duration_seconds=%.3f", step_name, duration)
    
    def _load_synergy_data(self):
        """Load synergy data from file."""
        self.synergy_df = pd.read_csv(self.synergy_data_file)
        print(f"Loaded synergy data: {self.synergy_df.shape}")
    
    def _get_chembl_ids(self):
        """Resolve ChEMBL IDs for drugs."""
        chembl_config = self.config.get_chembl_config()
        resolver = ChEMBLIDResolver(**chembl_config)
        self.drug_ids_df = resolver.resolve_chembl_ids()
        print(f"Got ChEMBL IDs for {len(self.drug_ids_df)} drugs")
    
    def _process_doses(self):
        """Process drug doses from synergy data."""
        concentration_a = self.config.columns.get('conc_A', 'conc_A')

        # Check if synergy data has concentration columns
        if self.synergy_df is None or concentration_a not in self.synergy_df.columns:
            print("No concentration data found in synergy data. Skipping dose processing.")
            self.drugdoses_df = None
            return
            
        doses_config = self.config.get_doses_config()
        doses_config.update({
            'drugscreen_df': self.synergy_df,
            'drugID_file': os.path.join(self.config.global_config['output_dir'], 'drug_ChEMBL_IDs.csv'),
            'column_drugname': self.config.columns.get('drug_name_A', 'drug_name_A'),
            'column_concentration': concentration_a,
            'column_chembl': self.config.columns.get('chembl_id', 'ChEMBL_ID')
        })
        processor = DoseProcessor(**doses_config)
        self.drugdoses_df = processor.get_drugdoses()
        print(f"Processed doses for {len(self.drugdoses_df)} unique drug entries")
    
    def _get_targets(self):
        """Get drug targets from database."""
        targets_config = self.config.get_targets_config()
        
        # For datasets without concentrations, use drug_ids_df instead of drugdoses_df
        if self.drugdoses_df is not None:
            # Use concentration-based processing
            manual_targets_column = [self.config.columns.get('targets_A', 'targets_A'), 
                                   self.config.columns.get('targets_B', 'targets_B')]
            
            # Check if synergy data has the required target columns
            if (self.synergy_df is not None and 
                all(col in self.synergy_df.columns for col in manual_targets_column)):
                manual_targets_df = self.synergy_df
            else:
                manual_targets_df = None
                manual_targets_column = None
            
            targets_config.update({
                'drugdoses_df': self.drugdoses_df,
                'column_drugname': 'drug_name',
                'column_chembl': 'ChEMBL_ID',
                'column_concentration': 'concentration',
                'manual_targets_df': manual_targets_df,
                'manual_targets_column': manual_targets_column,
                'merge_on': 'ChEMBL_conc'
            })
        else:
            # Use ChEMBL ID-based processing without concentrations
            # Check if synergy data has target columns for manual merging
            double_drug = self.config.options.get('double_drug_screen', False)
            manual_targets_column = ['targets_A', 'targets_B'] if double_drug else 'targets'
            
            # Check if synergy data exists and has the required target columns
            if (self.synergy_df is not None and 
                all(col in self.synergy_df.columns for col in 
                    (manual_targets_column if isinstance(manual_targets_column, list) else [manual_targets_column]))):
                manual_targets_df = self.synergy_df
            else:
                manual_targets_df = None
                manual_targets_column = None
            
            targets_config.update({
                'drugdoses_df': self.drug_ids_df,
                'column_drugname': 'drug_name',
                'column_chembl': 'ChEMBL_ID',
                'ic50_value': 10000,  # Default IC50 threshold (10 µM)
                'merge_on': 'ChEMBL_ID',  # Merge on ID only, no concentration
                'manual_targets_df': manual_targets_df,
                'manual_targets_column': manual_targets_column
            })
        
        processor = TargetProcessor(**targets_config)
        self.chembl_targets_df = processor.get_chembl_targets()
        print("Got targets for drugs")
    
    def _map_node_targets(self):
        """Map targets to consensus nodes."""
        node_config = self.config.get_node_targets_config()
        node_config.update({
            'chembl_targets_df': self.chembl_targets_df,
            'column_drugname': 'drug_name',
            'column_target': 'targets'
        })
        mapper = NodeTargetMapper(**node_config)
        self.node_targets_df = mapper.get_node_targets()
        print("Mapped targets to nodes")
    
    def _build_profiles(self):
        """Build drug profiles."""
        profiles_config = self.config.get_profiles_config()
        profiles_config.update({
            'node_targets_df': self.node_targets_df
        })
        builder = ProfileBuilder(**profiles_config)
        self.drugprofiles_df = builder.get_drugprofiles()
        print(f"Built drug profiles: {self.drugprofiles_df.shape}")
        
        # For combination screens without concentrations, create combinations_short_df now that profiles are available
        if (self.combinations_short_df is None and 
            self.synergy_df is not None and 
            self.config.options.get('double_drug_screen', False)):
            
            import pandas as pd
            print("Creating combinations data from synergy data (no concentrations)...")
            
            # Get list of drugs that have profiles
            drugs_with_profiles = set(self.drugprofiles_df['drug_name'].unique())
            
            # Filter synergy_df to only combinations where both drugs have profiles
            filtered_synergy_df = self.synergy_df[
                self.synergy_df['drug_name_A'].isin(drugs_with_profiles) & 
                self.synergy_df['drug_name_B'].isin(drugs_with_profiles)
            ].copy()
            
            if len(filtered_synergy_df) == 0:
                print("Warning: No combinations found where both drugs have profiles. Skipping combinations creation.")
                self.combinations_short_df = pd.DataFrame()
                self.synergy_obs_df = pd.DataFrame()
                return
            
            # Create mapping from drug names to PD_IDs
            # If a drug has multiple PD_profiles (e.g., different concentrations),
            # choose the first PD_profile to create a deterministic mapping for
            # no-concentration combination screens.
            drug_profiles = (
                self.drugprofiles_df.groupby('drug_name')['PD_profile']
                .first()
                .to_dict()
            )

            # Create combinations_short_df
            combinations_short_df = filtered_synergy_df.copy()
            combinations_short_df['anchor_pipeline_ID'] = combinations_short_df['drug_name_A'].map(drug_profiles)
            combinations_short_df['library_pipeline_ID'] = combinations_short_df['drug_name_B'].map(drug_profiles)
            
            # Rename columns to match expected format
            combinations_short_df = combinations_short_df.rename(columns={
                'tissue': 'TISSUE'
            })
            
            self.combinations_short_df = combinations_short_df
            self.synergy_obs_df = combinations_short_df.copy()
            print(f"Created combinations data: {self.combinations_short_df.shape}")
    
    def _create_panel(self):
        """Create drug panel."""
        panel_config = self.config.get_panel_config()
        panel_config.update({
            'drugprofiles_df': self.drugprofiles_df
        })
        maker = PanelMaker(**panel_config)
        self.panel_df = maker.get_drugpanel()
        print(f"Created drug panel: {self.panel_df.shape}")
    
    def _prepare_combinations(self):
        """Prepare combinations dataframe."""
        combos_config = self.config.get_combinations_config()
        combos_config.update({
            'chembl_targets_df': self.chembl_targets_df,
            'drugprofiles_df': self.drugprofiles_df,
            'drugscreen_df': self.synergy_df
        })
        processor = CombinationProcessor(**combos_config)
        self.combinations_short_df, self.synergy_obs_df = processor.prepare_combinations()
        print(f"Prepared combinations: {self.combinations_short_df.shape}")
    
    def _create_perturbations(self):
        """Create perturbation panels."""
        if self.combinations_short_df is None or self.combinations_short_df.empty:
            print("No combination data available. Skipping perturbation panel creation.")
            return
            
        pert_config = self.config.get_perturbations_config()
        pert_config.update({
            'combinations_short_df': self.combinations_short_df
        })
        builder = PerturbationPanelBuilder(**pert_config)
        builder.get_perturbation_panel()
        print("Created perturbation panels for all cell lines")
    
    def _process_synergies(self):
        """Process synergies for all cell lines."""
        if self.synergy_obs_df is None or self.synergy_obs_df.empty:
            print("No synergy observation data available. Skipping synergy processing.")
            return
            
        synergy_config = self.config.get_synergies_config()
        synergy_config.update({
            'synergies_obs_df': self.synergy_obs_df
        })
        processor = SynergyProcessor(**synergy_config)
        
        cell_lines = self.synergy_obs_df['cell_line'].unique()
        for cell_line in cell_lines:
            # print(f"Processing synergies for cell line: {cell_line}")
            processor.observed_synergies(cell_lines=cell_line)


def run_pipeline(config_dict=None, synergy_data_file=None, steps_to_run=None):
    """
    Convenience function to run the DREXPA pipeline.
    
    Args:
        config_dict (dict): Custom configuration dictionary. Uses default if None.
        synergy_data_file (str): Path to synergy data file.
        steps_to_run (str or list): Steps to run. Can be:
            - None: run all steps
            - "until_<step>": run until specified step (e.g., "until_profiles")
            - list: list of specific steps (e.g., ["profiles", "panel"])
    """
    config = Config(config_dict)
    pipeline = DrexpaPipeline(config, synergy_data_file)
    pipeline.run_pipeline(steps_to_run)