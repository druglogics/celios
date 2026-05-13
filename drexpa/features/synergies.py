# Synergies processing module
import os

class SynergyProcessor:
    """
    Synergies processing class.

    This class creates observed synergies files per cell line, handling both
    `synergies_obs_df` (with PD_IDs) and `synergy_df` (with drug names).

    Args:
        synergies_obs_df (pandas.DataFrame): DataFrame containing synergies with PD_IDs.
        synergy_df (pandas.DataFrame): DataFrame containing synergies with drug names.
        output_directory (str): Base output directory.
        tissue_cline_file (str): Path to tissue-cell line mapping CSV file.
        save (bool): Save output files. Default: False.
        threshold (float): Threshold for synergy score to consider as synergistic. Default: 0.0.
    """

    # Class-level constants for column names
    COLUMN_CELL_LINE = 'cell_line'
    COLUMN_SYNERGY = 'synergy'
    COLUMN_ANCHOR_ID = 'anchor_pipeline_ID'
    COLUMN_LIBRARY_ID = 'library_pipeline_ID'
    COLUMN_ANCHOR_NAME = 'drug_name_A'
    COLUMN_LIBRARY_NAME = 'drug_name_B'
    COLUMN_TISSUE = "TISSUE"
    COLUMN_CELL_LINE_NAME = "CELL_LINE_NAME"

    def __init__(self, synergies_obs_df=None, synergy_df=None, output_directory=None, tissue_cline_file=None, save: bool = False, threshold: float = 0.0):
        self.synergies_obs_df = synergies_obs_df
        self.synergy_df = synergy_df
        self.output_directory = output_directory
        self.tissue_cline_file = tissue_cline_file
        self.save = save
        self.threshold = threshold

        # Load tissue-cell line mapping if provided
        self.tissue_mapping = None
        if self.tissue_cline_file and os.path.exists(self.tissue_cline_file):
            import pandas as pd
            self.tissue_mapping = pd.read_csv(self.tissue_cline_file)
            self.tissue_mapping = dict(zip(
                self.tissue_mapping[self.COLUMN_CELL_LINE_NAME], 
                self.tissue_mapping[self.COLUMN_TISSUE]
            ))

        # Validate inputs
        self._validate_inputs()

    def _validate_inputs(self):
        """Validate that at least one input DataFrame is provided and contains required columns."""
        if self.synergies_obs_df is None and self.synergy_df is None:
            raise ValueError("At least one of `synergies_obs_df` or `synergy_df` must be provided.")

        if self.synergies_obs_df is not None:
            required_cols_obs = [
                self.COLUMN_CELL_LINE,
                self.COLUMN_SYNERGY,
                self.COLUMN_ANCHOR_ID,
                self.COLUMN_LIBRARY_ID
            ]
            missing_cols_obs = [col for col in required_cols_obs if col not in self.synergies_obs_df.columns]
            if missing_cols_obs:
                raise ValueError(f"synergies_obs_df is missing required columns: {missing_cols_obs}")
            if self.synergies_obs_df.empty:
                raise ValueError("synergies_obs_df is empty.")

        if self.synergy_df is not None:
            required_cols_synergy = [
                self.COLUMN_CELL_LINE,
                self.COLUMN_SYNERGY,
                self.COLUMN_ANCHOR_NAME,
                self.COLUMN_LIBRARY_NAME
            ]
            missing_cols_synergy = [col for col in required_cols_synergy if col not in self.synergy_df.columns]
            if missing_cols_synergy:
                raise ValueError(f"synergy_df is missing required columns: {missing_cols_synergy}")
            if self.synergy_df.empty:
                raise ValueError("synergy_df is empty.")

        if self.save and self.output_directory is None:
            raise ValueError("output_directory is required when save=True.")

    def observed_synergies(self, cell_lines=None):
        """
        Create observed synergies output for specific cell lines.

        Args:
            cell_lines (list or str): List of cell line names or a single cell line name.
                                      If None, process all cell lines.

        Returns:
            dict: A dictionary where keys are cell line names and values are either
                  DataFrames (if save=False) or file paths (if save=True).
        """
        # If no cell lines are specified, process all unique cell lines
        if cell_lines is None:
            cell_lines = set()
            if self.synergies_obs_df is not None:
                cell_lines.update(self.synergies_obs_df[self.COLUMN_CELL_LINE].unique())
            if self.synergy_df is not None:
                cell_lines.update(self.synergy_df[self.COLUMN_CELL_LINE].unique())
            cell_lines = list(cell_lines)
        elif isinstance(cell_lines, str):
            cell_lines = [cell_lines]

        results = {}
        for cell_line in cell_lines:
            # Process synergies_obs_df for PD_IDs
            if self.synergies_obs_df is not None:
                cell_line_synergies_obs = self._filter_synergies(self.synergies_obs_df, cell_line)
                if self.save:
                    cell_line_dir = self._create_cell_line_directory(cell_line)
                    output_file_obs = self._get_output_file(cell_line_dir, "observed_synergies")
                    self._write_synergies_with_PD_IDs(output_file_obs, cell_line_synergies_obs)
                    results[f"{cell_line}_PD_IDs"] = output_file_obs
                else:
                    results[f"{cell_line}_PD_IDs"] = cell_line_synergies_obs

            # Process synergy_df for drug names
            if self.synergy_df is not None:
                cell_line_synergies = self._filter_synergies(self.synergy_df, cell_line)
                if self.save:
                    cell_line_dir = self._create_cell_line_directory(cell_line)
                    output_file_synergy = self._get_output_file(cell_line_dir, "observed_synergies_with_drug_names")
                    self._write_synergies_with_drug_names(output_file_synergy, cell_line_synergies)
                    results[f"{cell_line}_drug_names"] = output_file_synergy
                else:
                    results[f"{cell_line}_drug_names"] = cell_line_synergies

        return results

    def _filter_synergies(self, df, cell_line):
        """Filter synergies dataframe for a specific cell line."""
        return df[(df[self.COLUMN_CELL_LINE] == cell_line) & (df[self.COLUMN_SYNERGY] > self.threshold)]

    def _create_cell_line_directory(self, cell_line):
        """Create a hierarchical directory structure: tissue/cell_line/"""
        if self.tissue_mapping and cell_line in self.tissue_mapping:
            tissue = self.tissue_mapping[cell_line]
            cell_line_dir = os.path.join(self.output_directory, tissue, cell_line)
        else:
            # Fallback to original structure if no tissue mapping
            cell_line_dir = os.path.join(self.output_directory, cell_line)
        
        os.makedirs(cell_line_dir, exist_ok=True)
        return cell_line_dir

    def _get_output_file(self, cell_line_dir, file_name):
        """Get the output file path for the synergies."""
        return os.path.join(cell_line_dir, file_name)

    def _write_synergies_with_PD_IDs(self, output_file, cell_line_synergies):
        """Write synergies with PD_IDs to the output file."""
        unique_combinations = set()
        try:
            with open(output_file, 'w') as f:
                for _, row in cell_line_synergies.iterrows():
                    anchor_id = row[self.COLUMN_ANCHOR_ID]
                    library_id = row[self.COLUMN_LIBRARY_ID]
                    sorted_ids = tuple(sorted([anchor_id, library_id]))
                    if sorted_ids not in unique_combinations:
                        unique_combinations.add(sorted_ids)
                        f.write(f"{sorted_ids[0]}~{sorted_ids[1]}\n")
        except Exception as e:
            raise RuntimeError(f"Error occurred while writing synergies with PD_IDs: {e}")

    def _write_synergies_with_drug_names(self, output_file, cell_line_synergies):
        """Write synergies with drug names to the output file."""
        try:
            with open(output_file, 'w') as f:
                for _, row in cell_line_synergies.iterrows():
                    anchor_name = row[self.COLUMN_ANCHOR_NAME]
                    library_name = row[self.COLUMN_LIBRARY_NAME]
                    f.write(f"{anchor_name}\t{library_name}\n")
        except Exception as e:
            raise RuntimeError(f"Error occurred while writing synergies with drug names: {e}")
