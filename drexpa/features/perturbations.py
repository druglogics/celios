# Perturbations panel creation module

import os
import logging


class PerturbationPanelBuilder:
    """
    Perturbation panel creation class.

    This class creates perturbation panel files per cell line using combinations_short_df.

    Args:
        combinations_short_df (pandas.DataFrame): DataFrame containing combinations.
        output_directory (str): Base output directory.
        tissue_cline_file (str): Path to tissue-cell line mapping CSV file.
        save (bool): Save output files. Default: False.
    """

    # Class-level constants for column names
    COLUMN_CELL_LINE = "cell_line"
    COLUMN_ANCHOR_ID = "anchor_pipeline_ID"
    COLUMN_LIBRARY_ID = "library_pipeline_ID"
    COLUMN_TISSUE = "TISSUE"
    COLUMN_CELL_LINE_NAME = "CELL_LINE_NAME"

    def __init__(self, combinations_short_df, output_directory=None, tissue_cline_file=None, save: bool = False):
        self.combinations_short_df = combinations_short_df
        self.output_directory = output_directory
        self.tissue_cline_file = tissue_cline_file
        self.save = save

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

        # Set up logging
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)

    def _validate_inputs(self):
        """Validate that all required columns are present in the combinations dataframe."""
        required_cols = [self.COLUMN_CELL_LINE, self.COLUMN_ANCHOR_ID, self.COLUMN_LIBRARY_ID]
        missing_cols = [col for col in required_cols if col not in self.combinations_short_df.columns]
        if missing_cols:
            raise ValueError(f"combinations_short_df is missing required columns: {missing_cols}")

        if self.save and self.output_directory is None:
            raise ValueError("output_directory is required when save=True.")

    def get_perturbation_panel(self, cell_lines=None):
        """
        Create perturbation panel files for specific cell lines.

        Args:
            cell_lines (list or str): List of cell line names or a single cell line name.
                                      If None, process all cell lines.

        Returns:
            dict: A dictionary where keys are cell line names and values are either
                  DataFrames (if save=False) or file paths (if save=True).
        """
        # If no cell lines are specified, process all unique cell lines
        if cell_lines is None:
            cell_lines = self.combinations_short_df[self.COLUMN_CELL_LINE].unique()
        elif isinstance(cell_lines, str):
            cell_lines = [cell_lines]

        results = {}
        for cell_line in cell_lines:
            # self.logger.info(f"Processing cell line: {cell_line}")
            cell_line_perturbation_df = self._get_cell_line_perturbation_df(cell_line)

            if not self.save:
                results[cell_line] = cell_line_perturbation_df
            else:
                cell_line_directory = self._create_cell_line_directory(cell_line)
                perturbation_output_file = self._get_perturbation_output_file(cell_line_directory)
                individual_drugs_list = self._get_individual_drugs_list(cell_line_perturbation_df)
                unique_perturbations_combi = self._get_unique_perturbations(cell_line_perturbation_df)

                self._write_perturbation_panel_file(
                    cell_line,
                    perturbation_output_file,
                    individual_drugs_list,
                    unique_perturbations_combi
                )
                results[cell_line] = perturbation_output_file

        return results

    def _get_cell_line_perturbation_df(self, cell_line):
        """Filter combinations dataframe for a specific cell line."""
        if cell_line not in self.combinations_short_df[self.COLUMN_CELL_LINE].unique():
            raise ValueError(f"Cell line {cell_line} not found in combinations_short_df")
        return self.combinations_short_df[self.combinations_short_df[self.COLUMN_CELL_LINE] == cell_line]

    def _create_cell_line_directory(self, cell_line):
        """Create a hierarchical directory structure: tissue/cell_line/"""
        if self.tissue_mapping and cell_line in self.tissue_mapping:
            tissue = self.tissue_mapping[cell_line]
            cell_line_directory = os.path.join(self.output_directory, tissue, cell_line)
        else:
            # Fallback to original structure if no tissue mapping
            cell_line_directory = os.path.join(self.output_directory, cell_line)
        
        os.makedirs(cell_line_directory, exist_ok=True)
        return cell_line_directory

    def _get_perturbation_output_file(self, cell_line_directory):
        """Get the output file path for the perturbation panel."""
        return os.path.join(cell_line_directory, "perturbations")

    def _get_individual_drugs_list(self, cell_line_perturbation_df):
        """Get a sorted list of individual drugs."""
        individual_drugs_list = list(set(
            cell_line_perturbation_df[self.COLUMN_ANCHOR_ID].unique().tolist() +
            cell_line_perturbation_df[self.COLUMN_LIBRARY_ID].unique().tolist()
        ))
        return sorted(individual_drugs_list)

    def _get_unique_perturbations(self, cell_line_perturbation_df):
        """Get unique perturbation combinations."""
        unique_perturbations_combi = set()
        for _, row in cell_line_perturbation_df.iterrows():
            sorted_ids = tuple(sorted([row[self.COLUMN_ANCHOR_ID], row[self.COLUMN_LIBRARY_ID]]))
            unique_perturbations_combi.add(sorted_ids)
        return unique_perturbations_combi

    def _write_perturbation_panel_file(self,
                                       cell_line,
                                       perturbation_output_file,
                                       individual_drugs_list,
                                       unique_perturbations_combi):
        """Write the perturbation panel file."""
        with open(perturbation_output_file, "w") as f:
            f.write(f"#cell_line_name: {cell_line}\n")
            f.write("#Individual drugs\n")
            for drug in individual_drugs_list:
                f.write(f"{drug}\n")
            f.write("#Combinations\n")
            for combination in unique_perturbations_combi:
                f.write(f"{combination[0]}\t{combination[1]}\n")
