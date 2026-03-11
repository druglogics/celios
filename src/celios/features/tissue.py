"""Tissue-aware DrugLogic training file writing.

This module provides functionality to organize and write DrugLogic pipeline
training files by tissue and cell line, using metadata from a cell line file.
"""

from pathlib import Path
import pandas as pd
from typing import Optional
import os

from . import files as files_mod


def write_tissue_files(
    activity_df: pd.DataFrame,
    tissue_dir: str,
    cell_line_file: str,
    proliferation_state: bool = True,
    dna_damage_state: bool = False,
    verbose: bool = True
) -> str:
    """Write DrugLogic training files organized by tissue/cell_line.

    Args:
        activity_df: DataFrame with 'symbol' index and cell line columns.
        tissue_dir: Root directory for tissue-organized output.
        cell_line_file: Path to CSV with 'tissue', 'SIDM', 'cell_line_name' columns.
        proliferation_state: Include proliferation state in training files.
        dna_damage_state: Include DNA damage state in training files.
        verbose: Enable logging.

    Returns:
        Root tissue directory path.

    Raises:
        ValueError: If required columns are missing in cell_line_file.
    """
    # Load CSV with case-insensitive column matching
    cell_df = pd.read_csv(cell_line_file)
    tissue_col = None
    cell_name_col = None
    for col in cell_df.columns:
        if col.lower() == 'tissue':
            tissue_col = col
        elif col.lower() == 'cell_line_name':
            cell_name_col = col

    if not tissue_col or not cell_name_col:
        raise ValueError(f"cell_line_file must contain 'tissue' and 'cell_line_name' columns (case-insensitive). Found: {list(cell_df.columns)}")

    # Build mappings: SIDM -> tissue and SIDM -> cell_line_name
    sidm_col = None
    cell_name_col = None
    for col in cell_df.columns:
        if col.lower() == 'sidm':
            sidm_col = col
        elif col.lower() == 'cell_line_name':
            cell_name_col = col
    if not sidm_col or not cell_name_col:
        raise ValueError(f"For tissue-organized output, cell_line_file must contain 'SIDM' and 'cell_line_name' columns (case-insensitive). Found: {list(cell_df.columns)}")

    tissue_mapping = dict(zip(cell_df[sidm_col], cell_df[tissue_col]))
    name_mapping = dict(zip(cell_df[sidm_col], cell_df[cell_name_col]))
    
    # Create reverse mapping from cell_line_name to SIDM
    name_to_sidm = dict(zip(cell_df[cell_name_col], cell_df[sidm_col]))

    # Create tissue-organized directories and write files
    tissue_root = Path(tissue_dir).expanduser().resolve()
    tissue_root.mkdir(parents=True, exist_ok=True)

    cell_dirs_created = []
    for cell_line in activity_df.columns:
        if cell_line == 'symbol':  # Skip the symbol column
            continue
        if cell_line not in name_to_sidm:
            if verbose:
                print(f"Warning: Cell line '{cell_line}' not found in cell_line_list.csv, skipping DL files")
            continue

        sidm = name_to_sidm[cell_line]
        tissue = tissue_mapping[sidm]
        cell_name = name_mapping[sidm]
        # Apply the same cleaning as DL_trainingfiles
        cell_name = ''.join(e for e in cell_name if e.isalnum()).upper()
        cell_dir = tissue_root / tissue / cell_name
        cell_dir.mkdir(parents=True, exist_ok=True)
        cell_dirs_created.append(str(cell_dir))

        # Write DL training files for this cell line
        try:
            cell_activity = pd.DataFrame({
                'symbol': activity_df.index,
                cell_line: activity_df[cell_line]
            })
            files_mod.DL_trainingfiles(cell_activity, str(cell_dir.parent),
                                     proliferation_state=proliferation_state,
                                     dna_damage_state=dna_damage_state)
        except Exception as e:
            if verbose:
                print(f"Error writing DL files for {cell_line}: {e}")

    if verbose:
        print(f"Wrote training files to {len(cell_dirs_created)} cell-line directories across {len(set(tissue_mapping.values()))} tissue folders")

    return str(tissue_root)