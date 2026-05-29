import json
import os
from importlib import resources
from pathlib import Path

import pandas as pd

from . import report as report_mod
from .cell_line_resolver import (
	resolve_identifiers_to_sidm,
	resolve_model_ids_to_sidm,
	resolve_sidm_from_dataframe,
)


def save_file(data, directory_path, file_name, file_type='csv', **kwargs):
	"""
	Save data to a file in the specified format.
	"""
	if directory_path:
		os.makedirs(directory_path, exist_ok=True)

	full_path = os.path.join(directory_path, file_name) if directory_path else file_name

	if file_type == 'csv':
		if isinstance(data, pd.DataFrame):
			data.to_csv(full_path, **kwargs)
			report_mod.add_log(f'Data saved to {full_path} as CSV.')
		else:
			raise ValueError('For CSV, data should be a pandas DataFrame.')
	elif file_type == 'json':
		if isinstance(data, dict):
			with open(full_path, 'w') as file:
				json.dump(data, file)
			report_mod.add_log(f'Data saved to {full_path} as JSON.')
		else:
			raise ValueError('For JSON, data should be a dictionary.')
	elif file_type == 'txt':
		if isinstance(data, str):
			with open(full_path, 'w') as file:
				file.write(data)
			report_mod.add_log(f'Data saved to {full_path} as TXT.')
		else:
			raise ValueError('For TXT, data should be a string.')
	else:
		raise ValueError('Unsupported file type. Use "csv", "json", or "txt".')




def load_csv_file(file_path, verbose=False, **kwargs):
	data = pd.read_csv(file_path, **kwargs)

	if verbose:
		report_mod.add_log(f'Data loaded from {file_path}.')
		report_mod.add_log(f'Columns in the data frame: {data.columns}')
	return data


def load_node_dict_from_csv(file_path, verbose=False):
	"""Load a node dictionary from a CSV file.
	
	Expects a CSV with at least two columns:
	- First column: node names
	- Second column: symbols (comma-separated or single value)
	
	Returns:
		dict: Mapping of node_name -> list of symbols
	"""
	if not os.path.exists(file_path):
		raise FileNotFoundError(f"Node dictionary file not found: {file_path}")
	
	df = pd.read_csv(file_path)
	node_dict = {}
	
	for _, row in df.iterrows():
		node_name = str(row.iloc[0])  # First column is node name
		symbols_str = str(row.iloc[1]) if len(row) > 1 else ""  # Second column is symbols
		
		# Parse symbols: split by comma and strip whitespace
		if isinstance(symbols_str, str) and symbols_str and symbols_str != "nan":
			symbols = [s.strip() for s in symbols_str.split(",") if s.strip()]
		else:
			symbols = []
		
		node_dict[node_name] = symbols
	
	if verbose:
		report_mod.add_log(f'Loaded node dictionary from {file_path} with {len(node_dict)} nodes.')
	
	return node_dict


def load_sidm_from_model_csv(cell_line_names, verbose=False):
	"""Load SIDM mapping for user-provided identifiers.
	
	Resolves identifiers through the online resolver (Sanger API primary, 
	Cellosaurus fallback). Function name retained for backwards compatibility.
	
	Args:
		cell_line_names (list): List of identifiers from user's cell_line_file.
		verbose (bool): If True, log detailed matching results.
	
	Returns:
		tuple: (sidm_dict, not_found)
		       sidm_dict: {SIDM -> cell_line_name} for matched entries
		       not_found: list of identifiers that couldn't be matched
	
	Raises:
		ValueError: If no identifiers from the input list can be resolved.
	"""
	sidm_dict, not_found, results = resolve_identifiers_to_sidm(cell_line_names)

	if verbose:
		for result in results:
			if result.status == "resolved":
				report_mod.add_log(f"Matched '{result.input_raw}' -> SIDM {result.sidm} ({result.matched_on})")
			else:
				report_mod.add_log(f"No match found for '{result.input_raw}' via online resolution")
	
	if not sidm_dict:
		raise ValueError(
			f"No cell lines from the provided list could be resolved. "
			f"Not found: {', '.join(not_found)}. "
			f"Please check your cell_line_file names or provide a 'SIDM' column explicitly."
		)
	
	if not_found and verbose:
		report_mod.add_log(
			f"Warning: {len(not_found)} cell line(s) could not be resolved and will be excluded: {', '.join(not_found)}"
		)
	
	return sidm_dict, not_found



def load_sidm_from_modelid(model_ids, model_registry=None, verbose=False):
	"""Load SIDM mapping from ModelID values.

	Maps DepMap ModelID values (ACH-*) to SIDM through the online resolver
	(Sanger API primary, Cellosaurus fallback).

	Args:
		model_ids (list): List of ModelID values (e.g., ['ACH-000001', 'ACH-000002', ...])
		model_registry (str): Kept for backward compatibility; currently unused (online resolution is used).
		verbose (bool): If True, log detailed matching results.

	Returns:
		tuple: (model_to_sidm, not_found)
		       model_to_sidm: {ModelID -> SIDM} for matched entries
		       not_found: unique model_ids that couldn't be matched

	Raises:
		ValueError: If no model IDs from the input list can be resolved.
	"""
	model_to_sidm, not_found = resolve_model_ids_to_sidm(model_ids, model_registry=model_registry)

	if verbose:
		for mid in model_ids:
			normalized = str(mid).strip()
			if normalized in model_to_sidm:
				report_mod.add_log(f"Matched ModelID '{normalized}' -> SIDM {model_to_sidm[normalized]}")
			else:
				report_mod.add_log(f"No match found for ModelID '{normalized}' via online resolution")
	
	if not model_to_sidm:
		raise ValueError(
			f"No ModelID values from the provided list could be resolved. "
			f"Not found: {', '.join(not_found)}. "
			f"Please check your ModelID values."
		)
	
	if not_found and verbose:
		report_mod.add_log(
			f"Warning: {len(not_found)} ModelID(s) could not be resolved and will be excluded: {', '.join(not_found)}"
		)
	
	return model_to_sidm, not_found


def build_alias_map_from_cellfile(cell_line_file, verbose=False):
	"""Resolve a user-provided cell_line_file and build alias -> SIDM mapping.

	Returns:
		tuple: (alias_to_sidm, sidm_dict, not_found, resolution_report)
	"""
	if not cell_line_file:
		raise ValueError("cell_line_file must be provided")

	df = pd.read_csv(cell_line_file)
	sidm_dict, not_found, resolution_report = resolve_sidm_from_dataframe(df)

	alias_map = {}
	if isinstance(resolution_report, dict):
		alias_map = resolution_report.get('alias_to_sidm', {}) or {}

	if verbose:
		report_mod.add_log(f"Built alias map for {len(sidm_dict)} SIDMs with {len(alias_map)} aliases")

	return alias_map, sidm_dict, not_found, resolution_report


def normalize_sample_ids_to_sidm(df, alias_map=None, sidm_list=None, collapse_method='max'):
	"""Normalize DataFrame sample identifiers (columns) to SIDM using provided alias_map.

	Args:
		df: pandas DataFrame with sample identifiers as columns
		alias_map: dict alias -> SIDM
		sidm_list: optional list of SIDMs to restrict output columns
		collapse_method: how to collapse multiple columns mapping to same SIDM ('max'|'mean'|'first')

	Returns:
		(df_normalized, unmapped_list)
	"""
	if df is None or df.empty:
		return df, []

	alias_map = alias_map or {}
	mapped_cols = {}
	unmapped = []

	for col in list(df.columns):
		key = str(col).strip()
		# direct SIDM
		if sidm_list and key in sidm_list:
			mapped_cols.setdefault(key, []).append(col)
			continue
		if key in alias_map:
			mapped_cols.setdefault(alias_map[key], []).append(col)
			continue
		# case-insensitive alias match
		low = key.lower()
		found = None
		for a, s in alias_map.items():
			if a.lower() == low:
				found = s
				break
		if found:
			mapped_cols.setdefault(found, []).append(col)
			continue
		unmapped.append(key)

	if not mapped_cols:
		return df, unmapped

	new_cols = {}
	for sidm, orig_cols in mapped_cols.items():
		if len(orig_cols) == 1:
			new_cols[sidm] = df[orig_cols[0]]
		else:
			if collapse_method == 'max':
				new_cols[sidm] = df[orig_cols].max(axis=1, skipna=True)
			elif collapse_method == 'mean':
				new_cols[sidm] = df[orig_cols].mean(axis=1, skipna=True)
			else:
				new_cols[sidm] = df[orig_cols].iloc[:, 0]

	new_df = pd.DataFrame(new_cols)

	# If sidm_list provided, restrict to those SIDMs (preserve order if possible)
	if sidm_list is not None:
		cols = [s for s in sidm_list if s in new_df.columns]
		new_df = new_df.loc[:, cols]

	return new_df, unmapped
