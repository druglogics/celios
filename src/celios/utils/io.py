import pandas as pd
import json
import os
from . import report as report_mod


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
	"""Load SIDM (SangerModelID) mapping from the bundled Model.csv file.
	
	This function provides an authoritative, stable mapping from cell-line names to
	SangerModelID values. It does not depend on activity file format, making it
	immune to changes in upstream data sources.
	
	Args:
		cell_line_names (list): List of cell-line names from user's cell_line_file.
		verbose (bool): If True, log detailed matching results.
	
	Returns:
		dict: Mapping of SIDM -> cell_line_name for matched cell lines.
		      Also returns dict with 'matched' and 'not_found' keys if verbose.
	
	Raises:
		FileNotFoundError: If Model.csv is not found in the package.
		ValueError: If no cell lines from the input list match any in Model.csv.
	"""
	import re
	from pathlib import Path
	
	# Locate Model.csv in the package
	package_dir = Path(__file__).parent.parent / 'features'
	model_csv_path = package_dir / 'Model.csv'
	
	if not model_csv_path.exists():
		raise FileNotFoundError(
			f"Model.csv not found at {model_csv_path}. "
			"Please reinstall celios to ensure the package includes the Model.csv registry."
		)
	
	# Load Model.csv
	model_df = pd.read_csv(model_csv_path)
	
	# Normalization function: uppercase and remove non-alphanumeric
	def _normalize_name(name):
		return re.sub(r'[^A-Za-z0-9]', '', str(name).upper()).strip()
	
	# Build a lookup map: normalized_name -> (original_name_from_model, sidm)
	model_lookup = {}
	for _, row in model_df.iterrows():
		cell_line = str(row['CellLineName']) if pd.notna(row['CellLineName']) else None
		sidm = str(row['SangerModelID']) if pd.notna(row['SangerModelID']) else None
		
		if cell_line and sidm:
			normalized = _normalize_name(cell_line)
			# Store both the original name and SIDM for reference
			model_lookup[normalized] = {'original_name': cell_line, 'sidm': sidm}
	
	# Match input cell-line names against the model lookup
	sidm_dict = {}  # Maps SIDM -> cell_line_name
	not_found = []
	
	for cname in cell_line_names:
		normalized = _normalize_name(cname)
		if normalized in model_lookup:
			entry = model_lookup[normalized]
			sidm = entry['sidm']
			sidm_dict[sidm] = cname  # Use the original input name
			if verbose:
				report_mod.add_log(f"Matched '{cname}' -> SIDM {sidm}")
		else:
			not_found.append(cname)
			if verbose:
				report_mod.add_log(f"No match found for '{cname}' in Model.csv")
	
	if not sidm_dict:
		raise ValueError(
			f"No cell lines from the provided list could be found in Model.csv. "
			f"Not found: {', '.join(not_found)}. "
			f"Please check your cell_line_file names or provide a 'SIDM' column explicitly."
		)
	
	if not_found and verbose:
		report_mod.add_log(
			f"Warning: {len(not_found)} cell line(s) not found in Model.csv and will be excluded: {', '.join(not_found)}"
		)
	
	return sidm_dict, not_found
