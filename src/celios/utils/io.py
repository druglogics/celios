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
