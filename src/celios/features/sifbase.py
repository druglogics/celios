"""SIF utilities (canonical implementation under `celios.features`).

Provides small helpers to parse SIF files and extract node lists.
"""

from typing import Set, Tuple
import pandas as pd
from ..utils import report as report_mod

#///////////////////////////////////////////////////////////////////////////////////////////////////////
# MAIN FEATURE: extract_nodes
#///////////////////////////////////////////////////////////////////////////////////////////////////////

#////////////////////////////////////////////////////////////////////////////
def _read_sif_to_dataframe(sif_path: str) -> pd.DataFrame:
	"""Read a SIF file into a DataFrame with columns ['Source','Interaction','Target'].

	The function skips blank lines and comments (lines starting with '#').
	Multiple targets on a single line are expanded into separate rows.
	"""
	data = []
	with open(sif_path, 'r', encoding='utf-8') as fh:
		for line in fh:
			line = line.strip()
			if not line or line.startswith('#'):
				continue
			parts = line.split()
			if len(parts) < 3:
				continue
			source = parts[0]
			interaction = parts[1]
			targets = parts[2:]
			for target in targets:
				data.append([source, interaction, target])

	df = pd.DataFrame(data, columns=['Source', 'Interaction', 'Target'])
	return df

#////////////////////////////////////////////////////////////////////////////
def _detect_input_output_nodes(df: pd.DataFrame) -> Tuple[Set[str], Set[str]]:
	"""Return (input_nodes, output_nodes) computed from a SIF DataFrame.

	Input nodes are sources that never appear as targets. Output nodes are
	targets that never appear as sources.
	"""
	sources = set(df['Source'].dropna().astype(str))
	targets = set(df['Target'].dropna().astype(str))
	input_nodes = sources - targets
	output_nodes = targets - sources
	return input_nodes, output_nodes

#////////////////////////////////////////////////////////////////////////////
    # MAIN FUNCTION
#///////////////////////////////////////////////////////////////////////////
def extract_nodes(sif_path: str) -> Set[str]:
	"""Extract all nodes (sources U targets) from a SIF file and report counts.

	The internal helpers are private (prefixed with '_'). Only
	``extract_nodes`` is part of the public pipeline API.
	"""
	df = _read_sif_to_dataframe(sif_path)
	input_nodes, output_nodes = _detect_input_output_nodes(df)
	report_mod.add_log(f'Identified {len(input_nodes)} input nodes and {len(output_nodes)} output nodes.')
	sources = set(df['Source'].dropna().astype(str))
	targets = set(df['Target'].dropna().astype(str))
	return sources.union(targets)
