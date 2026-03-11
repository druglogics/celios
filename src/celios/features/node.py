"""Node processing module (migrated into features).

This file contains the Node class implementation relocated into the
`celios.features` package. It mirrors the previous top-level `celios.node`
functionality but keeps imports local to the features package.
"""

import pandas
from pathlib import Path
import re
from ..utils.io import save_file
from ..utils import report as report_mod

class Node:
    """
    Node class.
    
        This class contains all the methods to process the file of node names and genes (or entities).
        Starting from a file with node names and genes, the class creates a dictionary with the node names as keys and the genes as values.
        Genes can be mapped to HGNC symbols using the HGNC symbols file, if specified. This will add more values to the dictionary.
        The node dictionary can be saved as a JSON or CSV file.

    Args:
    nodesymbols_file (str): file containing node names and genes in HGNC symbols
    hgnc_symbols_file (str): file containing HGNC symbols. Default: None. *
    directory_output (str): directory to save the node dictionary file.
    output_file (str): name of the output node dictionary file. Default: 'node_dict.csv'. **
    verbose (bool): print messages. Default: False.
        
        * If the HGNC symbols file is not specified, the node dictionary will be made only with the node names and symbols in the original file.
        ** The name of the node dictionary file can be changed. It is recommended to use the default value.
            The default name is 'node_dict.csv', when the hgnc_symbols_file is not specified.
            The default name is 'node_HGNC_dict.csv', when the hgnc_symbols_file is specified.

    Methods:
    """

    def __init__(self,
                node_input: object = None,
                hgnc_symbols_file: str = None,
                directory_output: str = None,
                verbose: bool = False,
                include_alias_prev: bool = True,
                ) -> None:
        """
        New constructor accepts a single `node_input` which may be:
        - a path (string) to a file (CSV/TSV or TXT)
        - a list of node names
        - a pandas.DataFrame with node names and symbols

        A separate helper `_detect_and_load_node_input` will inspect
        `node_input` and populate internal attributes used elsewhere in
        the class (`node_list`, `nodesymbols_df`, `node_list_file`,
        `nodesymbols_file`).
        """
        # Single unified input (path, list or DataFrame)
        self.node_input = node_input

        # Keep older-style attributes but they are now populated by the
        # detection loader if a single `node_input` was provided.
        self.node_list_file = None
        self.nodesymbols_file = None
        self.nodesymbols_df = None
        self.node_list = None

        # Use default headers internally (no longer accepted as ctor args)
        self.node_name_header = 'node_name'
        self.symbol_header = 'HGNC_symbol'

        # Other options
        self.hgnc_symbols_file = hgnc_symbols_file
        self.directory_output = directory_output
        self.verbose = verbose
        # Whether to include alias_symbol and prev_symbol in mapped results
        self.include_alias_prev = include_alias_prev
        self.node_dict = None
        # Fixed default filename; user no longer provides this argument
        self.node_dict_file = 'node_dict.csv'
        self.hgnc_node_dict = None
        # Cache for HGNC set to avoid reloading the file repeatedly
        self._HGNC_set_cache = None
        # Flag to avoid re-processing the same input repeatedly
        self._node_input_processed = False

        '''
        node_name_header (str): header of the node names column. Default: 'node_name' *
        symbol_header (str): header of the HGNC symbols column. Default: 'HGNC_symbol' *

        * The node and symbol headers will be used to rename the headers in the node names file.
        It is recommended to use the default values. This ensures that other modules can read and use the file correctly.
        '''

    #/////////////////////////////////////////////////////////////////////////////////////////
    # BRANCH: node_list_maker
    #/////////////////////////////////////////////////////////////////////////////////////////

    def _detect_and_load_node_input(self) -> None:
        """
        Inspect `self.node_input` (if provided) and populate the internal
        attributes the rest of the class expects: `node_list`,
        `nodesymbols_df`, `node_list_file`, or `nodesymbols_file`.

        This method is idempotent and sets `_node_input_processed` to
        avoid re-processing.
        """
        if getattr(self, '_node_input_processed', False):
            return

        inp = getattr(self, 'node_input', None)

        # Nothing to detect; leave attributes as-is
        if inp is None:
            self._node_input_processed = True
            return

        # If a list was provided, treat it as a node list
        if isinstance(inp, list):
            self.node_list = list(inp)
            self._node_input_processed = True
            return

        # If a DataFrame provided, use it directly
        if isinstance(inp, pandas.DataFrame):
            df_full = inp.copy()
            # If single-column DataFrame -> treat as list
            if df_full.shape[1] == 1:
                self.node_list = df_full.iloc[:, 0].astype(str).tolist()
                self._node_input_processed = True
                return

            # Use first two columns as node name and symbol column regardless of header names
            df2 = df_full.iloc[:, :2].copy()

            # Inspect sample values in second column to detect UniProt IDs
            sample = df2.iloc[:min(10, df2.shape[0]), 1].dropna().astype(str).tolist()
            uniprot_pattern = re.compile(r'^[A-Za-z][A-Za-z0-9]{5}$')
            if sample:
                matches = sum(1 for v in sample if uniprot_pattern.match(v))
                # If majority of sample rows look like UniProt IDs, mark as uniprot
                if matches >= max(1, len(sample) // 2 + (len(sample) % 2)):
                    self.symbol_header = 'uniprot'
                else:
                    self.symbol_header = 'HGNC_symbol'
            else:
                # No sample values; default to HGNC
                self.symbol_header = 'HGNC_symbol'

            df2.columns = [self.node_name_header, self.symbol_header]
            self.nodesymbols_df = df2
            self._node_input_processed = True
            return

        # If a path string was provided, inspect the file
        if isinstance(inp, str):
            path = inp
            # Try reading a small sample with pandas first (handles csv/tsv)
            try:
                sample_df = pandas.read_csv(path, header=0, nrows=10)

                # If only one column, treat as node list
                if sample_df.shape[1] == 1:
                    # Read full file to get all rows
                    df_full = pandas.read_csv(path, header=0)
                    self.node_list = df_full.iloc[:, 0].astype(str).tolist()
                    self.node_list_file = path
                    self._node_input_processed = True
                    return

                # If two or more columns: use first two as name/symbol regardless of header names
                if sample_df.shape[1] >= 2:
                    # Inspect sample values in second column to detect UniProt IDs
                    sample_vals = sample_df.iloc[:min(10, sample_df.shape[0]), 1].dropna().astype(str).tolist()
                    uniprot_pattern = re.compile(r'^[A-Za-z][A-Za-z0-9]{5}$')
                    if sample_vals:
                        matches = sum(1 for v in sample_vals if uniprot_pattern.match(v))
                        if matches >= max(1, len(sample_vals) // 2 + (len(sample_vals) % 2)):
                            symbol_header_guess = 'uniprot'
                        else:
                            symbol_header_guess = 'HGNC_symbol'
                    else:
                        symbol_header_guess = 'HGNC_symbol'

                    # Read full dataframe and pick first two columns
                    df_full = pandas.read_csv(path, header=0)
                    df2 = df_full.iloc[:, :2].copy()
                    df2.columns = [self.node_name_header, symbol_header_guess]
                    # If guessed uniprot, set symbol_header accordingly
                    self.symbol_header = symbol_header_guess
                    self.nodesymbols_df = df2
                    self.nodesymbols_file = path
                    self._node_input_processed = True
                    return

            except Exception:
                # Fall back to treating it as a plain text list file
                try:
                    with open(path, 'r') as f:
                        lines = [line.strip() for line in f.readlines() if line.strip()]
                    self.node_list = lines
                    self.node_list_file = path
                    self._node_input_processed = True
                    return
                except Exception:
                    # Give up and let other code raise more specific errors later
                    self._node_input_processed = True
                    return

        # Mark processed even if nothing matched to avoid loops
        self._node_input_processed = True


    #//////////////////////////////////////////////////////////
    def _map_node_list_to_HGNC(self) -> list:
        '''
        Function: map the node list to HGNC symbols.

        Returns:
            - hgnc_node_list: list containing HGNC symbols.
        '''
        # Load raw node list
        node_list = self._load_node_list()

        # Build a node_dict like the one used elsewhere (values as lists)
        node_dict = {n: [n] for n in node_list}

        # If an HGNC set is provided, map using existing helpers
        if self.hgnc_symbols_file:
            if self.symbol_header == 'HGNC_symbol':
                return self._map_nodes_to_HGNC(node_dict)
            elif self.symbol_header == 'uniprot':
                return self._map_uniprot_to_HGNC(node_dict)

        # If no HGNC file provided, return the simple identity mapping
        return node_dict


    #//////////////////////////////////////////////////////////
    def _load_node_list(self) -> list:
        '''
        Function: get the node list from the node names file.

        Returns:
            - node_list: list containing node names.
        '''
        # Ensure single `node_input` (if provided) is detected and loaded
        self._detect_and_load_node_input()

        if self.verbose:
            if self.node_list_file:
                print(f'Processing node list file: {self.node_list_file}\n')
            else:
                count = len(self.node_list) if self.node_list is not None else 0
                print(f'Processing node list (in-memory) with {count} entries\n')

        # If a direct node_list was provided to the constructor, use it
        if self.node_list is not None:
            return list(self.node_list)

        # Check if node_list_file is provided
        if self.node_list_file is None:
            raise ValueError("node_list_file is not provided.")
        
        # Inspect the node list file. It can be a csv or a list
        if self.node_list_file.endswith('.csv'):
            # Load the node list from a csv file
            node_list = self._load_node_names_list()
        else:
            # Load the node list from a txt file
            with open(self.node_list_file, 'r') as f:
                node_list = [line.strip() for line in f.readlines()]
            if self.verbose:
                report_mod.add_log('Node list created successfully')
                report_mod.add_log(f'Length of node list: {len(node_list)}')
        return node_list

     #//////////////////////////////////////////////////////////
    def _load_node_names_list(self) -> list:
        '''
        Function: load the node names file into a list.
        '''
        # Load the node names file
        # When used to load a plain node list (CSV), read from `node_list_file`.
        # Accept either a CSV with a named header column (default 'node_name')
        # or a single-column CSV without header.
        try:
            df = pandas.read_csv(self.node_list_file, header=0)
            # If the expected header exists, use it
            if self.node_name_header in df.columns:
                return df[self.node_name_header].astype(str).tolist()
            # If only one column present and there are data rows, use it
            if df.shape[1] == 1 and df.shape[0] > 0:
                return df.iloc[:, 0].astype(str).tolist()
            # Otherwise fall through to trying without header
        except Exception:
            pass

        # Try reading CSV without header and use first column
        df = pandas.read_csv(self.node_list_file, header=None)
        return df.iloc[:, 0].astype(str).tolist()


    #/////////////////////////////////////////////////////////////////////////////////////////
    # BRANCH: node_dict_maker
    #/////////////////////////////////////////////////////////////////////////////////////////

    def get_node_dict(self, manual_symbols_file: object = None) -> dict:
        '''
        Function: get the node dictionary from the node names file.

        Returns:
            - node_dict: dictionary containing node names and HGNC symbols.
        '''
        # Ensure `node_input` has been detected and loaded into internal attributes
        self._detect_and_load_node_input()

        # If a plain node list is provided (directly or as a file), build the node dict from it
        if (getattr(self, 'node_list', None) is not None) or self.node_list_file:

            # Load the list and construct a node_dict where each value is a list
            node_list = self._load_node_list()
            self.node_dict = {n: [n] for n in node_list}
            # Map to HGNC if an HGNC symbols file is provided
            hgnc_node_dict = None
            if self.hgnc_symbols_file:
                hgnc_node_dict = self._map_node_dict(self.node_dict)
        else:
            if self.verbose:
                report_mod.add_log(f'Processing node names file: {self.nodesymbols_file}')

            # Make the node dictionary from the nodenames file
            self.node_dict = self._make_node_dict()

            # Map the nodes to HGNC symbols
            hgnc_node_dict = None
            if self.hgnc_symbols_file:
                hgnc_node_dict = self._map_node_dict(self.node_dict)

        # If a manual symbols file was provided, apply overrides to the
        # appropriate dictionary (mapped HGNC dict if present, otherwise
        # the original node_dict).
        if manual_symbols_file:
            if hgnc_node_dict is not None:
                hgnc_node_dict = self._manual_symbols(hgnc_node_dict, manual_symbols_file)
            else:
                self.node_dict = self._manual_symbols(self.node_dict, manual_symbols_file)

        # Save the node dictionary (after manual overrides if any)
        if self.directory_output:
            # Make the node dictionary into a dataframe
            if self.hgnc_symbols_file and hgnc_node_dict is not None:
                self.node_dict_file = 'node_HGNC_dict.csv'
                node_dict_df = pandas.DataFrame(list(hgnc_node_dict.items()), columns=[self.node_name_header, self.symbol_header])
            else:
                node_dict_df = pandas.DataFrame(list(self.node_dict.items()), columns=[self.node_name_header, self.symbol_header])
            save_file(node_dict_df, self.directory_output, self.node_dict_file, file_type='csv', index=False)

        if self.hgnc_symbols_file:
            return hgnc_node_dict
        else:
            return self.node_dict

    #//////////////////////////////////////////////////////////
    def _make_node_dict(self) -> dict:
        '''
        Function: make a dictionary from the node names dataframe.

        Returns:
            - node_dict: dictionary containing node names and HGNC symbols.
        '''
        # Load the node names file
        nodensymbols_df = self._load_node_symbols_file()

        # Make the node dictionary
        node_dict = dict(zip(nodensymbols_df[self.node_name_header], nodensymbols_df[self.symbol_header]))

        # Process the node dictionary
        node_dict = self._process_node_dict(node_dict)

        if self.verbose:
            report_mod.add_log('\nNode dictionary created successfully')
            report_mod.add_log(f'Length of node dictionary: {len(node_dict)}')

        return node_dict
    
    #//////////////////////////////////////////////////////////
    def _process_node_dict(self, node_dict: dict) -> dict:
        '''
        Function: process the node dictionary. Split value strings with ", ".

        Returns:
            - node_dict: dictionary containing node names and HGNC symbols.
        '''
        # Split value strings with ", "
        for key, value in node_dict.items():
            if isinstance(value, str):
                val = value.split(', ')
                node_dict.update({key: val})
        
        # Check for missing values
        self._missing_values(node_dict)
        
        return node_dict
    
    #//////////////////////////////////////////////////////////
    def _missing_values(self, node_dict: str) -> None:
        '''
        Function: check for missing values in the node dictionary.

        Parameters:
            - node_dict: dictionary containing node names and HGNC symbols.
        '''
        missing = [key for key, value in node_dict.items() if value == '']
        
        if self.verbose:
            if missing:
                report_mod.add_log(f'Missing values in the node dictionary: {missing}')
                report_mod.add_log(f'Number of missing values: {len(missing)}')
            else:
                report_mod.add_log('No missing values in the node dictionary.')

    #//////////////////////////////////////////////////////////
    def _load_node_symbols_file(self) -> pandas.DataFrame:
        '''
        Function: load the node names csv file into a dataframe and rename the headers.
        '''
        # If a DataFrame was provided directly, use it
        if getattr(self, 'nodesymbols_df', None) is not None:
            nodesymbols_df = self.nodesymbols_df.copy()

            # If it already contains expected columns, select them
            if self.node_name_header in nodesymbols_df.columns and self.symbol_header in nodesymbols_df.columns:
                nodesymbols_df = nodesymbols_df[[self.node_name_header, self.symbol_header]]
            # If it's two columns, rename to expected headers
            elif nodesymbols_df.shape[1] == 2:
                nodesymbols_df.columns = [self.node_name_header, self.symbol_header]
            else:
                raise ValueError('Provided nodenames_df must contain two columns or columns matching node_name_header and symbol_header')
        else:
            # Load the node names file
            nodesymbols_df = pandas.read_csv(self.nodesymbols_file, header=0)

            # Rename the headers using the values of the node_name_header and HGNC_symbol_header
            nodesymbols_df.columns = [self.node_name_header, self.symbol_header]

        if self.verbose:
            report_mod.add_log('Node names file loaded successfully')
            report_mod.add_log(f'Shape of node names dataframe: {nodesymbols_df.shape}')
            report_mod.add_log(f'Node names dataframe columns: {nodesymbols_df.columns}')

        return nodesymbols_df
    

    #/////////////////////////////////////////////////////////////////////////////////////////
    # BRANCH: node_HGNC_dict_maker
    #/////////////////////////////////////////////////////////////////////////////////////////

    def _map_nodes_to_HGNC(self, node_dict: dict) -> dict:
        '''
        Function: map the nodes in the node dictionary to the HGNC symbols in the HGNC set.

        Returns:
            - hgnc_node_dict: dictionary containing node names and HGNC symbols.        
        '''
        # Get the HGNC set
        HGNC_set = self._get_HGNCset()

        hgnc_node_dict = {} 
        if self.verbose:
            report_mod.add_log('Mapping nodes to HGNC symbols ...')
        for key, values_list in node_dict.items():
            mapped_values = []

            for value in values_list:
                matching_rows = HGNC_set[HGNC_set.apply(lambda x: value in x.values, axis=1)]

                if not matching_rows.empty:
                    for _, row in matching_rows.iterrows():
                        mapped_symbols = [row['symbol']]
                        if getattr(self, 'include_alias_prev', True):
                            if pandas.notna(row['alias_symbol']):
                                mapped_symbols.extend(row['alias_symbol'].split(', '))
                            if pandas.notna(row['prev_symbol']):
                                mapped_symbols.extend(row['prev_symbol'].split(', '))
                        mapped_values.extend(mapped_symbols)

            hgnc_node_dict[key] = mapped_values
        
        if self.verbose:
            report_mod.add_log('Nodes mapped to HGNC symbols successfully')
            report_mod.add_log('Getting missing nodes ...')
        # Complete values for of missing nodes
        hgnc_node_dict = self._get_missing_nodes(hgnc_node_dict, node_dict)

        return hgnc_node_dict
    
    #//////////////////////////////////////////////////////////
    def _map_uniprot_to_HGNC(self, node_dict: dict) -> dict:
        '''
        Function: map the uniprot IDs in the node dictionary to the HGNC symbols in the HGNC set.
        Returns:
            - hgnc_node_dict: dictionary containing node names and HGNC symbols.
        '''
        # Get the HGNC set
        HGNC_set = self._get_HGNCset()
        hgnc_node_dict = {}
        if self.verbose:
            report_mod.add_log('Mapping uniprot IDs to HGNC symbols ...')
        for key, values_list in node_dict.items():
            mapped_values = []

            for value in values_list:
                # Look for the UniProt ID in the uniprot_ids column
                matching_rows = HGNC_set[HGNC_set['uniprot_ids'].notna() & 
                                        HGNC_set['uniprot_ids'].str.contains(value, na=False)]

                if not matching_rows.empty:
                    for _, row in matching_rows.iterrows():
                        mapped_symbols = [row['symbol']]
                        mapped_values.extend(mapped_symbols)

            hgnc_node_dict[key] = mapped_values

        if self.verbose:
            report_mod.add_log('Uniprot IDs mapped to HGNC symbols successfully')
            report_mod.add_log('Getting missing nodes ...')
        # Complete values for missing nodes
        hgnc_node_dict = self._get_missing_nodes(hgnc_node_dict, node_dict)

        return hgnc_node_dict

    #//////////////////////////////////////////////////////////
    def _get_HGNCset(self) -> pandas.DataFrame:
        '''
        Function: get the HGNC set from the HGNC symbols file.
        '''
        # Cache the processed HGNC set so repeated calls are cheap
        if getattr(self, '_HGNC_set_cache', None) is not None:
            return self._HGNC_set_cache

        # Always load the HGNC symbols file (the HGNC file is required for mapping)
        if not self.hgnc_symbols_file:
            raise ValueError('hgnc_symbols_file is not provided but is required to build HGNC set')

        # Load and process from file, then cache
        self._HGNC_set_cache = self._process_HGNC_set()
        return self._HGNC_set_cache

    #//////////////////////////////////////////////////////////
    def _map_node_dict(self, node_dict: dict) -> dict:
        '''
        Centralized mapping selector: maps a node_dict according to `self.symbol_header`.
        '''
        if not self.hgnc_symbols_file:
            return node_dict

        if self.symbol_header == 'uniprot':
            return self._map_uniprot_to_HGNC(node_dict)
        else:
            return self._map_nodes_to_HGNC(node_dict)
    
    #//////////////////////////////////////////////////////////
    def _keep_original_nodenames(self, hgnc_node_dict, node_dict) -> None:
        '''
        Function: keep the original node symbol (value) on node_dict if no HGNC symbol is found.
        '''
        for key, values in hgnc_node_dict.items():
            if not values:
                hgnc_node_dict[key] = node_dict[key]

        if self.verbose:
            report_mod.add_log('Original node symbol kept for missing HGNC symbols')

    #//////////////////////////////////////////////////////////
    def _get_missing_nodes(self, hgnc_node_dict, node_dict) -> dict:
        '''
        Function: get the list of missing nodes.
        '''
        missing_values = [key for key, values in hgnc_node_dict.items() if not values]
        if self.verbose:
            report_mod.add_log(f'Number of missing nodes: {len(missing_values)}')
            report_mod.add_log(f'Missing nodes: {missing_values}')

        # Keep original node symbol if no HGNC symbol is found
        if missing_values is not None:
            if self.verbose:
                report_mod.add_log('Keeping original node symbol if no HGNC symbol is found ...')
            self._keep_original_nodenames(hgnc_node_dict, node_dict)

        return hgnc_node_dict

    #//////////////////////////////////////////////////////////
    def _manual_symbols(self, node_dict: dict, manual_symbols: object) -> dict:
        '''
        Function: apply manual symbol overrides from a CSV file or DataFrame.

        Parameters:
            - node_dict: dictionary mapping node_name -> list of symbols
            - manual_symbols: path to a CSV file or a pandas.DataFrame containing
                              columns matching `self.node_name_header` and
                              `self.symbol_header` (e.g. 'node_name' and
                              'HGNC_symbol'). Values in the symbol column
                              may contain multiple symbols separated by commas.

        Returns:
            - node_dict: dictionary with values overridden by manual mappings
        '''
        if manual_symbols is None:
            return node_dict

        # Load manual mappings from provided input (path or DataFrame)
        if isinstance(manual_symbols, pandas.DataFrame):
            df = manual_symbols.copy()
        else:
            # Expect a file path
            try:
                df = pandas.read_csv(manual_symbols, header=0)
            except Exception as e:
                raise FileNotFoundError(f'Could not read manual symbols file: {manual_symbols}') from e

        # Expect columns with the node name header and symbol header
        expected_name_col = self.node_name_header
        expected_sym_col = self.symbol_header
        if expected_name_col not in df.columns or expected_sym_col not in df.columns:
            raise ValueError(f'Manual symbols file must contain columns: {expected_name_col}, {expected_sym_col}')

        # Iterate rows and override/add entries
        for _, row in df.iterrows():
            name = row[expected_name_col]
            val = row[expected_sym_col]
            if pandas.isna(name) or pandas.isna(val):
                continue

            # Split on comma (and optional whitespace) into a list
            if isinstance(val, str):
                symbols = [s.strip() for s in re.split(r',\s*', val) if s.strip()]
            else:
                symbols = [str(val)]

            if name in node_dict:
                node_dict[name] = symbols
                if self.verbose:
                    report_mod.add_log(f'Manual override applied for node: {name} -> {symbols}')
            else:
                # Skip manual overrides for names not present in the generated node_dict
                if self.verbose:
                    report_mod.add_log(f'Manual override skipped for unknown node: {name}')

        return node_dict

    

    #/////////////////////////////////////////////////////////////////////////////////////////
    # BRANCH: HGNCset_maker
    #/////////////////////////////////////////////////////////////////////////////////////////

    def _process_HGNC_set(self) -> pandas.DataFrame:
        '''
        Function: process the HGNC symbols file.

        Returns:
            - HGNC_set: dataframe containing 'symbol,' 'alias_symbol,', 'prev_symbol', and 'uniprot_ids' columns.
        '''
        # Load the HGNC symbols file
        HGNC_set = self._load_HGNC_set()

        # Replace values in columns from '|' separation to commas (',') separation
        # print('Replacing values in columns from "|" separation to commas (",") separation ...')
        HGNC_set['alias_symbol'] = HGNC_set['alias_symbol'].apply(lambda x: ', '.join(x.split('|')) if pandas.notna(x) and '|' in x else x)
        HGNC_set['prev_symbol'] = HGNC_set['prev_symbol'].apply(lambda x: ', '.join(x.split('|')) if pandas.notna(x) and '|' in x else x)
        HGNC_set['uniprot_ids'] = HGNC_set['uniprot_ids'].apply(lambda x: ', '.join(x.split('|')) if pandas.notna(x) and '|' in x else x)
        # print('Values replaced successfully\n')
        # print('\nHGNC symbols file processed successfully\n')

        return HGNC_set
    
    #//////////////////////////////////////////////////////////
    def _load_HGNC_set(self) -> pandas.DataFrame:
        '''
        Function: load HGNC symbols file into a dataframe.

        Returns:
            - HGNC_set: dataframe containing 'symbol,' 'alias_symbol,' 'prev_symbol', and 'uniprot_ids' columns.
        '''
        # Read the HGNC symbols csv file into a dataframe
        HGNC_set = pandas.read_csv(self.hgnc_symbols_file, sep='\t', low_memory=False)
        # print('HGNC symbols file loaded successfully\n')

        # Keep the columns of interest
        columns_to_keep = ['symbol', 'alias_symbol', 'prev_symbol', 'uniprot_ids']
        # print('Keeping columns:', columns_to_keep)
        HGNC_set = HGNC_set[columns_to_keep]
        # print('Shape of HGNC dataframe:', HGNC_set.shape ,'\n')

        return HGNC_set

    @classmethod
    def from_sif(cls,
                 sif_path: object = None,
                 hgnc_symbols_file: object = None,
                 directory_output: str = None,
                 include_alias_prev: bool = False,
                 verbose: bool = False,
                 manual_symbols_file: object = None):
        """
        Convenience constructor and runner: locate SIF and HGNC files (if not
        provided), extract unique node names from the SIF, run mapping, and
        optionally save results to `directory_output`.

        Returns a tuple `(node_instance, mapped_dict)`.
        """
        # Require explicit sif and hgnc paths from callers. These methods are
        # intended for pipeline use where the user provides the files.
        if not sif_path:
            raise ValueError('`sif_path` is required for from_sif()')
        if not hgnc_symbols_file:
            raise ValueError('`hgnc_symbols_file` is required for from_sif()')

        sif_p = Path(sif_path)
        if not sif_p.exists():
            raise FileNotFoundError(f'SIF file not found: {sif_p}')

        hgnc_p = Path(hgnc_symbols_file)
        if not hgnc_p.exists():
            raise FileNotFoundError(f'HGNC symbols file not found: {hgnc_p}')

        # Extract nodes from SIF using canonical parser; import lazy way
        try:
            from .sifbase import extract_nodes as sif_extract_nodes
            node_list = sorted(sif_extract_nodes(str(sif_p)))
        except Exception:
            # Fallback basic extraction
            nodes = set()
            with open(sif_p, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split()
                    if len(parts) >= 3:
                        nodes.add(parts[0])
                        nodes.add(parts[2])
            node_list = sorted(nodes)

        # Delegate to from_object to perform mapping and saving
        mapped = cls.from_object(node_input=node_list,
                 hgnc_symbols_file=str(hgnc_p),
                 directory_output=directory_output,
                 include_alias_prev=include_alias_prev,
                 verbose=verbose,
                 manual_symbols_file=manual_symbols_file)

        # Return Node instance and mapped dict for convenience
        node = cls(node_input=node_list,
                   hgnc_symbols_file=str(hgnc_p),
                   directory_output=directory_output,
                   verbose=verbose,
                   include_alias_prev=include_alias_prev)
        return node, mapped

    @classmethod
    def from_object(cls,
                    node_input: object = None,
                    hgnc_symbols_file: str = None,
                    directory_output: str = None,
                    include_alias_prev: bool = False,
                    verbose: bool = False,
                    manual_symbols_file: object = None):
        """
        Construct a Node from a single `node_input` and return the mapped
        node dictionary. `node_input` may be a list, a pandas.DataFrame, or
        a path to a file; `_detect_and_load_node_input()` will handle it.

        Additional mapping options (`hgnc_symbols_file`, `directory_output`,
        `include_alias_prev`, `verbose`) are forwarded to the Node.
        """
        if node_input is None:
            raise ValueError('No input provided. Pass the input via `node_input`.')

        node = cls(node_input=node_input,
                   hgnc_symbols_file=hgnc_symbols_file,
                   directory_output=directory_output,
                   verbose=verbose,
                   include_alias_prev=include_alias_prev)

        mapped = node.get_node_dict(manual_symbols_file=manual_symbols_file)
        return mapped


#/////////////////////////////////////////////////////////////////////////////////////////
