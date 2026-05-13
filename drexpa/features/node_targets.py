# Drug node targets mapping module

import pandas
import ast
from .helpers import save_file


class NodeTargetMapper:
    """
    Node targets mapping class.
    
    This class handles the mapping of HGNC gene symbols to model node names.
    
    Args:
        chembl_targets_df (pandas.DataFrame): DataFrame containing ChEMBL targets.
        node_dict (dict): Dictionary mapping node names to HGNC symbols.
        node_dict_file (str): Path to file containing node dictionary.
        column_drugname (str): Name of the column containing drug names. Default: 'drug_name'.
        column_target (str): Name of the column containing targets. Default: 'targets'.
        directory_output (str): Directory to save output files.
        verbose (bool): Print verbose output. Default: False.
        save (bool): Save output to file. Default: False.
    
    Methods:
        get_node_targets: Main method to get node targets dataframe.
    """
    
    def __init__(self,
                chembl_targets_df: pandas.DataFrame = None,
                node_dict: dict = None,
                node_dict_file: str = None,
                column_drugname: str = 'drug_name',
                column_target: str = 'targets',
                directory_output: str = None,
                verbose: bool = False,
                save: bool = False):
        
        self.chembl_targets_df = chembl_targets_df
        self.node_dict = node_dict
        self.node_dict_file = node_dict_file
        self.column_drugname = column_drugname
        self.column_target = column_target
        self.directory_output = directory_output
        self.verbose = verbose
        self.save = save
    
    def get_node_targets(self,
                        chembl_targets_df: pandas.DataFrame = None,
                        node_dict: dict = None,
                        column_target: str = None,
                        column_drugname: str = None) -> pandas.DataFrame:
        """
        Get the node targets dataframe by mapping HGNC symbols to node names.
        
        Args:
            chembl_targets_df (pandas.DataFrame): Override chembl_targets_df if provided.
            node_dict (dict): Override node_dict if provided.
            column_target (str): Override column_target if provided.
            column_drugname (str): Override column_drugname if provided.
            
        Returns:
            pandas.DataFrame: DataFrame with node_targets column added.
        """
        if chembl_targets_df is not None:
            self.chembl_targets_df = chembl_targets_df
            
        if node_dict is not None:
            self.node_dict = node_dict
            
        if column_target:
            self.column_target = column_target
        
        if column_drugname:
            self.column_drugname = column_drugname
        
        if self.chembl_targets_df is None:
            raise ValueError('No chembl_targets_df provided. Please provide a ChEMBL targets dataframe.')
        
        if self.verbose:
            print('\nMaking the node_targets dataframe ...')
            print('\nMapping drug HGNC targets to node names...')
            print(f'Using the target column {self.column_target} to map the targets to node names.\n')
        
        if self.node_dict_file and self.node_dict is None:
            # Read the node dictionary file
            self.node_dict = self._read_node_dict_file(self.node_dict_file)
        
        node_targets_df = self._map_target_nodes(
            self.chembl_targets_df,
            self.node_dict
        )
        
        # Drop duplicates to avoid multiple rows for same drug-concentration pair
        # This happens because expanded synergy data may have duplicate drug-concentration combinations
        if self.verbose:
            print(f'Dropping duplicate drug-concentration combinations...')
            print(f'Before deduplication: {len(node_targets_df)} rows')
        
        node_targets_df = node_targets_df.drop_duplicates(subset=[self.column_drugname, 'ChEMBL_conc'], keep='first').reset_index(drop=True)
        
        if self.verbose:
            print(f'After deduplication: {len(node_targets_df)} rows\n')
        
        if self.save and self.directory_output:
            save_file(node_targets_df, self.directory_output, 'drug_node_targets.csv', index=False)
        
        return node_targets_df
    
    def _read_node_dict_file(self, node_dict_file) -> dict:
        """
        Read the node dictionary file.
        
        Args:
            node_dict_file (str): Path to node dictionary file.
            
        Returns:
            dict: Dictionary mapping node names to HGNC symbols.
        """
        if self.verbose:
            print('\nReading the node-gene symbols dictionary from file ...')
            print('Node dictionary file:', node_dict_file)
        node_dict = pandas.read_csv(node_dict_file)
        node_dict = {k: ast.literal_eval(v) for k, v in zip(node_dict.node_name, node_dict.HGNC_symbol)}
        
        if self.verbose:
            print('Node dictionary:', node_dict, '\n')
        
        return node_dict
    
    def _map_target_nodes(self,
                        chembl_targets_df: pandas.DataFrame,
                        node_dict: dict) -> pandas.DataFrame:
        """
        Map each target to its corresponding node name.
        
        Args:
            chembl_targets_df (pandas.DataFrame): DataFrame with ChEMBL targets.
            node_dict (dict): Dictionary mapping node names to HGNC symbols.
            
        Returns:
            pandas.DataFrame: DataFrame with node_targets column added.
        """
        # Clean the ChEMBL targets dataframe
        chembl_targets_df = self._clean_ChEMBLtargets(chembl_targets_df)
        
        if self.verbose:
            print('Mapping HGNC targets to node names ...\n')
        
        def map_targets(column_target):
            matching_nodenames = []
            for node_key, hgnc_values in node_dict.items():
                for value in hgnc_values:
                    if value in column_target:
                        matching_nodenames.append(node_key)
            node_targets = list(set(matching_nodenames))  # Get unique node names
            return node_targets if node_targets else 'No matching targets'
        
        chembl_targets_df['node_targets'] = chembl_targets_df[self.column_target].apply(map_targets)
        
        if self.verbose:
            print('HGNC targets mapped to node names successfully.\n')
            print('Node targets dataframe:', chembl_targets_df.head(5))
            print('Shape of node_targets dataframe:', chembl_targets_df.shape)
        
        return chembl_targets_df
    
    def _clean_ChEMBLtargets(self,
                            chembl_targets_df: pandas.DataFrame,
                            column_drugname: str = None,
                            column_target: str = None) -> pandas.DataFrame:
        """
        Clean the ChEMBL targets dataframe.
        
        Args:
            chembl_targets_df (pandas.DataFrame): DataFrame with ChEMBL targets.
            column_drugname (str): Override column_drugname if provided.
            column_target (str): Override column_target if provided.
            
        Returns:
            pandas.DataFrame: Cleaned DataFrame.
        """
        if self.verbose:
            print('\nCleaning ChEMBLtargets dataframe ...\n')
            print('Keeping only the columns:', self.column_drugname, self.column_target)
        
        try:
            chembl_targets_df = chembl_targets_df[[self.column_drugname, self.column_target, 'ChEMBL_conc']]
        except KeyError:
            chembl_targets_df = chembl_targets_df[[self.column_drugname, self.column_target]]
        
        if chembl_targets_df[self.column_target].isnull().sum() > 0:
            chembl_targets_df.loc[:, self.column_target] = chembl_targets_df[self.column_target].replace('', None).fillna('not found')
        
        if self.verbose:
            # If there are 'not found' values in the targets column
            if chembl_targets_df[self.column_target].str.contains('not found').any():
                print('Drugs with not found values in targets column.')
                print('Number of drugs with not found targets:', 
                    chembl_targets_df[self.column_target].str.contains('not found').sum())
            else:
                print('No drugs with missing values in the targets column.\n')
        
        chembl_targets_df = chembl_targets_df[chembl_targets_df[self.column_target] != 'not found']
        
        if self.verbose:
            print('Shape of ChEMBLtargets dataframe:', chembl_targets_df.shape)
        
        return chembl_targets_df
