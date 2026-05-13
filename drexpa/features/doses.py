# Drug doses processing module

import pandas
from .helpers import save_file
from .chembl_ids import ChEMBLIDResolver


class DoseProcessor:
    """
    Drug dose processing class.
    
    This class handles the processing of drug doses and their integration with ChEMBL IDs.
    Supports both single-drug and dual-drug (synergy) data formats.
    
    Args:
        drugscreen_df (pandas.DataFrame): DataFrame containing drug screen data.
        drugdoses_file (str): Path to file containing drug doses.
        drugID_file (str): Path to file containing drug names and ChEMBL IDs.
        drugnames_file (str): Path to file containing drug names.
        column_drugname (str): Name of the column containing drug names. Default: 'drug_name'.
        column_concentration (str): Name of the column containing concentrations. Default: 'concentration'.
        column_chembl (str): Name of the column for ChEMBL IDs. Default: 'ChEMBL_ID'.
        dual_drug_mode (bool): Process both drug A and B from synergy data. Default: False.
        column_drugname_b (str): Name of column for drug B names (for dual_drug_mode). Default: None.
        column_concentration_b (str): Name of column for drug B concentrations (for dual_drug_mode). Default: None.
        directory_output (str): Directory to save output files.
        verbose (bool): Print verbose output. Default: False.
        save (bool): Save output to file. Default: False.
    
    Methods:
        get_drugdoses: Main method to get drug doses dataframe with ChEMBL IDs.
    """
    
    def __init__(self,
                 drugscreen_df: pandas.DataFrame = None,
                 drugdoses_file: str = None,
                 drugID_file: str = None,
                 drugnames_file: str = None,
                 column_drugname: str = None,
                 column_concentration: str = None,
                 column_chembl: str = None,
                 dual_drug_mode: bool = None,
                 column_drugname_b: str = None,
                 column_concentration_b: str = None,
                 directory_output: str = None,
                 verbose: bool = None,
                 save: bool = None,
                 config: dict = None):
        
        if config:
            global_config = config.get('global', {})
            paths = config.get('paths', {})
            columns = config.get('columns', {})
            options = config.get('options', {})
            self.drugscreen_df = drugscreen_df
            self.drugdoses_file = drugdoses_file
            self.drugID_file = drugID_file
            self.drugnames_file = paths.get('drug_names_file', drugnames_file)
            self.column_drugname = columns.get('drug_name', column_drugname or 'drug_name')
            self.column_concentration = columns.get('concentration', column_concentration or 'concentration')
            self.column_chembl = columns.get('chembl_id', column_chembl or 'ChEMBL_ID')
            self.dual_drug_mode = options.get('double_drug_screen', dual_drug_mode or False)
            self.column_drugname_b = columns.get('drug_name_B', column_drugname_b)
            self.column_concentration_b = columns.get('conc_B', column_concentration_b)
            self.directory_output = global_config.get('output_dir', directory_output)
            self.verbose = global_config.get('verbose', verbose or False)
            self.save = global_config.get('save', save or False)
        else:
            self.drugscreen_df = drugscreen_df
            self.drugdoses_file = drugdoses_file
            self.drugID_file = drugID_file
            self.drugnames_file = drugnames_file
            self.column_drugname = column_drugname or 'drug_name'
            self.column_concentration = column_concentration or 'concentration'
            self.column_chembl = column_chembl or 'ChEMBL_ID'
            self.dual_drug_mode = dual_drug_mode or False
            self.column_drugname_b = column_drugname_b
            self.column_concentration_b = column_concentration_b
            self.directory_output = directory_output
            self.verbose = verbose or False
            self.save = save or False
        
        self.drugscreen_df = drugscreen_df
        self.drugdoses_file = drugdoses_file
        self.drugID_file = drugID_file
        self.drugnames_file = drugnames_file
        self.column_drugname = column_drugname
        self.column_concentration = column_concentration
        self.column_chembl = column_chembl
        self.dual_drug_mode = dual_drug_mode
        self.column_drugname_b = column_drugname_b
        self.column_concentration_b = column_concentration_b
        self.directory_output = directory_output
        self.verbose = verbose
        self.save = save
        
        self.drugdoses_df = None
        self.ChEMBL_df = None
    
    def get_drugdoses(self) -> pandas.DataFrame:
        """
        Get the drug doses dataframe with ChEMBL IDs and concentrations.
        
        Returns:
            pandas.DataFrame: DataFrame with columns [drug_name, concentration, ChEMBL_ID, ChEMBL_conc].
        """
        
        # Load or create drugdoses dataframe
        if self.drugdoses_file:
            # Load the drug doses file
            self.drugdoses_df = self._load_drugdoses_file(self.drugdoses_file)
        elif self.drugscreen_df is not None:
            # Copy the drug screen dataframe
            self.drugdoses_df = self.drugscreen_df.copy()
            
            # If dual_drug_mode, expand to include both drug A and B
            if self.dual_drug_mode and self.column_drugname_b and self.column_concentration_b:
                if self.verbose:
                    print(f'\nExpanding dual-drug synergy data...')
                    print(f'Original shape: {self.drugdoses_df.shape}')
                
                # Extract Drug A side
                drug_a = self.drugdoses_df[[self.column_drugname, self.column_concentration]].copy()
                drug_a.rename(columns={
                    self.column_drugname: 'drug_name_temp',
                    self.column_concentration: 'concentration_temp'
                }, inplace=True)
                
                # Extract Drug B side
                drug_b = self.drugdoses_df[[self.column_drugname_b, self.column_concentration_b]].copy()
                drug_b.rename(columns={
                    self.column_drugname_b: 'drug_name_temp',
                    self.column_concentration_b: 'concentration_temp'
                }, inplace=True)
                
                # Combine both sides
                self.drugdoses_df = pandas.concat([drug_a, drug_b], axis=0, ignore_index=True)
                
                # Rename to standardized columns (IMPORTANT: use 'drug_name' and 'concentration' for rest of pipeline)
                self.drugdoses_df.rename(columns={
                    'drug_name_temp': 'drug_name',
                    'concentration_temp': 'concentration'
                }, inplace=True)
                
                # Drop duplicates to keep unique drug-concentration pairs
                self.drugdoses_df = self.drugdoses_df.drop_duplicates().reset_index(drop=True)
                
                # After expansion, update the column names for the rest of the pipeline
                self.column_drugname = 'drug_name'
                self.column_concentration = 'concentration'
                
                if self.verbose:
                    print(f'Expanded shape: {self.drugdoses_df.shape}')
                    print(f'Unique drugs: {self.drugdoses_df[self.column_drugname].nunique()}')
        else:
            # No drugdoses_file nor drugscreen_df, fill with 10000 nM for all drugs
            if self.drugnames_file:
                drugnames_list = self._load_drugnames_file(self.drugnames_file)
            elif self.ChEMBL_df is not None:
                drugnames_list = self.ChEMBL_df[self.column_drugname].tolist()
            else:
                raise ValueError('No drugdoses_file, drugscreen_df, or drugnames_file provided.')
            self.drugdoses_df = pandas.DataFrame({
                self.column_drugname: drugnames_list,
                self.column_concentration: [10000] * len(drugnames_list)
            })
        
        # Load or create drug ChEMBL IDs dictionary
        if self.drugID_file:
            # Load the drug ID file
            drugID_dict = self._load_drug_ChEMBL_IDs_file()
        else:
            # Make the drug ChEMBL_IDs df into a dictionary
            if self.ChEMBL_df is None:
                # Use ChEMBLIDResolver to get ChEMBL IDs
                resolver = ChEMBLIDResolver(
                    drugnames_file=self.drugnames_file,
                    drugnames_df=self.drugscreen_df if self.drugscreen_df is not None else None,
                    column_drugname=self.column_drugname,
                    column_chembl=self.column_chembl,
                    directory_output=self.directory_output,
                    verbose=self.verbose,
                    save=False
                )
                self.ChEMBL_df = resolver.resolve_chembl_ids()
            drugID_dict = self.ChEMBL_df.set_index(self.column_drugname)[self.column_chembl].to_dict()
        
        # Process doses (always use not_melted format)
        self.drugdoses_df = self._clean_drugdoses(self.drugdoses_df)
        
        # Merge ChEMBL IDs with drug doses
        self.drugdoses_df = self._merge_drugID_doses(self.drugdoses_df, drugID_dict)
        
        if self.save and self.directory_output:
            save_file(self.drugdoses_df, self.directory_output, 'drug_ChEMBL_doses.csv', index=False)
        
        return self.drugdoses_df
    
    def _merge_drugID_doses(self, drugdoses, drugID_dict) -> pandas.DataFrame:
        """
        Merge ChEMBL IDs with drug doses.
        
        Args:
            drugdoses (pandas.DataFrame): DataFrame with drug doses.
            drugID_dict (dict): Dictionary mapping drug names to ChEMBL IDs.
            
        Returns:
            pandas.DataFrame: Merged DataFrame with ChEMBL_conc column.
        """
        # Map drug names to ChEMBL IDs using the drug ID dictionary
        drugdoses_IDs = self._map_drugID(drugdoses, drugID_dict)
        # Merge ChEMBL IDs with drug doses
        drugdoses_IDs['ChEMBL_conc'] = (
            drugdoses_IDs[self.column_chembl] + '_' + 
            drugdoses_IDs[self.column_concentration].astype(str) + '_nM'
        )
        # Drop rows with missing ChEMBL IDs
        drugdoses_IDs = drugdoses_IDs.dropna(subset=[self.column_chembl]).reset_index(drop=True)
        
        if self.verbose:
            print('Shape of Drug doses dataframe after dropping missing ChEMBL IDs:', drugdoses_IDs.shape)
            print('Columns in drugdoses_merged:', drugdoses_IDs.columns)
            print('Head of drugdoses_merged:', drugdoses_IDs.head(3), '\n')
            print('ChEMBL IDs merged with drug doses successfully\n')
        
        return drugdoses_IDs
    
    def _map_drugID(self, drugdoses, drugID_dict):
        """
        Map drug names to ChEMBL IDs using the drug ID dictionary.
        
        Args:
            drugdoses (pandas.DataFrame): DataFrame with drug doses.
            drugID_dict (dict): Dictionary mapping drug names to ChEMBL IDs.
            
        Returns:
            pandas.DataFrame: DataFrame with ChEMBL IDs mapped.
        """
        drugdoses_IDs = drugdoses.copy()
        drugdoses_IDs[self.column_chembl] = drugdoses_IDs[self.column_drugname].map(drugID_dict)
        
        if self.verbose:
            # Find drug names not mapped to ChEMBL IDs
            print('Columns in drugdoses_IDs:', drugdoses_IDs.columns)
            print('Head of drugdoses_IDs:', drugdoses_IDs.head(3), '\n')
            print('Drug names not mapped to ChEMBL IDs:', 
                  drugdoses_IDs[drugdoses_IDs[self.column_chembl].isnull()][self.column_drugname].unique())
        
        return drugdoses_IDs
    
    def _melt_drugdoses(self, drugdoses) -> pandas.DataFrame:
        """
        Melt drug doses dataframe to long format (for files with multiple dose columns).
        
        Args:
            drugdoses (pandas.DataFrame): DataFrame with drug doses.
            
        Returns:
            pandas.DataFrame: Melted DataFrame in long format.
        """
        if self.verbose:
            print('\nMelting the drug doses dataframe ...\n')
            print('Drug doses dataframe before melting:', drugdoses.shape, '\n')
        
        # Melt the dataframe for each concentration column to be in a separate row
        # Concentration columns are all except the first column (drug names)
        drugdoses_melted = drugdoses.melt(
            id_vars=[drugdoses.columns[0]],  # drug names
            value_vars=drugdoses.columns[1:],  # concentration columns
            var_name='dose_type',  # this can be any placeholder name
            value_name=self.column_concentration
        )
        
        # Keep only the necessary columns
        drugdoses_melted = drugdoses_melted[[drugdoses.columns[0], self.column_concentration]]
        drugdoses_melted = drugdoses_melted.rename(columns={drugdoses.columns[0]: self.column_drugname})
        
        if self.verbose:
            print('Drug doses dataframe after melting:', drugdoses_melted.shape, drugdoses_melted.columns, '\n')
            print('Head of drug doses dataframe after melting:', drugdoses_melted.head(3), '\n')
        
        # Clean the melted dataframe
        drugdoses_melted = self._clean_drugdoses(drugdoses_melted)
        
        return drugdoses_melted
    
    def _clean_drugdoses(self, drugdoses) -> pandas.DataFrame:
        """
        Clean the drug doses dataframe (remove duplicates and zero concentrations).
        
        Args:
            drugdoses (pandas.DataFrame): DataFrame with drug doses.
            
        Returns:
            pandas.DataFrame: Cleaned DataFrame.
        """
        if self.verbose:
            print('Cleaning the drug doses dataframe ...')
            print('Dropping duplicates doses and doses in 0 concentration...')
            print('Columns in Drug doses dataframe before dropping duplicates:', drugdoses.columns)
            print('Shape of Drug doses dataframe before dropping duplicates:', drugdoses.shape)
        
        if self.column_concentration is None:
            self.column_concentration = 'concentration'
        
        drugdoses = drugdoses[[self.column_drugname, self.column_concentration]].drop_duplicates().reset_index(drop=True)
        drugdoses = drugdoses[drugdoses[self.column_concentration] != 0]
        
        return drugdoses
    
    def _load_drug_ChEMBL_IDs_file(self) -> dict:
        """
        Load the drug ID file (dictionary mapping drug names to ChEMBL IDs).
        
        Returns:
            dict: Dictionary mapping drug names to ChEMBL IDs.
        """
        if self.drugID_file is not None:
            if self.verbose:
                print('Loading the drug ID file ...')
                print('Drug ID file:', self.drugID_file)
            drugID = pandas.read_csv(self.drugID_file)
            # Rename columns to match expected names
            drugID = drugID.rename(columns={
                drugID.columns[0]: self.column_drugname,
                drugID.columns[1]: self.column_chembl
            })
            drugID_dict = drugID.set_index(self.column_drugname)[self.column_chembl].to_dict()
        else:
            # Use ChEMBLIDResolver to get ChEMBL IDs
            resolver = ChEMBLIDResolver(
                drugnames_file=self.drugnames_file,
                drugnames_df=self.drugscreen_df,
                column_drugname=self.column_drugname,
                column_chembl=self.column_chembl,
                directory_output=self.directory_output,
                verbose=self.verbose,
                save=False
            )
            drugID = resolver.resolve_chembl_ids()
            drugID_dict = drugID.set_index(self.column_drugname)[self.column_chembl].to_dict()
        
        return drugID_dict
    
    def _load_drugdoses_file(self, drugdoses_file) -> pandas.DataFrame:
        """
        Load and read the drug doses file (CSV file).
        
        Args:
            drugdoses_file (str): Path to drug doses file.
            
        Returns:
            pandas.DataFrame: DataFrame with drug doses.
        """
        drugdoses = pandas.read_csv(drugdoses_file)
        
        if self.verbose:
            print('Drug doses file:', drugdoses_file)
            print('Drug doses dataframe:', drugdoses.head(), '\n')
        
        return drugdoses
    
    def _load_drugnames_file(self, drugnames_file) -> list:
        """
        Load a file containing drug names (CSV file).
        
        Args:
            drugnames_file (str): Path to drug names file.
            
        Returns:
            list: List of drug names.
        """
        # Use ChEMBLIDResolver to load drug names
        resolver = ChEMBLIDResolver(
            drugnames_file=drugnames_file,
            column_drugname=self.column_drugname,
            verbose=self.verbose
        )
        return resolver._load_drugnames_file(drugnames_file)
