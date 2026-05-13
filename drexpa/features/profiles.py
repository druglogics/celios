# Drug profiles creation module

import pandas
from .helpers import save_file


class ProfileBuilder:
    """
    Drug profiles builder class.
    
    This class handles the creation of drug profiles by grouping drugs with the same node targets.
    
    Args:
        node_targets_df (pandas.DataFrame): DataFrame containing node targets.
        node_targets_file (str): Path to file containing node targets.
        profile_values (str): Column name to use for profile grouping ('ChEMBL_conc' or 'drug_name'). Default: 'ChEMBL_conc'.
        directory_output (str): Directory to save output files.
        verbose (bool): Print verbose output. Default: False.
        save (bool): Save output to file. Default: False.
    
    Methods:
        get_drugprofiles: Main method to get drug profiles dataframe.
    """
    
    def __init__(self,
                 node_targets_df: pandas.DataFrame = None,
                 node_targets_file: str = None,
                 profile_values: str = 'ChEMBL_conc',
                 directory_output: str = None,
                 verbose: bool = False,
                 save: bool = False):
        
        self.node_targets_df = node_targets_df
        self.node_targets_file = node_targets_file
        self.profile_values = profile_values
        self.directory_output = directory_output
        self.verbose = verbose
        self.save = save
        
        self.drugprofiles_df = None
        self.PD_profiles_dict = None
        self.drugprofiles_output_file = 'drug_profiles.csv'
        self.drugprofiles_dict_output_file = 'drug_profiles_dict.csv'
    
    def get_drugprofiles(self, 
                        node_targets_df: pandas.DataFrame = None) -> pandas.DataFrame:
        """
        Get the pipeline drug profiles by grouping drugs with the same node targets.
        
        Args:
            node_targets_df (pandas.DataFrame): Override node_targets_df if provided.
            
        Returns:
            pandas.DataFrame: DataFrame with PD_profile column added.
        """
        if node_targets_df is not None:
            self.node_targets_df = node_targets_df
        
        if self.verbose:
            print('\nCreating pipeline drug profiles ...')
        
        if self.node_targets_df is None:
            if self.node_targets_file:
                if self.verbose:
                    print('Reading the node targets dataframe from file ...')
                    print('Node targets file:', self.node_targets_file, '\n')
                self.node_targets_df = pandas.read_csv(self.node_targets_file)
            else:
                raise ValueError('No node_targets_df or node_targets_file provided. Please provide node targets data.')
        
        # Make the pipeline drug profiles
        self.drugprofiles_df, self.PD_profiles_dict = self._make_drugprofiles()
        
        if self.verbose:
            print('Pipeline drug profiles created successfully.')
            print('Columns in the updated node_targets dataframe:', self.drugprofiles_df.columns, '\n')
        
        if self.save and self.directory_output:
            # Save the updated node_targets_df to a new file
            save_file(self.drugprofiles_df, self.directory_output, self.drugprofiles_output_file, index=False)
            drugprofiles_dict = pandas.DataFrame(
                self.PD_profiles_dict.items(), 
                columns=['PD_profile', self.profile_values]
            )
            save_file(drugprofiles_dict, self.directory_output, self.drugprofiles_dict_output_file, index=False)
        
        return self.drugprofiles_df
    
    def _make_drugprofiles(self) -> tuple:
        """
        Make the pipeline drug profiles by grouping drugs with the same node targets.
        
        Returns:
            tuple: (drugprofiles_df, PD_profiles_dict)
        """
        # Clean the node targets dataframe
        self.node_targets_df = self._clean_node_targets_df(self.node_targets_df)
        
        self.PD_profiles_dict = {}
        self.node_targets_to_pd = {}  # NEW: Map from node_targets tuple to PD_profile
        
        # Group the drugs by the node targets and create the pipeline drug profiles
        for node_targets in self.node_targets_df['node_targets'].apply(tuple).unique():
            # Filter the DataFrame to only include rows with the current node_targets
            subset = self.node_targets_df[self.node_targets_df['node_targets'].apply(tuple) == node_targets]
            
            # Generate a new name for the group
            group_name = f"PD_{len(self.PD_profiles_dict) + 1:02d}"
            
            # Add the group to the dictionary, removing duplicate drug names while preserving order
            drug_list = list(subset[self.profile_values])
            unique_drugs = []
            seen = set()
            for drug in drug_list:
                if drug not in seen:
                    unique_drugs.append(drug)
                    seen.add(drug)
            self.PD_profiles_dict[group_name] = unique_drugs
            
            # NEW: Map this node_targets combination to its PD_profile
            self.node_targets_to_pd[node_targets] = group_name
        
        if self.verbose:
            print('Number of pipeline drug profiles:', len(self.PD_profiles_dict))
            print('Pipeline drug profiles created successfully.\n')
        
        # Add the pipeline drug profiles as a new column to the node targets dataframe
        self.drugprofiles_df = self._add_PD_profiles_df(self.node_targets_df, self.PD_profiles_dict)
        
        return self.drugprofiles_df, self.PD_profiles_dict
    
    def _clean_node_targets_df(self, node_targets_df) -> pandas.DataFrame:
        """
        Clean the node targets dataframe.
        
        Args:
            node_targets_df (pandas.DataFrame): DataFrame with node targets.
            
        Returns:
            pandas.DataFrame: Cleaned DataFrame.
        """
        # Check the column node_targets for 'No matching targets' values and fill them with NaN
        node_targets_df['node_targets'] = node_targets_df['node_targets'].apply(
            lambda x: x if x != 'No matching targets' else None
        )
        
        if self.verbose:
            print('Checking for empty values in the node_targets column ...')
            print('Number of drugs with no targets in the model:', node_targets_df['node_targets'].isnull().sum())
        
        # Clean node_targets to remove NaN values
        node_targets_df = node_targets_df.dropna(subset=['node_targets'])
        
        if self.verbose:
            print('Removed drugs with no targets from the dataframe successfully.')
            print('Shape of the updated node_targets dataframe:', node_targets_df.shape, '\n')
        
        return node_targets_df
    
    def _add_PD_profiles_df(self, node_targets_df, pipeline_ID_dict) -> pandas.DataFrame:
        """
        Add the pipeline drug profiles as a new column in the node_targets_df.
        Match each row based on its node_targets tuple (not just drug name).
        
        Args:
            node_targets_df (pandas.DataFrame): DataFrame with node targets.
            pipeline_ID_dict (dict): Dictionary of pipeline profiles.
            
        Returns:
            pandas.DataFrame: DataFrame with PD_profile column added.
        """
        # Add the pipeline drug profiles to the node targets dataframe
        # Use the node_targets_to_pd mapping created during _make_drugprofiles
        drugprofiles_df = node_targets_df.copy()
        
        # Map each row's node_targets tuple to its PD_profile
        drugprofiles_df['PD_profile'] = drugprofiles_df['node_targets'].apply(
            lambda targets: self.node_targets_to_pd.get(targets if isinstance(targets, tuple) else tuple(targets), None)
        )
        
        if self.verbose:
            print('Adding pipeline drug profiles to the node_targets dataframe ...')
            print('Columns in the updated node_targets dataframe:', drugprofiles_df.columns)
            print('Shape of the updated node_targets dataframe:', drugprofiles_df.shape, '\n')
        
        return drugprofiles_df
