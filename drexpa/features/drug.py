# Drug processing module

import pandas
import ast
from drexpa.features.helpers import save_file
import drexpa.features.target_checker as target_checker
from chembl_webresource_client.new_client import new_client

#/////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////
class Drug:
    """
    Drug class.

        This class contains all the methods to process the drug data.
        Starting from a list of drug names, or the data from a drug screen with drug concentrations,
        the class will create the drug profiles to be used in the pipeline.
        
    Args:
        profile_values (str): type of identifiers to create the drug profiles. Default is 'ChEMBL_conc'. *

            * Profile values can be 'ChEMBL_conc' or 'drug_name'.
                It is recommended to use 'ChEMBL_conc' when working with different drug concentrations.
                If no concentration data is available, use 'drug_name' instead.
    Methods:
        get_drugprofiles: get the pipeline drug profiles dataframe [drug_name, ChEMBL_conc, targets, node_targets, PD_profile].
        get_node_targets: get the node targets dataframe [drug_name, ChEMBL_conc, targets, node_targets].
        get_chembl_targets: get the ChEMBL targets dataframe [drug_name, targets, ChEMBL_ID, ChEMBL_conc].
        get_drugdoses: get the drug doses dataframe [drug_name, concentration, ChEMBL_ID, ChEMBL_conc].
    
    """

    def __init__(self,
                drugscreen_df: pandas.DataFrame = None,
                drugnames_file: str = None,         # text file with drug names separated by new lines. *
                node_targets_df: pandas.DataFrame = None,
                node_targets_file: str = None,
                drugID_file: str = None,        # text file with drug names and ChEMBL IDs separated by new lines.
                drugdoses_file: str = None,     # text file with drug names and concentrations separated by new lines.
                doses_type: str = 'not_melted',     # or 'not_melted' **
                node_dict: dict = None,
                node_dict_file: str = None,                
                profile_values: str = 'ChEMBL_conc',  # or 'drug_name'
                db_file: str = None,
                ic50_value: float = None,
                directory_output: str = None,
                verbose: bool = False
                ):

        self.profile_values = profile_values
        self.directory_output = directory_output
        self.drugprofiles_output_file = 'drug_profiles.csv'
        self.drugprofiles_dict_output_file = 'drug_profiles_dict.csv'
        self.verbose = verbose

        self.node_targets_df = node_targets_df    # dataframe containing the node targets.
        self.node_targets_file = node_targets_file
        self.drugprofiles_df = None       # dataframe containing the pipeline drug profiles (drugs with the same node targets).
        self.PD_profiles_dict = None      # pipeline drug profiles dictionary.

        self.node_dict = node_dict
        self.node_dict_file = node_dict_file
        self.chembl_targets_df = None     # dataframe containing drug ChEMBL IDs and targets.

        self.db_file = db_file    # database file containing drug-target interactions.
        self.ChEMBL_df = None    # dataframe containing drug names and ChEMBL IDs.
        self.column_chembl='ChEMBL_ID'
        self.column_drugname='drug_name'      # or 'anchor_name', 'library_name'
        self.column_concentration='concentration'
        self.column_target='targets'    # or 'combined_targets'
        self.ic50_value = ic50_value
        self.merge_on = 'ChEMBL_conc' # Or ChEMBL_ID

        self.drugdoses_df = None
        self.drugID_file = drugID_file
        self.drugdoses_file = drugdoses_file
        self.doses_type = doses_type

        self.drugscreen_df = drugscreen_df
        self.drugnames_file = drugnames_file

        # * Header should match the drugname_column value. If not, it will be renamed to this value.
        # ** Use 'melted' doses_type, when the drugdose file has multiple doses columns. [Example: oncologics drug doses data]
        #       Otherwise, use 'not_melted', when the dataframe has a unique column with all doses.

    #/////////////////////////////////////////////////////////////////////////////////////////
    # BRANCH: drugprofiles
    #/////////////////////////////////////////////////////////////////////////////////////////

    def get_drugprofiles(self,) -> pandas.DataFrame:
        '''
        Function: get the pipeline drug profiles for each drug in the node targets dataframe.
        
        '''

        if self.verbose:
            print('\nCreating pipeline drug profiles ...')

        if self.node_targets_df is None:
            # Get the node targets dataframe from branch drug_node_targets
            self.node_targets_df = self.get_node_targets()

        if self.node_targets_file:
            if self.verbose:
                print('Reading the node targets dataframe from file ...')
                print('Node targets file:', self.node_targets_file, '\n')
            self.node_targets_df = pandas.read_csv(self.node_targets_file)

        # Make the pipeline drug profiles
        self.drugprofiles_df, self.PD_profiles_dic = self._make_drugprofiles()

        if self.verbose:
            print('Pipeline drug profiles created successfully.')
            print('Columns in the updated node_targets dataframe:', self.drugprofiles_df.columns, '\n')

        if self.directory_output:
            # Save the updated node_targets_df to a new file
            save_file(self.drugprofiles_df, self.directory_output, self.drugprofiles_output_file, index=False)
            drugprofiles_dict = pandas.DataFrame(self.PD_profiles_dic.items(), columns=['PD_profile', self.profile_values])
            save_file(drugprofiles_dict, self.directory_output, self.drugprofiles_dict_output_file, index=False)

        return self.drugprofiles_df
                        
    #//////////////////////////////////////////////////////////
    def _make_drugprofiles(self) -> pandas.DataFrame:
        '''
        Function: make the pipeline drug profiles by grouping drugs with the same node targets.
        '''


        # Clean the node targets dataframe
        self.node_targets_df = self._clean_node_targets_df(self.node_targets_df) 

        self.PD_profiles_dict = {}

        # Group the drugs by the node targets and create the pipeline drug profiles
        for node_targets in self.node_targets_df['node_targets'].apply(tuple).unique():
            # Filter the DataFrame to only include rows with the current node_targets
            subset = self.node_targets_df[self.node_targets_df['node_targets'].apply(tuple) == node_targets]

            # Generate a new name for the group
            group_name = f"PD_{len(self.PD_profiles_dict) + 1:02d}"

            # Add the group to the dictionary
            self.PD_profiles_dict[group_name] = list(subset[self.profile_values])

        if self.verbose:
            print('Number of pipeline drug profiles:', len(self.PD_profiles_dict))
            print('Pipeline drug profiles created successfully.\n')
                
        # Add the pipeline drug profiles as a new column to the node targets dataframe
        self.drugprofiles_df = self._add_PD_profiles_df(self.node_targets_df, self.PD_profiles_dict)
        
        return self.drugprofiles_df, self.PD_profiles_dict
    
    #//////////////////////////////////////////////////////////
    def _clean_node_targets_df(self, node_targets_df) -> pandas.DataFrame:
        '''
        Function: clean the node targets dataframe.
        '''
        # Check the column node_targets for 'No matching targets' values and fill them with NaN
        node_targets_df['node_targets'] = node_targets_df['node_targets'].apply(lambda x: x if x != 'No matching targets' else None)
        
        if self.verbose:
            print('Checking for empty values in the node_targets column ...')
            print('Number of drugs with no targets in the model:', node_targets_df['node_targets'].isnull().sum())

        # Clean node_targets to remove NaN values
        node_targets_df = node_targets_df.dropna(subset=['node_targets'])

        if self.verbose:
            print('Removed drugs with no targets from the dataframe successfully.')
            print('Shape of the updated node_targets dataframe:', node_targets_df.shape, '\n')
        
        return node_targets_df
    
    #//////////////////////////////////////////////////////////
    def _add_PD_profiles_df(self, node_targets_df, pipeline_ID_dict) -> pandas.DataFrame:
        '''
        Function: add the pipeline drug profiles as a new column in the node_targets_df.
        '''

        # Match the drugs to the pipeline drug profiles
        pipeline_ID_inverted = self._match_drugs_to_profiles(pipeline_ID_dict)

        # Add the pipeline drug profiles to the node targets dataframe
        drugprofiles_df = node_targets_df.copy()
        drugprofiles_df['PD_profile'] = drugprofiles_df[self.profile_values].apply(lambda drug_id: pipeline_ID_inverted.get(drug_id, None))

        if self.verbose:
            print('Adding pipeline drug profiles to the node_targets dataframe ...')
            print('Columns in the updated node_targets dataframe:', self.node_targets_df.columns,)
            print('Shape of the updated node_targets dataframe:', self.node_targets_df.shape, '\n')

        return drugprofiles_df
    
    #//////////////////////////////////////////////////////////
    def _match_drugs_to_profiles(self, drugprofiles_dict ) -> dict:
        '''
        Function: make a dictionary of drugs to pipeline drug profiles.
        '''
        drugprofiles_inverted = {drug_id: group_name for group_name, drug_ids in drugprofiles_dict.items() for drug_id in drug_ids}

        if self.verbose:
            print('Matching drugs to pipeline drug profiles ...')
            print('Length of pipeline drug profiles dictionary:', len(drugprofiles_inverted), '\n')

        return drugprofiles_inverted


    #/////////////////////////////////////////////////////////////////////////////////////////
    # BRANCH: drug_node_targets
    #/////////////////////////////////////////////////////////////////////////////////////////

    def get_node_targets(self,
                        column_target: str = None,
                        column_drugname: str = None
                        ) -> pandas.DataFrame:
        '''
        Function: get the node_targets dataframe.
        '''
        if column_target:
            self.column_target = column_target

        if column_drugname:
            self.column_drugname = column_drugname

        if self.verbose:
            print('\nMaking the node_targets dataframe ...')
            print('\nMapping drug HGNC targets to node names...')
            print(f'Using the target column {self.column_target} to map the targets to node names.\n')

        if self.node_dict_file:
            # Read the node dictionary file
            self.node_dict = self._read_node_dict_file(self.node_dict_file)

        if self.chembl_targets_df is None:
            # Get the ChEMBL targets dataframe from branch drug_ChEMBL_targets
            self.chembl_targets_df = self.get_chembl_targets()

        node_targets_df = self._map_target_nodes(
            self.chembl_targets_df,
            self.node_dict
        )

        if self.directory_output:
            save_file(node_targets_df, self.directory_output, 'drug_node_targets.csv', index=False)

        return node_targets_df

    #//////////////////////////////////////////////////////////
    def _read_node_dict_file(self, node_dict_file) -> dict:
        '''
        Function: read the node dictionary file.
        '''
        if self.verbose:
            print('\nReading the node-gene symbols dictionary from file ...')
            print('Node dictionary file:', node_dict_file)
        node_dict = pandas.read_csv(node_dict_file)
        node_dict = {k: ast.literal_eval(v) for k, v in zip(node_dict.node_name, node_dict.HGNC_symbol)}
        
        if self.verbose:
            print('Node dictionary:', node_dict, '\n')

        return node_dict
    #//////////////////////////////////////////////////////////
    def _map_target_nodes(self,
                        chembl_targets_df: pandas.DataFrame,
                        node_dict: dict,
                        ) -> pandas.DataFrame:
        '''
        Function: map each target to its corresponding HGNC symbol, for each drug ChEMBL ID.

        Args:
        - chembl_targets_df: the ChEMBL targets dataframe.
        - node_dict: the HGNC node dictionary. Key: node name, Value: list of HGNC symbols.
        '''

        # Clean the ChEMBL targets dataframe
        chembl_targets_df = self._clean_ChEMBLtargets(chembl_targets_df)

        if self.verbose:
            print('Mapping HGNC targets to node names ...\n')

        def map_targets(column_target):    
            matching_nodenames = []
            for node_key, hgnc_values in node_dict.items():
                # print(f'Node name: {node_key}, HGNC value: {hgnc_values}')
                for value in hgnc_values:
                    if value in column_target:
                        matching_nodenames.append(node_key)
                    # print(f'Node {node_key} mapped to {value} HGNC symbol')
            node_targets = list(set(matching_nodenames)) # Get unique node names
            return node_targets if node_targets else 'No matching targets'
    
        chembl_targets_df['node_targets'] = chembl_targets_df[self.column_target].apply(map_targets)

        if self.verbose:
            print('HGNC targets mapped to node names successfully.\n')
            print('Node targets dataframe:', chembl_targets_df.head(5))
            print('Shape of node_targets dataframe:', chembl_targets_df.shape)   

        return chembl_targets_df
    
    #//////////////////////////////////////////////////////////
    def _clean_ChEMBLtargets(self,
                            chembl_targets_df: pandas.DataFrame,
                            column_drugname: str = None,
                            column_target: str = None
                            ) -> pandas.DataFrame:
        '''
        Function: clean the ChEMBL targets dataframe.
        '''
        if self.verbose:
            print('\nCleaning ChEMBLtargets dataframe ...\n')
            print('Keeping only the columns:', self.column_drugname, self.column_target)
        
        try:
            chembl_targets_df = chembl_targets_df[[self.column_drugname, self.column_target, 'ChEMBL_conc']]
        except KeyError:
            chembl_targets_df = chembl_targets_df[[self.column_drugname, self.column_target]] # Take this part out when testing drug screen with no doses
        
        if chembl_targets_df[self.column_target].isnull().sum() > 0:
            chembl_targets_df.loc[:, self.column_target] = chembl_targets_df[self.column_target].replace('', None).fillna('not found')
            
        if self.verbose:
            # If there are 'not found' values in the targets column
            if chembl_targets_df[self.column_target].str.contains('not found').any():
                print('Drugs with not found values in targets column.')
                print('Number of drugs with not found targets:', chembl_targets_df[self.column_target].str.contains('not found').sum())
            else:
                print('No drugs with missing values in the targets column.\n')

        chembl_targets_df = chembl_targets_df[chembl_targets_df[self.column_target] != 'not found']

        if self.verbose:
            print('Shape of ChEMBLtargets dataframe:', chembl_targets_df.shape)

        return chembl_targets_df


    #/////////////////////////////////////////////////////////////////////////////////////////
    # BRANCH: drug_ChEMBL_targets
    #/////////////////////////////////////////////////////////////////////////////////////////

    def get_chembl_targets(self,
                            merge_on: str = None,
                            # drugdoses_screen: str = None, # or 'main_screen'
                            ) -> pandas.DataFrame:
        '''
        Function: get the ChEMBL targets dataframe.
        '''
        if merge_on:
            self.merge_on = merge_on

        if self.verbose:
            print('\nMaking the ChEMBL targets dataframe ...')
            print('Merging on:', self.merge_on, '\n')
        
        # if self.ChEMBL_df is None:
        #     # Get the ChEMBL IDs dataframe from branch drug_ChEMBL_IDs
        #     if self.verbose:
        #         print('Making the ChEMBL IDs dataframe ...\n')
        #     self.ChEMBL_df = self._process_drugnames_to_ChEMBL()

        # if drugdoses_screen is None:
            # Get the drug doses dataframe from branch drug_doses
        drugdoses_screen = self.get_drugdoses()

        # Find targets for the ChEMBL IDs
        target_query = self._find_targets(drugdoses_screen)

        # Merge the ChEMBL IDs dataframe with the target query dataframe
        chembl_targets_df = self._merge_ChEMBL_targetsquery(target_query, drugdoses_screen)

        # Check the ChEMBL targets dataframe for found and missing targets
        chembl_targets_df = self._check_ChEMBL_targets(chembl_targets_df)

        if self.directory_output:
            save_file(chembl_targets_df, self.directory_output, 'drug_ChEMBL_targets.csv', index=False)

        if self.verbose:
            print('ChEMBL targets dataframe created successfully.\n')

        return chembl_targets_df

    #//////////////////////////////////////////////////////////
    def _check_ChEMBL_targets(self, chembl_targets_df) -> pandas.DataFrame:
        '''
        Function: check the ChEMBL targets dataframe for found and missing targets.
        '''
        # Fill drugs with no targets found (missing values) with NAN
        chembl_targets_df['targets'] = chembl_targets_df['targets'].replace('', None)

        # Count the number of drugs with found (No. of targets > 0) and missing targets (No. of targets = None)
        missing_targets = chembl_targets_df['targets'].isnull().sum()
        found_targets = len(chembl_targets_df) - missing_targets

        if self.verbose:
            print(f'Number of drugs with missing targets: {missing_targets}')
            print(f'Number of drugs with found targets: {found_targets}\n')

        # Get the list of drugs with missing targets
        self._get_missing_ChEMBL_targets(chembl_targets_df)

        return chembl_targets_df

    #//////////////////////////////////////////////////////////
    def _get_missing_ChEMBL_targets(self, chembl_targets_df) -> list:
        '''
        Function: get the list of drugs with missing targets.
        
        parameters:
            - chembl_targets_df: dataframe containing drug names, ChEMBL IDs, drug IDs, HGNC symbols, and IC50 values.
    '''
        # Get the list of drugs with missing targets
        missing_targets = chembl_targets_df[chembl_targets_df['targets'].isnull()]
        # Create a list of drugs with missing targets
        missing_targets = missing_targets[self.column_drugname].tolist()
        # Drop the duplicates in the list
        missing_targets = list(set(missing_targets))
        
        if self.verbose:
            print('List of drugs with missing targets:', missing_targets, '\n')
        
        return missing_targets

    #//////////////////////////////////////////////////////////
    def _merge_ChEMBL_targetsquery(self, target_query, drugdoses_screen) -> pandas.DataFrame:
        '''
        Function: merge the ChEMBL IDs dataframe with the target query dataframe.
        
        Args:
            - target_query: dataframe containing drug IDs and HGNC symbols.
        
        '''
        target_query = self._split_drugIDconcentration(target_query)

        # Merge the ChEMBL IDs dataframe with the target query dataframe on ChEMBL_ID column
        if self.verbose:
            print('\nMerging ChEMBL targets dataframe with target query dataframe ...\n')
        chembl_targets_df = pandas.merge(drugdoses_screen, target_query, on=self.merge_on, how='left')

        if self.verbose:
            print('Length of merged dataframe:', len(chembl_targets_df))
            print('Columns in the merged dataframe:', chembl_targets_df.columns, '\n')
            print('Head of the merged dataframe:', chembl_targets_df.head(), '\n')
            print('Dataframes merged successfully\n')
        
        # Structure the ChEMBL targets dataframe
        if self.verbose:
            print('Counting the number of targets for each ChEMBL ID ...')
            print('Creating a new column with the number of targets ...')
        chembl_targets_df['no_targets'] = chembl_targets_df['targets'].apply(self._count_targets).astype('Int64')

        return chembl_targets_df

    #//////////////////////////////////////////////////////////
    def _split_drugIDconcentration(self, target_query) -> pandas.DataFrame:
        '''
        Function: rename the columns in the target query dataframe to get drugID and IC50 columns.

        parameters:
            - target_query: dataframe containing drug IDs and HGNC symbols.
        '''
        if self.verbose:
            print('\nSplitting and renaming columns in target query dataframe ...\n')
        # Add a new column with the ChEMBL_ID
        target_query['ChEMBL_ID'] = target_query['drugID'].str.split('_').str[0]
        # Add a new column with the IC50 Concentration
        target_query['IC50'] = target_query['drugID'].str.split('_').str[1]
        # Rename the columns
        target_query = target_query.rename(columns={'drugID': 'ChEMBL_conc'})

        if self.verbose:
            print('Columns in target query dataframe:', target_query.columns)
            
        return target_query

    #//////////////////////////////////////////////////////////
    def _count_targets(self, targets) -> int:
        '''
        Function: count the number of targets for each ChEMBL ID and create a new column with the number of targets.
        '''
        if pandas.isna(targets):
            return None
        return len(targets.split(','))

    #//////////////////////////////////////////////////////////
    def _find_targets(self, drugdoses_screen) -> pandas.DataFrame:
        '''
        Function: find targets for the ChEMBL IDs using the DrugTargetInteractionDB database.
        Uses the target_checker.py script to find targets for a list of drugs.

        Args:
        - db_file: the DrugTargetInteractionDB database file, containing containing drug-target interactions.
        - drugdoses_screen: dataframe containing drug names, ChEMBL IDs, and IC50 values.
        - chembl_column: the column in the ChEMBL targets dataframe containing the ChEMBL IDs.
        - ic50_column: column name of the IC50 values in the ChEMBL dataframe. Default: 'concentration'.
        - ic50_limit: upper limit for IC50 value for binding affinity, when only one concentration is available. Default: None.
        '''
        target_holder = []

        # Iterate over the ChEMBL IDs dataframe
        if self.verbose:
            print('\nFinding targets for ChEMBL IDs ...\n')
        for index, row in drugdoses_screen.iterrows():        
            chembl_id = row[self.column_chembl]
            if self.ic50_value is not None:
                target_frames = target_checker.targetProfileBA(
                    self.db_file,
                    'drugpanel_supplementary.txt',
                    id_list=[chembl_id],
                    ic50_limit = self.ic50_value)
            else:
                target_frames = target_checker.targetProfileBA(
                    self.db_file,
                    'drugpanel_supplementary.txt',
                    id_list=[chembl_id],
                    ic50_limit = row[self.column_concentration])

            for frame in target_frames:
                if len(frame)!=0:
                    target_holder.append(frame)
                else:
                    errorFrame = pandas.DataFrame(
                        [[chembl_id, None]],
                        columns=['drugID', 'HGNC'])
                    target_holder.append(errorFrame)
        
        # Concatenate the target frames        
        target_query = pandas.concat(target_holder)

        if self.verbose:
            print('\nTargets found successfully\n')
            print('Length of targetquery_df:', len(target_query))
            print('Number of drugs with no found targets:', target_query.isnull().sum().sum(), '\n')

        target_query = self._aggregate_targets(target_query)

        return target_query

    #//////////////////////////////////////////////////////////
    def _aggregate_targets(self, target_query)-> pandas.DataFrame:
        '''
        Function: aggregate the targets in the target query dataframe.
        '''
        target_query = target_query.rename(columns={'HGNC': 'targets'})
        target_query = target_query.groupby(
            ['drugID'], as_index=False).agg({
                'drugID': 'first',
                'targets': lambda x: ', '.join(filter(None, x))
            }
        )
        return target_query


    #/////////////////////////////////////////////////////////////////////////////////////////
    # BRANCH: drug_ChEMBL_IDs
    #/////////////////////////////////////////////////////////////////////////////////////////

    def _process_drugnames_to_ChEMBL(self,) -> pandas.DataFrame:
        '''
        Function: process drug names from a file to ChEMBL IDs using ChEMBL web services API.
        '''

        if self.verbose:
            print('\nProcessing drug names to ChEMBL IDs ...')
        
        if self.drugnames_file:
            if self.verbose:
                print('Using drug names file:', self.drugnames_file, '\n')
            drugnames_list = self._load_drugnames_file(self.drugnames_file)
        elif self.drugscreen_df is not None:
            if self.verbose:
                print('Using drug screen dataframe ...')
            drugscreen_df = self.drugscreen_df.drop_duplicates(subset=[self.column_drugname])
            drugnames_list = drugscreen_df[self.column_drugname].tolist()
        else:
            # Raise error
            raise ValueError('No drug names file or drug screen dataframe provided.')
        
        # Make the ChEMBL IDs dataframe
        ChEMBL_df = self._make_df_ChEMBL_IDs(drugnames_list)

        return ChEMBL_df

    #//////////////////////////////////////////////////////////
    def _make_df_ChEMBL_IDs(self, drugnames_list) -> pandas.DataFrame:
        '''
        Function: make a dataframe with drug names and ChEMBL IDs after mapping drug names to ChEMBL IDs using ChEMBL web services API.
        '''

        # Call the _count_ChEMBL_IDs function
        IDs = self._count_ChEMBL_IDs(drugnames_list)

        ChEMBL_df = pandas.DataFrame({
            self.column_drugname: drugnames_list,
            self.column_chembl: IDs
        })

        if self.verbose:
            print('Processed drug names to ChEMBL symbols.')
            print('Columns in the ChEMBL IDs dataframe:', ChEMBL_df.columns)
        
        # Clean NaN from ChEMBL dataframe
        ChEMBL_df = self._clean_nan_ChEMBL_IDs(ChEMBL_df)

        return ChEMBL_df

            
    #//////////////////////////////////////////////////////////
    def _clean_nan_ChEMBL_IDs(self, ChEMBL_df):
        '''
        Function: clean NaN from ChEMBL dataframe after mapping drug names to ChEMBL IDs using ChEMBL web services API.
        '''

        ChEMBL_df = ChEMBL_df.dropna().reset_index(drop=True)

        if self.verbose:
            print('\nCleaning NaN from ChEMBL dataframe ...')
            print(f'Shape of the cleaned ChEMBL dataframe: {ChEMBL_df.shape} \n')

        return ChEMBL_df

    #//////////////////////////////////////////////////////////
    def _count_ChEMBL_IDs(self, drugnames_list) -> int:
        '''
        Function: count ChEMBL IDs after mapping drug names to ChEMBL IDs using ChEMBL web services API.
        '''
        # Call the _check_ChEMBL_IDs function
        IDs = self._check_ChEMBL_IDs(drugnames_list)

        # Count the number of ChEMBL IDs
        missing = IDs.count(None)
        found = len(IDs) - missing
        # Missing ChEMBL IDs list
        missing_IDs = [drugnames_list[i] for i in range(len(IDs)) if IDs[i] is None]

        if self.verbose:
            print(f'Number of missing ChEMBL IDs: {missing}')
            print(f'Number of found ChEMBL IDs:{found}')
            print(f'Missing ChEMBL IDs: {missing_IDs}\n\n')

        return IDs

    #//////////////////////////////////////////////////////////
    def _check_ChEMBL_IDs(self, drugnames_list) -> list:
        '''
        Function: check ChEMBL IDs after mapping drug names to ChEMBL IDs using ChEMBL web services API.
        
        Return:
        - IDs: python list of ChEMBL IDs
        '''

        # Call the _get_ChEMBL_IDs function
        ChEMBL_ID = self._get_ChEMBL_IDs(drugnames_list)

        # Extract the ChEMBL IDs from the list
        IDs = []
        for mols in ChEMBL_ID:
            if mols:
                try:
                    value = mols[0]['molecule_chembl_id']
                except KeyError:
                    value = None   # This happens if the compound ID is not found, due to spelling of the name or different reasons
            else:
                value = None
            IDs.append(value)

        return IDs

    #//////////////////////////////////////////////////////////
    def _get_ChEMBL_IDs(self, drugnames_list):
        '''
        Function: get ChEMBL IDs from a list of drug names using ChEMBL web services API.
        Maps to molecule_synonyms_iexact.
        '''

        # Connect to the ChEMBL API
        molecule = self._connect_to_ChEMBL_API()

        # List to store the ChEMBL IDs
        ChEMBL_ID_list = []

        # Iterate over the drug names in the list
        for compound_name in drugnames_list:
            # Get the ChEMBL ID
            compound_id = molecule.filter(molecule_synonyms__molecule_synonym__iexact=compound_name).only('molecule_chembl_id')
            ChEMBL_ID_list.append(compound_id)

        return ChEMBL_ID_list

    #//////////////////////////////////////////////////////////
    def _connect_to_ChEMBL_API(self,) -> object:
        '''
        Function: connect to the ChEMBL web services API.
        '''
        # Connect to the ChEMBL API
        print('Connecting to ChEMBL API ...')
        molecule = new_client.molecule
        print('ChEMBL API connected successfully\n')

        return molecule

    #//////////////////////////////////////////////////////////
    def _load_drugnames_file(self, drugnames_file) -> list:
        '''
        Function: load a file containing drug names (csv file).
        '''
        # Try reading with header first, if it fails or has no header, try without
        try:
            drugnames_df = pandas.read_csv(drugnames_file)
            # If it has no proper columns, it might have been saved without header
            if len(drugnames_df.columns) == 0 or drugnames_df.columns[0].startswith('Unnamed'):
                drugnames_df = pandas.read_csv(drugnames_file, header=None, names=[self.column_drugname])
        except Exception:
            drugnames_df = pandas.read_csv(drugnames_file, header=None, names=[self.column_drugname])

        # Check if the header is 'drug_name', if not, rename the column to self.drugname_column
        if self.column_drugname not in drugnames_df.columns:
            drugnames_df = drugnames_df.rename(columns={drugnames_df.columns[0]: self.column_drugname})

        # Drop the duplicates in the drug names
        drugnames_df = drugnames_df.drop_duplicates(subset=[self.column_drugname])
        drugnames_list = drugnames_df[self.column_drugname].tolist()

        if self.verbose:
            print('Drug names file loaded successfully\n')
            print(f'Drug names in the list: {drugnames_list}\n')
            print(f'Number of drug names: {len(drugnames_list)}\n')

        return drugnames_list
    

    #/////////////////////////////////////////////////////////////////////////////////////////
    # BRANCH: drug_doses
    #/////////////////////////////////////////////////////////////////////////////////////////
    
    def get_drugdoses(self, doses_type: str = None) -> pandas.DataFrame:
        '''
        Function: get the drug doses in a dataframe of [drug_name, concentration, ChEMBL_ID, ChEMBL_conc] columns.
        '''
        if doses_type:
            self.doses_type = doses_type

        if self.verbose:
            print('\nGetting the drug doses dataframe ...')
            print('Drug doses type:', self.doses_type)
        
        if self.drugdoses_file:
            # Load the drug doses file
            self.drugdoses_df = self._load_drugdoses_file(self.drugdoses_file)
        elif self.drugscreen_df is not None:
            # Copy the drug screen dataframe
            self.drugdoses_df = self.drugscreen_df.copy()
        else:
            # No drugdoses_file nor drugscreen_df, fill with 10000 for all drugs
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

        if self.drugID_file:
            # Load the drug ID file
            drugID_dict = self._load_drug_ChEMBL_IDs_file()
        else:
            # Make the drug ChEMBL_IDs df into a dictionary
            if self.ChEMBL_df is None:
                self.ChEMBL_df = self._process_drugnames_to_ChEMBL()
            drugID_dict = self.ChEMBL_df.set_index(self.column_drugname)[self.column_chembl].to_dict()
        
        if self.doses_type == 'melted':
            # Melt the drug doses dataframe
            self.drugdoses_df = self._melt_drugdoses(self.drugdoses_df)

        elif self.doses_type == 'not_melted':
            # Make and clean the drugdoses dataframe
            self.drugdoses_df = self._clean_drugdoses(self.drugdoses_df)

        # Merge ChEMBL IDs with drug doses
        self.drugdoses_df = self._merge_drugID_doses(self.drugdoses_df, drugID_dict)

        if self.directory_output:
            save_file(self.drugdoses_df, self.directory_output, 'drug_ChEMBL_doses.csv', index=False)

        return self.drugdoses_df
    
    #//////////////////////////////////////////////////////////
    def _merge_drugID_doses(self, drugdoses, drugID_dict) -> pandas.DataFrame:
        '''
        Function: merge ChEMBL IDs with drug doses.
        '''
        # Map drug names to ChEMBL IDs using the drug ID dictionary
        drugdoses_IDs = self._map_drugID(drugdoses, drugID_dict)
        # Merge ChEMBL IDs with drug doses
        drugdoses_IDs['ChEMBL_conc'] = drugdoses_IDs[self.column_chembl] + '_' + drugdoses_IDs[self.column_concentration].astype(str) + '_nM'
        # Drop rows with missing ChEMBL IDs
        drugdoses_IDs = drugdoses_IDs.dropna(subset=[self.column_chembl]).reset_index(drop=True)
            
        if self.verbose:
            print('Shape of Drug doses dataframe after dropping missing ChEMBL IDs:', drugdoses_IDs.shape)
            print('Columns in drugdoses_merged:', drugdoses_IDs.columns)
            print('Head of drugdoses_merged:', drugdoses_IDs.head(3), '\n')
            print('ChEMBL IDs merged with drug doses successfully\n')

        return drugdoses_IDs
    
    #//////////////////////////////////////////////////////////
    def _map_drugID(self, drugdoses, drugID_dict):
        '''
        Function: filter or map drug names to ChEMBL IDs using the drug ID dictionary.
        '''
        drugdoses_IDs = drugdoses.copy()
        drugdoses_IDs[self.column_chembl] = drugdoses_IDs[self.column_drugname].map(drugID_dict)

        if self.verbose:
            # Find drug names not mapped to ChEMBL IDs
            print('Columns in drugdoses_IDs:', drugdoses_IDs.columns)
            print('Head of drugdoses_IDs:', drugdoses_IDs.head(3), '\n')
            print('Drug names not mapped to ChEMBL IDs:', drugdoses_IDs[drugdoses_IDs[self.column_chembl].isnull()][self.column_drugname].unique())

        return drugdoses_IDs

    #//////////////////////////////////////////////////////////
    def _melt_drugdoses(self, drugdoses) -> pandas.DataFrame:
        '''
        Function: for drugdose files with multiple doses columns, melt the dataframe to have a long format where each concentration is in a separate row.
        '''
        if self.verbose:
            print('\nMelting the drug doses dataframe ...\n')
            print('Drug doses dataframe before melting:', drugdoses.shape, '\n')
        
        # Melt the dataframe for each concentration column to be in a separate row
        # Concentration columns are all except the first column (drug names)
        drugdoses_melted = drugdoses.melt(
            id_vars=[drugdoses.columns[0]], # drug names
            value_vars=drugdoses.columns[1:], # concentration columns
            # var_name=self.column_drugname,
            var_name='dose_type',  # this can be any placeholder name
            value_name=self.column_concentration)
        
        # Keep only the necessary columns
        drugdoses_melted = drugdoses_melted[[drugdoses.columns[0], self.column_concentration]]
        drugdoses_melted = drugdoses_melted.rename(columns={drugdoses.columns[0]: self.column_drugname})

        if self.verbose:
            print('Drug doses dataframe after melting:', drugdoses_melted.shape, drugdoses_melted.columns, '\n')
            print('Head of drug doses dataframe after melting:', drugdoses_melted.head(3), '\n')
        # Clean the melted dataframe
        drugdoses_melted = self._clean_drugdoses(drugdoses_melted)

        return drugdoses_melted


    #//////////////////////////////////////////////////////////
    def _clean_drugdoses(self, drugdoses) -> pandas.DataFrame:
        '''
        Function: clean the drug doses dataframe.
        '''
        if self.verbose:
            print('Cleaning the drug doses dataframe ...')
            print('Dropping duplicates doses and doses in 0 concentration...')
            print('Columns in Drug doses dataframe before dropping duplicates:', drugdoses.columns)
            print('Shape of Drug doses dataframe before dropping duplicates:', drugdoses.shape)
            # print('Head of Drug doses dataframe before dropping duplicates:', drugdoses, '\n')

        if self.column_concentration is None:
            self.column_concentration = 'concentration'

        drugdoses = drugdoses[[self.column_drugname, self.column_concentration]].drop_duplicates().reset_index(drop=True)
        drugdoses = drugdoses[drugdoses[self.column_concentration] != 0]

        return drugdoses

    #//////////////////////////////////////////////////////////
    def _load_drug_ChEMBL_IDs_file(self,) -> dict:
        '''
        Function: load or make the drug ID file (dictionary).
        '''
        if self.drugID_file is not None:
            if self.verbose:
                print('Loading the drug ID file ...')
                print('Drug ID file:', self.drugID_file)
            drugID = pandas.read_csv(self.drugID_file)
            # Rename index to match self.
            drugID = drugID.rename(columns={drugID.columns[0]: self.column_drugname, drugID.columns[1]: self.column_chembl})
            drugID_dict = drugID.set_index(self.column_drugname)[self.column_chembl].to_dict()
        else:
            drugID = self._process_drugnames_to_ChEMBL()
            drugID_dict = drugID.set_index(self.column_drugname)[self.column_chembl].to_dict()
        
        return drugID_dict
    
    #//////////////////////////////////////////////////////////
    def _load_drugdoses_file(self, drugdoses_file) -> pandas.DataFrame:
        '''
        Function: load and read the drug doses file (csv file).
        '''
        drugdoses = pandas.read_csv(drugdoses_file)

        if self.verbose:
            print('Drug doses file:', drugdoses_file)
            print('Drug doses dataframe:', drugdoses.head(), '\n')

        return drugdoses

#/////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////

# TEST THE CLASS FUNCTION

# node_targets_file = 'C:\\Users\\viviamsb\\OneDrive - NTNU\\PhD Folder\\Pipeline\\DrugLogics_pipeline_modules\\module_scripts\\output\\jaaks_data_output_cell_fate\\supplementary_files\\drug_node_targets.csv'
# directory_output1 = 'C:\\Users\\viviamsb\\OneDrive - NTNU\\PhD Folder\\Pipeline\\DrugLogics_pipeline_modules\\module_scripts\\output\\test_class\\anchor\\'
# directory_output2 = 'C:\\Users\\viviamsb\\OneDrive - NTNU\\PhD Folder\\Pipeline\\DrugLogics_pipeline_modules\\module_scripts\\output\\test_class\\librar\\'
# db_file = 'C:\\Users\\viviamsb\\OneDrive - NTNU\\PhD Folder\\Pipeline\\DrugLogics_pipeline_modules\\module_scripts\\input\\DrugTargetInteractionDB.db'
# jaaks_drugscreen_file = 'C:\\Users\\viviamsb\\OneDrive - NTNU\\PhD Folder\\Pipeline\\DrugLogics_pipeline_modules\\module_scripts\\input\\jaaks_dataset\\jaaks_drugscreen.csv'
# jaaks_drugscreen = pandas.read_csv(jaaks_drugscreen_file)

# cell_fate_node_dict_file = 'C:\\Users\\viviamsb\\OneDrive - NTNU\\PhD Folder\\Pipeline\\DrugLogics_pipeline_modules\\module_scripts\\output\\jaaks_data_output_cell_fate\\supplementary_files\\node_dict.csv'


# jaaks_cell_fate_drugprofile_anchor = Drug(
#     drugscreen_df=jaaks_drugscreen,
#     db_file=db_file,
#     node_dict_file=cell_fate_node_dict_file,
#     directory_output=directory_output1,
#     verbose=False,
# )
# jaaks_cell_fate_drugprofile_anchor.column_drugname = 'ANCHOR_NAME'
# jaaks_cell_fate_drugprofile_anchor.column_concentration = 'ANCHOR_CONC'

# jaaks_cell_fate_drugprofile_anchor.get_drugprofiles()

# jaaks_cell_fate_drugprofile_library = Drug(
#     drugscreen_df=jaaks_drugscreen,
#     db_file=db_file,
#     node_dict_file=cell_fate_node_dict_file,
#     # directory_output=directory_output2,
#     verbose=False,
# )
# jaaks_cell_fate_drugprofile_library.column_drugname = 'LIBRARY_NAME'
# jaaks_cell_fate_drugprofile_library.column_concentration = 'LIBRARY_CONC'

# library_profiles = jaaks_cell_fate_drugprofile_library.get_drugprofiles()
# print(library_profiles.head(5))
