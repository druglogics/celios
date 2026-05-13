# Drug ChEMBL targets processing module

import pandas
from .helpers import save_file
from ...features import target_checker


class TargetProcessor:
    """
    ChEMBL targets processing class.
    
    This class handles the retrieval and processing of drug targets from ChEMBL IDs
    using the DrugTargetInteractionDB database. Optionally merges with manual targets
    from the original dataset.
    
    Args:
        drugdoses_df (pandas.DataFrame): DataFrame containing drug doses with ChEMBL IDs.
        db_file (str): Path to DrugTargetInteractionDB database file.
        column_drugname (str): Name of the column containing drug names. Default: 'drug_name'.
        column_chembl (str): Name of the column for ChEMBL IDs. Default: 'ChEMBL_ID'.
        column_concentration (str): Name of the column containing concentrations. Default: 'concentration'.
        ic50_value (float): Upper limit for IC50 value for binding affinity. Default: None.
        merge_on (str): Column name to merge on ('ChEMBL_conc' or 'ChEMBL_ID'). Default: 'ChEMBL_conc'.
        manual_targets_df (pandas.DataFrame): DataFrame with manual targets from original dataset. Default: None.
        manual_targets_column (str or list): Column name(s) containing manual targets. Default: None.
        merge_strategy (str): How to merge database and manual targets ('fill_missing', 'combine', 'override'). Default: 'fill_missing'.
        directory_output (str): Directory to save output files.
        verbose (bool): Print verbose output. Default: False.
        save (bool): Save output to file. Default: False.
    
    Methods:
        get_chembl_targets: Main method to get ChEMBL targets dataframe.
    """
    
    def __init__(self,
                 drugdoses_df: pandas.DataFrame = None,
                 db_file: str = None,
                 column_drugname: str = 'drug_name',
                 column_chembl: str = 'ChEMBL_ID',
                 column_concentration: str = 'concentration',
                 ic50_value: float = None,
                 merge_on: str = 'ChEMBL_conc',
                 manual_targets_df: pandas.DataFrame = None,
                 manual_targets_column: str = None,
                 merge_strategy: str = 'fill_missing',
                 directory_output: str = None,
                 verbose: bool = False,
                 save: bool = False):
        
        self.drugdoses_df = drugdoses_df
        self.db_file = db_file
        self.column_drugname = column_drugname
        self.column_chembl = column_chembl
        self.column_concentration = column_concentration
        self.ic50_value = ic50_value
        self.merge_on = merge_on
        self.manual_targets_df = manual_targets_df
        self.manual_targets_column = manual_targets_column
        self.merge_strategy = merge_strategy
        self.directory_output = directory_output
        self.verbose = verbose
        self.save = save
    
    def get_chembl_targets(self, 
                          drugdoses_df: pandas.DataFrame = None,
                          merge_on: str = None) -> pandas.DataFrame:
        """
        Get the ChEMBL targets dataframe.
        
        Args:
            drugdoses_df (pandas.DataFrame): Override drugdoses_df if provided.
            merge_on (str): Override merge_on parameter if provided.
            
        Returns:
            pandas.DataFrame: DataFrame with columns [drug_name, ChEMBL_ID, ChEMBL_conc, targets, no_targets].
        """
        if drugdoses_df is not None:
            self.drugdoses_df = drugdoses_df
            
        if merge_on:
            self.merge_on = merge_on
        
        if self.drugdoses_df is None:
            raise ValueError('No drugdoses_df provided. Please provide a drug doses dataframe.')
        
        if self.verbose:
            print('\nMaking the ChEMBL targets dataframe ...')
            print('Merging on:', self.merge_on, '\n')
        
        # Find targets for the ChEMBL IDs
        target_query = self._find_targets(self.drugdoses_df)
        
        # Merge the ChEMBL IDs dataframe with the target query dataframe
        chembl_targets_df = self._merge_ChEMBL_targetsquery(target_query, self.drugdoses_df)
        
        # Clean up duplicate columns and deduplicate rows before checking and merging with manual targets
        # This follows the same logic as the original notebook: drop extra columns, rename, then deduplicate
        if 'ChEMBL_ID_x' in chembl_targets_df.columns:
            chembl_targets_df = chembl_targets_df.drop(columns=['ChEMBL_ID_x'])
        if 'ChEMBL_ID_y' in chembl_targets_df.columns:
            chembl_targets_df = chembl_targets_df.rename(columns={'ChEMBL_ID_y': 'ChEMBL_ID'})
        if 'IC50' in chembl_targets_df.columns:
            chembl_targets_df = chembl_targets_df.drop(columns=['IC50'])
        
        # Drop duplicates on all columns to ensure one row per drug-concentration-target combination
        chembl_targets_df = chembl_targets_df.drop_duplicates().reset_index(drop=True)
        
        if self.verbose:
            print(f'After cleanup and deduplication: {len(chembl_targets_df)} rows\n')
        
        # Merge with manual targets if provided (BEFORE checking/counting)
        if self.manual_targets_df is not None and self.manual_targets_column is not None:
            chembl_targets_df = self._merge_manual_targets(chembl_targets_df, self.manual_targets_df)
        
        # Check the ChEMBL targets dataframe for found and missing targets (AFTER manual merge)
        chembl_targets_df = self._check_ChEMBL_targets(chembl_targets_df)
        
        # Count targets AFTER manual merge so it includes both database and manually merged targets
        if self.verbose:
            print('Counting the number of targets for each drug-concentration combination ...')
        chembl_targets_df['number_targets'] = chembl_targets_df['targets'].apply(self._count_targets).astype('Int64')
        
        if self.save and self.directory_output:
            save_file(chembl_targets_df, self.directory_output, 'drug_ChEMBL_targets.csv', index=False)
        
        if self.verbose:
            print('ChEMBL targets dataframe created successfully.\n')
        
        return chembl_targets_df
    
    def _merge_manual_targets(self, chembl_targets_df: pandas.DataFrame, manual_targets_df: pandas.DataFrame) -> pandas.DataFrame:
        """
        Merge manual targets from the original dataset with database targets.
        
        Args:
            chembl_targets_df (pandas.DataFrame): DataFrame with targets from database.
            manual_targets_df (pandas.DataFrame): DataFrame with manual targets from original dataset.
            
        Returns:
            pandas.DataFrame: Merged DataFrame with targets from database and manual sources.
        """
        if self.verbose:
            print('\nMerging manual targets with database targets ...')
            print(f'Merge strategy: {self.merge_strategy}')
            print(f'Manual target column(s): {self.manual_targets_column}\n')
        
        # Handle single or multiple manual target columns
        if isinstance(self.manual_targets_column, str):
            target_columns = [self.manual_targets_column]
        else:
            target_columns = self.manual_targets_column
        
        # OPTIMIZATION: Select only needed columns and deduplicate EARLY to avoid processing 1M+ rows
        cols_to_keep = [self.column_drugname] + target_columns
        if self.column_concentration in manual_targets_df.columns:
            cols_to_keep.insert(1, self.column_concentration)
        
        # For synergy data with drug_name_A/B, we need to handle both sides
        # Try to get columns using the standardized names first, then fall back to A/B variants
        actual_cols_to_keep = []
        for col in cols_to_keep:
            if col in manual_targets_df.columns:
                actual_cols_to_keep.append(col)
            elif col == self.column_drugname and 'drug_name_A' in manual_targets_df.columns:
                # Skip for now, we'll handle both A and B separately
                pass
            elif col == self.column_concentration and 'conc_A' in manual_targets_df.columns:
                # Skip for now, we'll handle both A and B separately
                pass
        
        # If we have drug_name_A/B with side-specific targets, expand and keep each side's targets separate.
        if 'drug_name_A' in manual_targets_df.columns and 'drug_name_B' in manual_targets_df.columns:
            if self.verbose:
                print(f'Manual targets data has dual drugs (A/B), expanding...')

            target_col_a = None
            target_col_b = None
            if isinstance(target_columns, list):
                target_col_a = next((col for col in target_columns if col.endswith('_A') and col in manual_targets_df.columns), None)
                target_col_b = next((col for col in target_columns if col.endswith('_B') and col in manual_targets_df.columns), None)
                # Fallback for custom names without _A/_B suffixes
                if target_col_a is None and len(target_columns) >= 1 and target_columns[0] in manual_targets_df.columns:
                    target_col_a = target_columns[0]
                if target_col_b is None and len(target_columns) >= 2 and target_columns[1] in manual_targets_df.columns:
                    target_col_b = target_columns[1]

            if target_col_a and target_col_b:
                # Extract Drug A using only the A-side manual targets
                manual_a = manual_targets_df[
                    ['drug_name_A'] +
                    (['conc_A'] if 'conc_A' in manual_targets_df.columns else []) +
                    [target_col_a]
                ].copy()
                manual_a.rename(columns={target_col_a: 'manual_targets'}, inplace=True)
                manual_a.rename(columns={'drug_name_A': 'drug_name'}, inplace=True)
                if 'conc_A' in manual_targets_df.columns:
                    manual_a.rename(columns={'conc_A': 'concentration'}, inplace=True)

                # Extract Drug B using only the B-side manual targets
                manual_b = manual_targets_df[
                    ['drug_name_B'] +
                    (['conc_B'] if 'conc_B' in manual_targets_df.columns else []) +
                    [target_col_b]
                ].copy()
                manual_b.rename(columns={target_col_b: 'manual_targets'}, inplace=True)
                manual_b.rename(columns={'drug_name_B': 'drug_name'}, inplace=True)
                if 'conc_B' in manual_targets_df.columns:
                    manual_b.rename(columns={'conc_B': 'concentration'}, inplace=True)

                manual_targets_df = pandas.concat([manual_a, manual_b], axis=0, ignore_index=True)
            else:
                # Extract Drug A targets
                manual_a = manual_targets_df[['drug_name_A'] + (['conc_A'] if 'conc_A' in manual_targets_df.columns else []) + target_columns].copy()
                manual_a.rename(columns={'drug_name_A': 'drug_name'}, inplace=True)
                if 'conc_A' in manual_targets_df.columns:
                    manual_a.rename(columns={'conc_A': 'concentration'}, inplace=True)
                
                # Extract Drug B targets  
                manual_b = manual_targets_df[['drug_name_B'] + (['conc_B'] if 'conc_B' in manual_targets_df.columns else []) + target_columns].copy()
                manual_b.rename(columns={'drug_name_B': 'drug_name'}, inplace=True)
                if 'conc_B' in manual_targets_df.columns:
                    manual_b.rename(columns={'conc_B': 'concentration'}, inplace=True)
                
                # Combine both sides
                manual_targets_df = pandas.concat([manual_a, manual_b], axis=0, ignore_index=True)
        else:
            # Use standardized columns if already in that format
            manual_targets_df = manual_targets_df[actual_cols_to_keep].copy()
        
        # Drop duplicates EARLY by drug name and concentration (if available)
        dedup_cols = [self.column_drugname]
        if self.column_concentration in manual_targets_df.columns:
            dedup_cols.append(self.column_concentration)
        manual_targets_df = manual_targets_df.drop_duplicates(subset=dedup_cols, keep='first')
        
        if self.verbose:
            print(f'Manual targets after deduplication: {len(manual_targets_df)} rows')
        
        # Combine all manual targets into a single column when not already normalized.
        if 'manual_targets' not in manual_targets_df.columns:
            if len(target_columns) > 1:
                # Concatenate multiple columns
                manual_targets_df['manual_targets'] = manual_targets_df[target_columns].apply(
                    lambda row: ', '.join(filter(None, row.astype(str).str.strip())), 
                    axis=1
                )
            else:
                # Use single column as-is
                manual_targets_df['manual_targets'] = manual_targets_df[target_columns[0]]
        
        # Keep only drug name, concentration (if present), and manual targets
        cols_to_keep_final = [self.column_drugname]
        if self.column_concentration in manual_targets_df.columns:
            cols_to_keep_final.append(self.column_concentration)
        cols_to_keep_final.append('manual_targets')
        
        manual_targets_df = manual_targets_df[cols_to_keep_final]
        
        # Remove rows where manual_targets is empty
        manual_targets_df = manual_targets_df[manual_targets_df['manual_targets'].notna()]
        manual_targets_df = manual_targets_df[manual_targets_df['manual_targets'].str.strip() != '']
        
        # Merge with chembl_targets_df
        chembl_targets_df = chembl_targets_df.copy()
        
        # Determine merge keys: use both drug name and concentration if available
        merge_keys = [self.column_drugname]
        if self.column_concentration in manual_targets_df.columns and self.column_concentration in chembl_targets_df.columns:
            merge_keys.append(self.column_concentration)
        
        if self.merge_strategy == 'fill_missing':
            # Use database targets, fill missing with manual targets
            merged = chembl_targets_df.merge(manual_targets_df, on=merge_keys, how='left')
            merged['targets'] = merged['targets'].fillna(merged['manual_targets'])
            chembl_targets_df = merged.drop(columns=['manual_targets'])
            
        elif self.merge_strategy == 'combine':
            # Combine both database and manual targets (union)
            merged = chembl_targets_df.merge(manual_targets_df, on=merge_keys, how='left')
            # Combine targets: keep database targets, add manual targets that aren't already there
            def combine_targets(row):
                db_targets = set(filter(None, str(row['targets']).split(', '))) if pandas.notna(row['targets']) else set()
                manual = set(filter(None, str(row['manual_targets']).split(', '))) if pandas.notna(row['manual_targets']) else set()
                combined = db_targets.union(manual)
                return ', '.join(sorted(combined)) if combined else None
            
            merged['targets'] = merged.apply(combine_targets, axis=1)
            chembl_targets_df = merged.drop(columns=['manual_targets'])
            
        elif self.merge_strategy == 'override':
            # Use manual targets if available, otherwise use database targets
            merged = chembl_targets_df.merge(manual_targets_df, on=merge_keys, how='left')
            merged['targets'] = merged['manual_targets'].fillna(merged['targets'])
            chembl_targets_df = merged.drop(columns=['manual_targets'])
        
        # Deduplicate to ensure one row per drug-concentration pair (important when merging data with multiple rows per drug)
        cols_to_check = [self.column_drugname, self.merge_on] if self.merge_on in chembl_targets_df.columns else [self.column_drugname]
        cols_to_check = [col for col in cols_to_check if col in chembl_targets_df.columns]
        if cols_to_check:
            chembl_targets_df = chembl_targets_df.drop_duplicates(subset=cols_to_check, keep='first')
        
        if self.verbose:
            print(f'Number of rows with merged targets: {len(chembl_targets_df)}')
            print('Manual targets merged successfully.\n')
        
        return chembl_targets_df
    
    def _check_ChEMBL_targets(self, chembl_targets_df) -> pandas.DataFrame:
        """
        Check the ChEMBL targets dataframe for found and missing targets.
        
        Args:
            chembl_targets_df (pandas.DataFrame): DataFrame with ChEMBL targets.
            
        Returns:
            pandas.DataFrame: Validated DataFrame.
        """
        # Fill drugs with no targets found (missing values) with NAN
        chembl_targets_df['targets'] = chembl_targets_df['targets'].replace('', None)
        
        # Count the number of drugs with found and missing targets
        missing_targets = chembl_targets_df['targets'].isnull().sum()
        found_targets = len(chembl_targets_df) - missing_targets
        
        if self.verbose:
            print(f'Number of drugs with missing targets: {missing_targets}')
            print(f'Number of drugs with found targets: {found_targets}\n')
        
        # Get the list of drugs with missing targets
        self._get_missing_ChEMBL_targets(chembl_targets_df)
        
        return chembl_targets_df
    
    def _get_missing_ChEMBL_targets(self, chembl_targets_df) -> list:
        """
        Get the list of drugs with missing targets.
        
        Args:
            chembl_targets_df (pandas.DataFrame): DataFrame with ChEMBL targets.
            
        Returns:
            list: List of drug names with missing targets.
        """
        # Get the list of drugs with missing targets
        missing_targets = chembl_targets_df[chembl_targets_df['targets'].isnull()]
        # Create a list of drugs with missing targets
        missing_targets = missing_targets[self.column_drugname].tolist()
        # Drop the duplicates in the list
        missing_targets = list(set(missing_targets))
        
        if self.verbose:
            print('List of drugs with missing targets:', missing_targets, '\n')
        
        return missing_targets
    
    def _merge_ChEMBL_targetsquery(self, target_query, drugdoses_screen) -> pandas.DataFrame:
        """
        Merge the ChEMBL IDs dataframe with the target query dataframe.
        
        Args:
            target_query (pandas.DataFrame): DataFrame containing drug IDs and HGNC symbols.
            drugdoses_screen (pandas.DataFrame): DataFrame containing drug doses.
            
        Returns:
            pandas.DataFrame: Merged DataFrame with targets.
        """
        target_query = self._split_drugIDconcentration(target_query)
        
        # Merge the ChEMBL IDs dataframe with the target query dataframe
        if self.verbose:
            print('\nMerging ChEMBL targets dataframe with target query dataframe ...\n')
        chembl_targets_df = pandas.merge(drugdoses_screen, target_query, on=self.merge_on, how='left')
        
        if self.verbose:
            print('Length of merged dataframe:', len(chembl_targets_df))
            print('Columns in the merged dataframe:', chembl_targets_df.columns, '\n')
            print('Head of the merged dataframe:', chembl_targets_df.head(), '\n')
            print('Dataframes merged successfully\n')
        
        return chembl_targets_df
    
    def _split_drugIDconcentration(self, target_query) -> pandas.DataFrame:
        """
        Rename the columns in the target query dataframe to get drugID and IC50 columns.
        
        Args:
            target_query (pandas.DataFrame): DataFrame with drug IDs and HGNC symbols.
            
        Returns:
            pandas.DataFrame: DataFrame with split columns.
        """
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
    
    def _count_targets(self, targets) -> int:
        """
        Count the number of targets for each ChEMBL ID.
        
        Args:
            targets: String of comma-separated targets or None.
            
        Returns:
            int: Number of targets.
        """
        if pandas.isna(targets):
            return None
        return len(targets.split(','))
    
    def _find_targets(self, drugdoses_screen) -> pandas.DataFrame:
        """
        Find targets for the ChEMBL IDs using the DrugTargetInteractionDB database.
        
        Args:
            drugdoses_screen (pandas.DataFrame): DataFrame containing drug names, ChEMBL IDs, and concentrations.
            
        Returns:
            pandas.DataFrame: DataFrame with drug IDs and targets.
        """
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
                    ic50_limit=self.ic50_value
                )
            else:
                target_frames = target_checker.targetProfileBA(
                    self.db_file,
                    'drugpanel_supplementary.txt',
                    id_list=[chembl_id],
                    ic50_limit=row[self.column_concentration]
                )
            
            for frame in target_frames:
                if len(frame) != 0:
                    target_holder.append(frame)
                else:
                    errorFrame = pandas.DataFrame(
                        [[chembl_id, None]],
                        columns=['drugID', 'HGNC']
                    )
                    target_holder.append(errorFrame)
        
        # Concatenate the target frames
        target_query = pandas.concat(target_holder)
        
        if self.verbose:
            print('\nTargets found successfully\n')
            print('Length of targetquery_df:', len(target_query))
            print('Number of drugs with no found targets:', target_query.isnull().sum().sum(), '\n')
        
        target_query = self._aggregate_targets(target_query)
        
        return target_query
    
    def _aggregate_targets(self, target_query) -> pandas.DataFrame:
        """
        Aggregate the targets in the target query dataframe.
        
        Args:
            target_query (pandas.DataFrame): DataFrame with targets.
            
        Returns:
            pandas.DataFrame: Aggregated DataFrame.
        """
        target_query = target_query.rename(columns={'HGNC': 'targets'})
        target_query = target_query.groupby(
            ['drugID'], as_index=False
        ).agg({
            'drugID': 'first',
            'targets': lambda x: ', '.join(filter(None, x))
        })
        return target_query
