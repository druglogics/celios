# Combinations processing module

import pandas


class CombinationProcessor:
    """
    Combinations processing class.
    
    This class prepares combinations and synergies dataframes by integrating
    ChEMBL targets, drug profiles, and drug screen data.
    
    Args:
        chembl_targets_df (pandas.DataFrame): ChEMBL targets dataframe.
            Required columns: ['drug_name', 'concentration', 'ChEMBL_conc'].
        drugprofiles_df (pandas.DataFrame): Drug profiles dataframe.
            Required columns: ['ChEMBL_conc', 'PD_profile'] (legacy 'pipeline_ID' supported).
        drugscreen_df (pandas.DataFrame): Drug screen dataframe.
            Required columns: ['cell_line', 'drug_name_A', 'drug_name_B', 'conc_A', 'conc_B']
            plus the synergy column (default 'Synergy').
        column_cell_line_name (str): Column for cell line name. Default: 'cell_line'.
        column_synergy (str): Column for synergy values. Default: 'Synergy'.
    """
    
    def __init__(self,
                 chembl_targets_df: pandas.DataFrame = None,
                 drugprofiles_df: pandas.DataFrame = None,
                 drugscreen_df: pandas.DataFrame = None,
                 column_cell_line_name: str = 'cell_line',
                 column_synergy: str = 'synergy'):
        self.chembl_targets_df = chembl_targets_df
        self.drugprofiles_df = drugprofiles_df
        self.drugscreen_df = drugscreen_df
        self.column_cell_line_name = column_cell_line_name
        self.column_synergy = column_synergy
        self.column_anchorID = 'anchor_pipeline_ID'
        self.column_libraryID = 'library_pipeline_ID'
        self.column_anchorcon = 'conc_A'
        self.column_librarycon = 'conc_B'
        self.column_anchorname = 'drug_name_A'
        self.column_libraryname = 'drug_name_B'
        self.column_anchorname_con = 'anchor_name_conc'
        self.column_libraryname_con = 'library_name_conc'
    
    def prepare_combinations(self) -> tuple[pandas.DataFrame, pandas.DataFrame]:
        if self.chembl_targets_df is None:
            raise ValueError('chembl_targets_df is required.')
        if self.drugprofiles_df is None:
            raise ValueError('drugprofiles_df is required.')
        if self.drugscreen_df is None:
            raise ValueError('drugscreen_df is required.')

        required_cols = [
            self.column_cell_line_name,
            self.column_anchorname,
            self.column_libraryname,
            self.column_anchorcon,
            self.column_librarycon,
        ]
        missing = [col for col in required_cols if col not in self.drugscreen_df.columns]
        if missing:
            raise ValueError(f"drugscreen_df is missing required columns: {missing}")

        combinations_df = self.drugscreen_df[required_cols].copy()
        combinations_df = self._add_chembl_ids(
            combinations_df,
            self.chembl_targets_df,
            self.column_anchorcon,
            self.column_librarycon,
            self.column_anchorname,
            self.column_libraryname
        )
        combinations_df = self._add_pipeline_IDs(
            combinations_df,
            self.drugprofiles_df,
            self.column_anchorID,
            self.column_libraryID
        )
        combinations_short_df = self._shorten_combinations_df(
            combinations_df,
            self.column_cell_line_name,
            self.column_anchorID,
            self.column_libraryID
        )
        synergies_obs_df = self._create_synergies_df(
            combinations_df,
            self.drugscreen_df,
            self.column_cell_line_name,
            self.column_synergy,
            self.column_anchorcon,
            self.column_librarycon,
            self.column_anchorname,
            self.column_libraryname,
            self.column_anchorname_con,
            self.column_libraryname_con
        )
        return combinations_short_df, synergies_obs_df
    
    def _add_chembl_ids(self, combinations_df, chembl_targets_df,
                        column_anchorcon, column_librarycon,
                        column_anchorname, column_libraryname):

        combinations_df = combinations_df.merge(
            chembl_targets_df[['drug_name', 'concentration', 'ChEMBL_conc']],
            left_on=[column_anchorname, column_anchorcon],
            right_on=['drug_name', 'concentration'],
            how='left'
        )
        combinations_df = combinations_df.rename(
            columns={'ChEMBL_conc': 'anchor_ChEMBL_conc'}
        ).drop(columns=['drug_name', 'concentration'])
        
        combinations_df = combinations_df.merge(
            chembl_targets_df[['drug_name', 'concentration', 'ChEMBL_conc']],
            left_on=[column_libraryname, column_librarycon],
            right_on=['drug_name', 'concentration'],
            how='left'
        )
        combinations_df = combinations_df.rename(
            columns={'ChEMBL_conc': 'library_ChEMBL_conc'}
        ).drop(columns=['drug_name', 'concentration'])
        
        return combinations_df
    
    def _add_pipeline_IDs(self, combinations_df, drugprofiles_df,
                          column_anchorID, column_libraryID):
        profile_col = 'PD_profile'
        if profile_col not in drugprofiles_df.columns:
            if 'pipeline_ID' in drugprofiles_df.columns:
                profile_col = 'pipeline_ID'
            else:
                raise ValueError(
                    "drugprofiles_df is missing 'PD_profile' (or legacy 'pipeline_ID')."
                )

        combinations_df = combinations_df.merge(
            drugprofiles_df[['ChEMBL_conc', profile_col]],
            left_on='anchor_ChEMBL_conc',
            right_on='ChEMBL_conc',
            how='left'
        )
        combinations_df = combinations_df.rename(
            columns={profile_col: column_anchorID}
        ).drop(columns=['ChEMBL_conc'])

        combinations_df = combinations_df.merge(
            drugprofiles_df[['ChEMBL_conc', profile_col]],
            left_on='library_ChEMBL_conc',
            right_on='ChEMBL_conc',
            how='left'
        )
        combinations_df = combinations_df.rename(
            columns={profile_col: column_libraryID}
        ).drop(columns=['ChEMBL_conc'])

        return combinations_df
    
    def _shorten_combinations_df(self, combinations_df,
                                 column_cell_line_name,
                                 column_anchorID,
                                 column_libraryID):
        combinations_short_df = combinations_df[[
            column_cell_line_name,
            column_anchorID,
            column_libraryID
        ]]
        combinations_short_df = combinations_short_df.dropna(
            subset=[column_anchorID, column_libraryID]
        )
        combinations_short_df = combinations_short_df.drop_duplicates(
            keep='first'
        ).reset_index(drop=True)
        return combinations_short_df
    
    def _create_synergies_df(self, combinations_df, drugscreen_df,
                             column_cell_line_name,
                             column_synergy,
                             column_anchorcon,
                             column_librarycon,
                             column_anchorname,
                             column_libraryname,
                             column_anchorname_con,
                             column_libraryname_con):
        synergies_obs_df = combinations_df.copy()
        synergies_obs_df = synergies_obs_df.drop_duplicates(keep='first').reset_index(drop=True)

        synergies_obs_df = synergies_obs_df.merge(
            drugscreen_df[[
                column_cell_line_name,
                column_anchorname,
                column_libraryname,
                column_anchorcon,
                column_librarycon,
                column_synergy
            ]],
            left_on=[
                column_cell_line_name,
                column_anchorname,
                column_libraryname,
                column_anchorcon,
                column_librarycon
            ],
            right_on=[
                column_cell_line_name,
                column_anchorname,
                column_libraryname,
                column_anchorcon,
                column_librarycon
            ],
            how='left'
        )

        synergies_obs_df = synergies_obs_df.dropna(subset=['anchor_pipeline_ID', 'library_pipeline_ID'])

        synergies_obs_df[column_anchorname_con] = (
            synergies_obs_df[column_anchorname] + '_' +
            synergies_obs_df[column_anchorcon].astype(str)
        )
        synergies_obs_df[column_libraryname_con] = (
            synergies_obs_df[column_libraryname] + '_' +
            synergies_obs_df[column_librarycon].astype(str)
        )

        return synergies_obs_df
