# Drug panel creation module

import pandas
import os


class PanelMaker:
    """
    Drug panel creation class.
    
    This class creates a drug panel file and panel dataframe from drug profiles.
    
    Args:
        directory_main_output (str): Output directory for drugpanel file.
        directory_suppl_output (str): Output directory for supplementary files.
        activation_drugs (list): List of PD_profiles corresponding to activation drugs.
        drugprofiles_df (pandas.DataFrame): Drug profiles dataframe.
        drugprofiles_file (str): Path to drug profiles file.
        save (bool): Save output files. Default: False.
    """
    
    def __init__(self,
                 directory_main_output: str = None,
                 directory_suppl_output: str = None,
                 activation_drugs=None,
                 drugprofiles_df: pandas.DataFrame = None,
                 drugprofiles_file: str = None,
                 save: bool = False):
        
        self.directory_main_output = directory_main_output
        self.directory_suppl_output = directory_suppl_output
        self.activation_drugs = activation_drugs
        self.drugprofiles_df = drugprofiles_df
        self.drugprofiles_file = drugprofiles_file
        self.save = save
    
    def get_drugpanel(self) -> pandas.DataFrame:
        """
        Create the drug panel dataframe and optionally write files.
        
        Returns:
            pandas.DataFrame: Drug panel dataframe.
        """
        if self.drugprofiles_df is None:
            if self.drugprofiles_file:
                self.drugprofiles_df = pandas.read_csv(self.drugprofiles_file)
            else:
                raise ValueError('drugprofiles_df or drugprofiles_file is required.')
        
        panel_df = self._write_drugpanel(self.drugprofiles_df)
        return panel_df
    
    def _write_drugpanel(self, drugprofiles_df) -> pandas.DataFrame:
        """
        Write the drug panel and return the panel dataframe.
        """
        panel_df = self._make_drugpanel_df(drugprofiles_df)
        
        if not self.save:
            return panel_df
        
        if self.directory_main_output is None:
            raise ValueError('directory_main_output is required when save=True.')
        
        # Identify activation drugs if provided
        activation_PD_profiles = []
        if self.activation_drugs is not None:
            activation_PD_profiles = panel_df[
                panel_df['PD_profile'].isin(self.activation_drugs)
            ]['PD_profile'].tolist()
        
        # Write the drug panel in a file
        os.makedirs(self.directory_main_output, exist_ok=True)
        with open(os.path.join(self.directory_main_output, 'drugpanel'), 'w') as f:
            f.write("#name\taction\ttarget\ttarget_n\n")
            try:
                for _, row in panel_df.iterrows():
                    PD_profile_str = row["PD_profile"]
                    targets = "\t".join(row["node_targets"])
                    action_type = "inhibits"
                    if self.activation_drugs is not None and PD_profile_str in activation_PD_profiles:
                        action_type = "activates"
                    f.write(f"{PD_profile_str}\t{action_type}\t{targets}\n")
            except Exception as e:
                print(f"Error occurred while processing the DataFrame: {e}")
        
        return panel_df
    
    def _make_drugpanel_df(self, drugprofiles_df) -> pandas.DataFrame:
        """
        Make the drug panel dataframe and optionally save it.
        """
        # Keep the columns of interest
        panel_df = drugprofiles_df[['PD_profile', 'node_targets']].copy()
        
        # Convert list columns to strings to make them hashable
        panel_df['node_targets'] = panel_df['node_targets'].apply(
            lambda x: str(x) if isinstance(x, list) else x
        )
        
        # Drop duplicates
        panel_df = panel_df.drop_duplicates().reset_index(drop=True)
        
        # Modify the value format from the 'node_targets' column
        panel_df['node_targets'] = panel_df['node_targets'].apply(
            lambda x: eval(x) if isinstance(x, str) and x.startswith('[') and x.endswith(']') else x
        )
        
        if self.save:
            self._save_drugpanel_df(panel_df)
        
        return panel_df
    
    def _save_drugpanel_df(self, panel_df):
        """
        Save the drug panel dataframe in a file.
        """
        if self.directory_suppl_output is None:
            raise ValueError('directory_suppl_output is required when save=True.')
        os.makedirs(self.directory_suppl_output, exist_ok=True)
        panel_df.to_csv(os.path.join(self.directory_suppl_output, 'drug_panel_df.csv'), index=False)
