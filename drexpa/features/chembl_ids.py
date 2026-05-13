# ChEMBL IDs resolution module

import json
import os
import re
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

import pandas
from chembl_webresource_client.new_client import new_client

from .helpers import save_file


class ChEMBLIDResolver:
    """
    ChEMBL ID resolution class.
    
    This class handles the resolution of drug names or PubChemIDs to ChEMBL IDs 
    using the ChEMBL web services API.
    
    Args:
        drugnames_file (str): Path to file containing drug names (CSV). Optional if using dataframe or pubchem_ids.
        drugnames_df (pandas.DataFrame): DataFrame containing drug names. Optional if using file or pubchem_ids.
        column_drugname (str): Name of the column containing drug names. Default: 'drug_name'.
        column_chembl (str): Name of the column for ChEMBL IDs. Default: 'ChEMBL_ID'.
        directory_output (str): Directory to save output files.
        verbose (bool): Print verbose output. Default: False.
        save (bool): Save output to file. Default: False.
        manual_chembl_file (str): Path to manual ChEMBL mapping CSV for fallback resolution.
        pubchem_ids (list): List of PubChemIDs for direct resolution. Optional.
        column_pubchem (str): Name of the column for PubChemID. Default: 'PubChemID'.
    
    Workflow Options:
        1. Drug names only: Pass drugnames_file or drugnames_df with drug names
        2. PubChemID only: Pass pubchem_ids list for direct PubChemID→ChEMBL resolution
        3. Hybrid (recommended): Pass drugnames_df with both drug_name and PubChemID columns
           - First tries drug name lookup
           - Falls back to PubChemID if drug name not found
    
    Methods:
        resolve_chembl_ids: Main method to resolve drug names/PubChemIDs to ChEMBL IDs.
    """
    
    def __init__(self,
                 drugnames_file: str = None,
                 drugnames_df: pandas.DataFrame = None,
                 column_drugname: str = 'drug_name',
                 column_chembl: str = 'ChEMBL_ID',
                 directory_output: str = None,
                 verbose: bool = False,
                 save: bool = False,
                 manual_chembl_file: str = None,
                 pubchem_ids: list = None,
                 column_pubchem: str = 'PubChemID'):
        
        self.drugnames_file = drugnames_file
        self.drugnames_df = drugnames_df
        self.column_drugname = column_drugname
        self.column_chembl = column_chembl
        self.directory_output = directory_output
        self.verbose = verbose
        self.save = save
        self.manual_chembl_file = manual_chembl_file
        self.pubchem_ids = pubchem_ids
        self.column_pubchem = column_pubchem
    
    def resolve_chembl_ids(self) -> pandas.DataFrame:
        """
        Main method to resolve drug names or PubChemIDs to ChEMBL IDs.
        
        Returns:
            pandas.DataFrame: DataFrame with drug names/PubChemIDs and ChEMBL IDs.
        """
        if self.verbose:
            print('\nProcessing to ChEMBL IDs ...')
        
        # Handle PubChemID resolution
        if self.pubchem_ids is not None and len(self.pubchem_ids) > 0:
            if self.verbose:
                print('Using PubChemID list for resolution ...')
            ChEMBL_df = self._make_df_ChEMBL_IDs_from_pubchem(self.pubchem_ids)
        # Handle drug names from dataframe
        elif self.drugnames_df is not None:
            if self.verbose:
                print('Using drug names dataframe ...')
            drugnames_df = self.drugnames_df.drop_duplicates(subset=[self.column_drugname])
            drugnames_list = drugnames_df[self.column_drugname].tolist()
            pubchem_list = drugnames_df[self.column_pubchem].tolist() if self.column_pubchem in drugnames_df.columns else None
            ChEMBL_df = self._make_df_ChEMBL_IDs(drugnames_list, pubchem_list)
        # Handle drug names from file
        elif self.drugnames_file:
            if self.verbose:
                print('Using drug names file:', self.drugnames_file, '\n')
            drugnames_list = self._load_drugnames_file(self.drugnames_file)
            ChEMBL_df = self._make_df_ChEMBL_IDs(drugnames_list)
        else:
            raise ValueError('No drug names file, drug names dataframe, or PubChemID list provided.')
        
        if self.save and self.directory_output:
            save_file(ChEMBL_df, self.directory_output, 'drug_ChEMBL_IDs.csv', index=False)
        
        return ChEMBL_df
    
    def _make_df_ChEMBL_IDs(self, drugnames_list, pubchem_list=None) -> pandas.DataFrame:
        """
        Make a dataframe with drug names and ChEMBL IDs.
        
        Args:
            drugnames_list (list): List of drug names.
            pubchem_list (list, optional): List of PubChemIDs. If provided, used as fallback when drug name lookup fails.
            
        Returns:
            pandas.DataFrame: DataFrame with drug names and ChEMBL IDs.
        """
        # Call the _count_ChEMBL_IDs function
        IDs = self._count_ChEMBL_IDs(drugnames_list, pubchem_list)
        
        ChEMBL_df = pandas.DataFrame({
            self.column_drugname: drugnames_list,
            self.column_chembl: IDs
        })
        
        if self.verbose:
            print('Processed drug names to ChEMBL symbols.')
            print('Columns in the ChEMBL IDs dataframe:', ChEMBL_df.columns)
        
        # Fill missing IDs from manual CSV if available before cleaning
        manual_dict = self._load_manual_chembl_dict()
        if manual_dict:
            for idx, row in ChEMBL_df.iterrows():
                if pandas.isna(row[self.column_chembl]) and row[self.column_drugname] in manual_dict:
                    ChEMBL_df.at[idx, self.column_chembl] = manual_dict[row[self.column_drugname]]
                    if self.verbose:
                        print(f"Filled ChEMBL ID for '{row[self.column_drugname]}' from manual CSV: {manual_dict[row[self.column_drugname]]}")
        
        # Clean NaN from ChEMBL dataframe
        ChEMBL_df = self._clean_nan_ChEMBL_IDs(ChEMBL_df)
        
        return ChEMBL_df
    
    def _clean_nan_ChEMBL_IDs(self, ChEMBL_df):
        """
        Clean NaN from ChEMBL dataframe.
        
        Args:
            ChEMBL_df (pandas.DataFrame): DataFrame with ChEMBL IDs.
            
        Returns:
            pandas.DataFrame: Cleaned DataFrame.
        """
        ChEMBL_df = ChEMBL_df.dropna().reset_index(drop=True)
        
        if self.verbose:
            print('\nCleaning NaN from ChEMBL dataframe ...')
            print(f'Shape of the cleaned ChEMBL dataframe: {ChEMBL_df.shape} \n')
        
        return ChEMBL_df
    
    def _count_ChEMBL_IDs(self, drugnames_list, pubchem_list=None) -> list:
        """
        Count ChEMBL IDs after mapping drug names, with PubChemID fallback.
        
        Args:
            drugnames_list (list): List of drug names.
            pubchem_list (list, optional): List of PubChemIDs for fallback.
            
        Returns:
            list: List of ChEMBL IDs.
        """
        # Call the _check_ChEMBL_IDs function
        IDs = self._check_ChEMBL_IDs(drugnames_list, pubchem_list)
        
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
    
    def _check_ChEMBL_IDs(self, drugnames_list, pubchem_list=None) -> list:
        """
        Check ChEMBL IDs after mapping drug names, with PubChemID fallback.
        
        Args:
            drugnames_list (list): List of drug names.
            pubchem_list (list, optional): List of PubChemIDs for fallback.
            
        Returns:
            list: List of ChEMBL IDs.
        """
        # Call the _get_ChEMBL_IDs function
        ChEMBL_ID = self._get_ChEMBL_IDs(drugnames_list, pubchem_list)
        
        # Extract the ChEMBL IDs from the list
        IDs = []
        for mols in ChEMBL_ID:
            if mols:
                try:
                    value = mols[0]['molecule_chembl_id']
                except KeyError:
                    value = None  # This happens if the compound ID is not found
            else:
                value = None
            IDs.append(value)
        
        return IDs
    
    def _get_ChEMBL_IDs(self, drugnames_list, pubchem_list=None):
        """
        Get ChEMBL IDs from drug names using ChEMBL web services API.
        If drug name lookup fails, tries PubChemID lookup as fallback.
        
        Args:
            drugnames_list (list): List of drug names.
            pubchem_list (list, optional): List of PubChemIDs for fallback.
            
        Returns:
            list: List of ChEMBL ID query results.
        """
        # Connect to the ChEMBL API
        molecule = self._connect_to_ChEMBL_API()
        
        # List to store the ChEMBL IDs
        ChEMBL_ID_list = []
        
        # Iterate over the drug names in the list
        for idx, compound_name in enumerate(drugnames_list):
            compound_id = self._resolve_identifier_to_chembl(molecule, compound_name)
            
            # If primary lookup failed and a PubChemID fallback list is available, try it
            if (not compound_id) and pubchem_list and idx < len(pubchem_list):
                pubchem_id = pubchem_list[idx]
                if pandas.notna(pubchem_id):
                    if self.verbose:
                        print(f"Identifier '{compound_name}' not found in ChEMBL, trying PubChemID fallback: {pubchem_id}")
                    compound_id = self._resolve_pubchem_identifier(molecule, pubchem_id)
            
            ChEMBL_ID_list.append(compound_id)
        
        return ChEMBL_ID_list

    def _resolve_identifier_to_chembl(self, molecule, identifier):
        """Resolve a single identifier to a ChEMBL query result."""
        if pandas.isna(identifier):
            return []

        identifier_text = str(identifier).strip()
        if not identifier_text:
            return []

        if self._looks_like_chembl_id(identifier_text):
            if self.verbose:
                print(f"Identifier '{identifier_text}' already looks like a ChEMBL ID")
            return [{'molecule_chembl_id': identifier_text}]

        if self._looks_like_pubchem_cid(identifier_text):
            if self.verbose:
                print(f"Identifier '{identifier_text}' looks like a PubChem CID")
            resolved = self._resolve_pubchem_identifier(molecule, identifier_text)
            if resolved:
                return resolved

        # Default: treat as a compound name and query ChEMBL first.
        compound_id = molecule.filter(
            molecule_synonyms__molecule_synonym__iexact=identifier_text
        ).only('molecule_chembl_id')

        if compound_id:
            return compound_id

        # If the ChEMBL synonym lookup failed, ask PubChem for candidate CIDs by name.
        if self.verbose:
            print(f"Identifier '{identifier_text}' not found in ChEMBL synonym lookup, querying PubChem ...")

        return self._resolve_pubchem_identifier(molecule, identifier_text)

    def _resolve_pubchem_identifier(self, molecule, identifier):
        """Resolve a PubChem CID or name through PubChem PUG REST, then map to ChEMBL."""
        pubchem_cids = []
        identifier_text = str(identifier).strip()

        if self._looks_like_pubchem_cid(identifier_text):
            pubchem_cids = [str(int(float(identifier_text)))]
        else:
            pubchem_cids = self._query_pubchem_cids_by_name(identifier_text)

        for pubchem_cid in pubchem_cids:
            inchikey = self._query_pubchem_inchikey(pubchem_cid)
            if not inchikey:
                continue

            try:
                compound_id = molecule.filter(
                    molecule_structures__standard_inchi_key__iexact=inchikey
                ).only('molecule_chembl_id')
                if compound_id:
                    if self.verbose:
                        print(f"PubChem CID {pubchem_cid} (InChIKey {inchikey}) → {compound_id[0]['molecule_chembl_id']}")
                    return compound_id
            except Exception as e:
                if self.verbose:
                    print(f"PubChem→ChEMBL mapping failed for CID {pubchem_cid}: {e}")

        return []

    def _query_pubchem_cids_by_name(self, compound_name):
        """Query PubChem for candidate CIDs by compound name."""
        encoded_name = quote(compound_name)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded_name}/cids/JSON"

        try:
            payload = self._fetch_json(url)
        except Exception as e:
            if self.verbose:
                print(f"PubChem name lookup failed for '{compound_name}': {e}")
            return []

        try:
            identifiers = payload['IdentifierList']['CID']
        except (KeyError, TypeError):
            return []

        return [str(cid) for cid in identifiers[:3]]

    def _query_pubchem_inchikey(self, pubchem_cid):
        """Query PubChem PUG REST for the InChIKey of a CID."""
        url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
            f"{quote(str(pubchem_cid))}/property/InChIKey/JSON"
        )

        try:
            payload = self._fetch_json(url)
        except Exception as e:
            if self.verbose:
                print(f"PubChem InChIKey lookup failed for CID {pubchem_cid}: {e}")
            return None

        try:
            properties = payload['PropertyTable']['Properties']
            if properties:
                return properties[0].get('InChIKey')
        except (KeyError, TypeError, IndexError):
            return None

        return None

    def _fetch_json(self, url, timeout=15):
        """Fetch JSON from a URL using the standard library."""
        request = Request(url, headers={'Accept': 'application/json'})
        try:
            with urlopen(request, timeout=timeout) as response:
                return json.loads(response.read().decode('utf-8'))
        except (HTTPError, URLError, TimeoutError, ValueError) as e:
            raise RuntimeError(f"Unable to fetch JSON from {url}: {e}") from e

    def _looks_like_chembl_id(self, identifier_text):
        """Return True when the identifier already looks like a ChEMBL ID."""
        return bool(re.fullmatch(r'CHEMBL\d+', identifier_text.strip().upper()))

    def _looks_like_pubchem_cid(self, identifier_text):
        """Return True when the identifier is a numeric PubChem CID."""
        return identifier_text.strip().isdigit()
    
    def _connect_to_ChEMBL_API(self) -> object:
        """
        Connect to the ChEMBL web services API.
        
        Returns:
            object: ChEMBL molecule client.
        """
        print('Connecting to ChEMBL API ...')
        molecule = new_client.molecule
        print('ChEMBL API connected successfully\n')
        
        return molecule
    
    def _load_drugnames_file(self, drugnames_file) -> list:
        """
        Load a file containing drug names (CSV file).
        
        Args:
            drugnames_file (str): Path to drug names file.
            
        Returns:
            list: List of drug names.
        """
        # Try reading with header first, if it fails or has no header, try without
        try:
            drugnames_df = pandas.read_csv(drugnames_file)
            # If it has no proper columns, it might have been saved without header
            if len(drugnames_df.columns) == 0 or drugnames_df.columns[0].startswith('Unnamed'):
                drugnames_df = pandas.read_csv(drugnames_file, header=None, names=[self.column_drugname])
        except Exception:
            drugnames_df = pandas.read_csv(drugnames_file, header=None, names=[self.column_drugname])
        
        # Check if the header is 'drug_name', if not, rename the column
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
    
    def _make_df_ChEMBL_IDs_from_pubchem(self, pubchem_list) -> pandas.DataFrame:
        """
        Make a dataframe with PubChemIDs and resolved ChEMBL IDs.
        Used when only PubChemID list is provided.
        
        Args:
            pubchem_list (list): List of PubChemIDs.
            
        Returns:
            pandas.DataFrame: DataFrame with PubChemIDs and ChEMBL IDs.
        """
        if self.verbose:
            print('Resolving ChEMBL IDs from PubChemIDs ...')
        
        IDs = self._get_ChEMBL_IDs_from_pubchem(pubchem_list)
        
        ChEMBL_df = pandas.DataFrame({
            self.column_pubchem: pubchem_list,
            self.column_chembl: IDs
        })
        
        if self.verbose:
            print('Processed PubChemIDs to ChEMBL IDs.')
            print('Columns in the ChEMBL IDs dataframe:', ChEMBL_df.columns)
        
        # Clean NaN from ChEMBL dataframe
        ChEMBL_df = self._clean_nan_ChEMBL_IDs(ChEMBL_df)
        
        return ChEMBL_df
    
    def _get_ChEMBL_IDs_from_pubchem(self, pubchem_list) -> list:
        """
        Get ChEMBL IDs directly from PubChemIDs using ChEMBL web services API.
        
        Args:
            pubchem_list (list): List of PubChemIDs.
            
        Returns:
            list: List of ChEMBL IDs.
        """
        molecule = self._connect_to_ChEMBL_API()
        
        IDs = []
        for pubchem_id in pubchem_list:
            if pandas.notna(pubchem_id):
                try:
                    # Query ChEMBL by PubChemID (stored as external_id)
                    compound_id = molecule.filter(
                        external_ids__external_id=str(int(pubchem_id))
                    ).only('molecule_chembl_id')
                    
                    if compound_id:
                        chembl_id = compound_id[0]['molecule_chembl_id']
                        if self.verbose:
                            print(f"PubChemID {pubchem_id} → {chembl_id}")
                        IDs.append(chembl_id)
                    else:
                        if self.verbose:
                            print(f"PubChemID {pubchem_id} → Not found")
                        IDs.append(None)
                except Exception as e:
                    if self.verbose:
                        print(f"Error querying PubChemID {pubchem_id}: {e}")
                    IDs.append(None)
            else:
                IDs.append(None)
        
        return IDs
    
    def _load_manual_chembl_dict(self) -> dict:
        """
        Load manual ChEMBL IDs from CSV file into a dictionary.
        
        Returns:
            dict: Dictionary mapping drug names to ChEMBL IDs.
        """
        if not self.manual_chembl_file or not os.path.exists(self.manual_chembl_file):
            return {}
        
        try:
            manual_df = pandas.read_csv(self.manual_chembl_file)
            if self.column_drugname not in manual_df.columns or self.column_chembl not in manual_df.columns:
                if self.verbose:
                    print(f"Warning: Manual ChEMBL CSV missing required columns '{self.column_drugname}' or '{self.column_chembl}'. Skipping.")
                return {}
            manual_dict = dict(zip(manual_df[self.column_drugname], manual_df[self.column_chembl]))
            if self.verbose:
                print(f"Loaded {len(manual_dict)} manual ChEMBL IDs from {self.manual_chembl_file}")
            return manual_dict
        except Exception as e:
            if self.verbose:
                print(f"Error loading manual ChEMBL CSV: {e}. Skipping.")
            return {}
