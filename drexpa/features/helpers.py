# module with helper functions

import os
import pandas
import json

#//////////////////////////////////////////////////////////

def save_file(data, directory_path, file_name, file_type='csv', **kwargs):
    """
    Save data to a file in the specified format.

    Parameters:
    - data: The data to be saved. Should be a DataFrame for 'csv', or a dict for 'json'.
    - directory_path: Directory where the file will be saved.
    - file_name: Name of the file to save.
    - file_type: Type of the file to save. Can be 'csv' or 'json'.
    - kwargs: Additional arguments passed to the saving function. E.g., index=False for CSV.
    """
    # Ensure directory exists
    if directory_path and not os.path.exists(directory_path):
        os.makedirs(directory_path, exist_ok=True)
    
    # Build full file path
    file_path = os.path.join(directory_path, file_name) if directory_path else file_name
    
    if file_type == 'csv':
        if isinstance(data, pandas.DataFrame):
            data.to_csv(file_path, **kwargs)
            print(f'Data saved to {file_name} as CSV, in {directory_path}.')
        else:
            raise ValueError('For CSV, data should be a pandas DataFrame.')
    elif file_type == 'json':
        if isinstance(data, dict):
            with open(file_path, 'w') as file:
                json.dump(data, file)
            print(f'Data saved to {file_name} as JSON, in {directory_path}.')
        else:
            raise ValueError('For JSON, data should be a dictionary.')
    elif file_type == 'txt':
        if isinstance(data, str):
            with open(file_path, 'w') as file:
                file.write(data)
            print(f'Data saved to {file_name} as TXT, in {directory_path}.')
        else:
            raise ValueError('For TXT, data should be a string.')
    else:
        raise ValueError('Unsupported file type. Use "csv", "json", or "txt.')
    

#//////////////////////////////////////////////////////////
def load_csv_file(file_path, verbose=False, **kwargs):
    """
    Load a CSV file and return the data as a DataFrame.

    Parameters:
    - file_path: Path to the CSV file.

    Returns:
    - data: DataFrame with the data from the CSV file.
    """
    data = pandas.read_csv(file_path, **kwargs)

    if verbose:
        print(f'Data loaded from {file_path}.')
        print(f'Columns in the data frame: {data.columns}')
    return data