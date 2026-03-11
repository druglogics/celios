import os
import logging

logger = logging.getLogger(__name__)

#/////////////////////////////////////////////////////////////////////////////////////////

def DL_trainingfiles(
        activity_matrix,
        directory_cell_lines,
        proliferation_state = True,
        dna_damage_state = False,
        report=False
        ):
    '''
    Function: write the training files for each cell line
    Input: activity_matrix - dataframe with the activity values for each node and cell line
                directory_cell_lines - directory to save the training files
                proliferation_state - boolean to indicate if the proliferation state is included
    '''
    # Allow passing a path to a CSV file instead of a DataFrame
    if isinstance(activity_matrix, (str, os.PathLike)):
        import pandas as pd
        if not os.path.exists(activity_matrix):
            raise FileNotFoundError(f"activity_matrix path not found: {activity_matrix}")
        activity_matrix = pd.read_csv(activity_matrix)
    elif not hasattr(activity_matrix, 'columns'):
        raise TypeError("activity_matrix must be a pandas DataFrame or a path to a CSV file")

    # Call reporting only when requested
    if report:
        try:
            from celios.utils import report as report_mod
            report_mod.cellfiles_report(activity_matrix, directory_cell_lines, proliferation_state=proliferation_state, dna_damage_state=dna_damage_state)
        except Exception as e:
            logger.debug("Reporting skipped or failed: %s", e)

    # Create the training text file for each cell line
    logger.info('Extracting the cell line names and creating the training files in directory...')
    logger.info('Cell line directory: %s', directory_cell_lines)
    for column in activity_matrix.columns[1:]:
        try:
            cell_line_name = column  # Extract column name as cell_line_name
            # Remove non-alphanumeric characters (e.g. spaces, slashes) from cell_line_name and make all uppercase
            cell_line_name = ''.join(e for e in cell_line_name if e.isalnum()).upper()
            # Create the directory for the specific cell_line if it doesn't exist
            cell_line_dir = os.path.join(directory_cell_lines, cell_line_name)
            os.makedirs(cell_line_dir, exist_ok=True)
            logger.info('Cell line directory created for: %s', cell_line_name)

            file_name = os.path.join(cell_line_dir, f"{cell_line_name}_training")

            with open(file_name, 'w') as file:
                file.write(f"# {column}\n")
                file.write("Condition\n-\nResponse\n")
                for index, row in activity_matrix.iterrows():
                    # The activity matrix uses the node name as the DataFrame index
                    node_name = index
                    # get value for this cell line/column
                    try:
                        value = row[column]
                    except Exception:
                        # Column missing for this row; skip
                        continue
                    # skip missing values (NaN)
                    try:
                        import pandas as _pd
                        if _pd.isna(value):
                            continue
                    except Exception:
                        if value is None or value != value:
                            continue
                    file.write(f"{node_name}:{value}\t")
                file.write("\nWeight:1\n")

                # Add the proliferation state if it is included
                if proliferation_state:
                    file.write('\n# Proliferation state\n')
                    file.write("Condition\n-\nResponse\n")
                    file.write("globaloutput:1")
                    file.write("\nWeight:1\n")

                # Add a dna damage state
                if dna_damage_state:
                    file.write('\n# DNA damage state\n')
                    file.write("Condition\n-\nResponse\n")
                    file.write("DSB_event:1\tSSB_event:1")
                    file.write("\nWeight:1\n")

                logger.info("Training file created for cell line: %s", cell_line_name)
        except Exception as e:
            logger.exception("Error processing column %s: %s", column, e)
    logger.info('-------------------------------------------------------------------')
    logger.info('\nTraining files created for each cell line successfully\n')


#/////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////