"""
File upload validation and preview utilities
"""

import pandas as pd
import streamlit as st
from io import StringIO


def validate_and_preview_files(files_dict):
    """
    Validate uploaded files and return preview + validation messages.

    Returns:
        tuple: (validation_ok: bool, messages: list)
    """
    messages = []
    validation_ok = True

    # Check required files
    required = ["drug_names_file", "synergy_data_file", "node_dict_file"]
    for file_key in required:
        if not files_dict.get(file_key):
            messages.append(f"✗ Missing required: {file_key.replace('_file', '')}")
            validation_ok = False

    # Validate drug_names.txt
    if files_dict["drug_names_file"]:
        try:
            content = files_dict["drug_names_file"].read().decode("utf-8")
            lines = [line.strip() for line in content.split("\n") if line.strip()]
            messages.append(f"✓ drug_names.txt: {len(lines)} drugs")
            files_dict["drug_names_file"].seek(0)  # Reset for later read
        except Exception as e:
            messages.append(f"✗ Error reading drug_names.txt: {str(e)}")
            validation_ok = False

    # Validate synergy_data.csv
    if files_dict["synergy_data_file"]:
        try:
            df = pd.read_csv(files_dict["synergy_data_file"])
            messages.append(f"✓ synergy_data.csv: {len(df)} rows, {len(df.columns)} columns")
            files_dict["synergy_data_file"].seek(0)

            # Preview
            with st.expander("📋 synergy_data.csv Preview"):
                st.dataframe(df.head(5))
        except Exception as e:
            messages.append(f"✗ Error reading synergy_data.csv: {str(e)}")
            validation_ok = False

    # Validate node_dict.csv
    if files_dict["node_dict_file"]:
        try:
            df = pd.read_csv(files_dict["node_dict_file"])
            messages.append(f"✓ node_dict.csv: {len(df)} rows, {len(df.columns)} columns")
            files_dict["node_dict_file"].seek(0)

            with st.expander("📋 node_dict.csv Preview"):
                st.dataframe(df.head(5))
        except Exception as e:
            messages.append(f"✗ Error reading node_dict.csv: {str(e)}")
            validation_ok = False

    # Optional files
    if files_dict["tissue_cline_file"]:
        try:
            df = pd.read_csv(files_dict["tissue_cline_file"])
            messages.append(f"✓ tissue_cline.csv: {len(df)} rows")
            files_dict["tissue_cline_file"].seek(0)
        except Exception as e:
            messages.append(f"⚠️ Warning reading tissue_cline.csv: {str(e)}")

    if files_dict["manual_chembl_file"]:
        try:
            df = pd.read_csv(files_dict["manual_chembl_file"])
            messages.append(f"✓ manual_chembl.csv: {len(df)} rows")
            files_dict["manual_chembl_file"].seek(0)
        except Exception as e:
            messages.append(f"⚠️ Warning reading manual_chembl.csv: {str(e)}")

    return validation_ok, messages


def get_file_summary(files_dict):
    """Get summary of uploaded files"""
    summary = {}
    for key, file in files_dict.items():
        if file:
            summary[key] = f"{file.name} ({file.size} bytes)"
        else:
            summary[key] = "Not provided"
    return summary
