"""
Pipeline execution and configuration utilities
"""

import os
import sys
import time
import pandas as pd
import streamlit as st
from pathlib import Path
import traceback
from contextlib import redirect_stdout, redirect_stderr
from io import StringIO


def build_config(app_config, files_info, work_dir):
    """
    Build a DREXPA config dict from app configuration and file paths.

    Args:
        app_config: Configuration from Streamlit sidebar
        files_info: Dictionary with file paths
        work_dir: Working directory for output

    Returns:
        dict: DREXPA configuration
    """
    config = {
        "global": {
            "output_dir": os.path.join(work_dir, app_config["global"]["output_dir"]),
            "base_data_dir": work_dir,
            "verbose": app_config["global"]["verbose"],
            "save": app_config["global"]["save"],
        },
        "paths": {
            "drug_names_file": files_info.get("drug_names_file", "drug_names.txt"),
            "synergy_data_file": files_info.get("synergy_data_file", "synergy_data.csv"),
            "node_dict_file": files_info.get("node_dict_file", "node_dict.csv"),
            "tissue_cline_file": files_info.get("tissue_cline_file", "tissue_cline.csv"),
            "manual_chembl_csv": files_info.get("manual_chembl_file", "manual_chembl.csv"),
        },
    }

    # Add options if specified
    if "options" in app_config:
        config["options"] = app_config["options"]

    # Add custom column names if specified
    if "columns" in app_config:
        config["columns"] = app_config["columns"]

    return config


def run_pipeline_with_progress(app_config, files_info, work_dir, progress_placeholder, log_placeholder):
    """
    Run the DREXPA pipeline with progress tracking.

    Args:
        app_config: Streamlit configuration
        files_info: Dictionary with uploaded file paths
        work_dir: Working directory
        progress_placeholder: Streamlit placeholder for progress
        log_placeholder: Streamlit placeholder for logs

    Returns:
        tuple: (success: bool, results: dict or None, error_msg: str or None)
    """
    try:
        from drexpa import run_pipeline
        from drexpa.config import get_default_config

        # Build config
        config = build_config(app_config, files_info, work_dir)

        # Get selected steps
        steps_to_run = app_config.get("selected_steps", None)

        # Capture output
        output_buffer = StringIO()
        start_time = time.time()

        try:
            # Run pipeline with output capture
            with redirect_stdout(output_buffer):
                run_pipeline(
                    config_dict=config,
                    steps_to_run=steps_to_run,
                )

            elapsed_time = time.time() - start_time

            # Update logs
            logs = output_buffer.getvalue()
            with log_placeholder.container():
                st.text_area(
                    "Pipeline Logs",
                    logs,
                    height=200,
                    disabled=True
                )

            # Load results from output directory
            output_dir = config["global"]["output_dir"]
            results = load_pipeline_results(output_dir, elapsed_time)

            return True, results, None

        except Exception as e:
            elapsed_time = time.time() - start_time
            error_trace = traceback.format_exc()
            error_msg = str(e)

            with log_placeholder.container():
                st.text_area(
                    "Pipeline Logs & Error",
                    output_buffer.getvalue() + "\n\nERROR:\n" + error_trace,
                    height=300,
                    disabled=True
                )

            return False, None, error_msg

    except ImportError as e:
        return False, None, f"Failed to import DREXPA package: {str(e)}"


def load_pipeline_results(output_dir, elapsed_time):
    """
    Load pipeline output files into DataFrames.

    Args:
        output_dir: Path to pipeline output directory
        elapsed_time: Elapsed execution time

    Returns:
        dict: Results with DataFrames and metadata
    """
    results = {
        "execution_time": elapsed_time,
        "output_dir": output_dir,
    }

    output_path = Path(output_dir)
    if not output_path.exists():
        return results

    # Load key output files
    files_to_load = {
        "drug_profiles": "drug_profiles.csv",
        "drug_panel": "drug_panel_df.csv",
        "drug_targets": "drug_ChEMBL_targets.csv",
        "drug_node_targets": "drug_node_targets.csv",
        "drug_chembl_ids": "drug_ChEMBL_IDs.csv",
        "drug_doses": "drug_ChEMBL_doses.csv",
    }

    for key, filename in files_to_load.items():
        filepath = output_path / filename
        if filepath.exists():
            try:
                df = pd.read_csv(filepath)
                results[key] = df

                # Calculate statistics
                if key == "drug_profiles":
                    results["n_profiles"] = len(df)
                    results["n_drugs"] = df["drug_name"].nunique() if "drug_name" in df.columns else len(df)
                elif key == "drug_targets":
                    results["n_drugs"] = df["drug_name"].nunique() if "drug_name" in df.columns else len(df)

            except Exception as e:
                st.warning(f"Could not load {filename}: {str(e)}")

    # Count perturbation files
    try:
        perturbation_dir = output_path / "perturbations"
        if perturbation_dir.exists():
            perturbation_files = list(perturbation_dir.rglob("*"))
            results["n_perturbations"] = len([f for f in perturbation_files if f.is_file()])
    except:
        pass

    return results
