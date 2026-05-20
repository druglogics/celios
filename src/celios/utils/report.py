"""Reporting utilities for CELIOS.

This module contains two helpers:
- `write_run_report` (migrated from `features/report.py`) which writes a run
  report to disk given a metadata dict.
- `print_training_report` which encapsulates the simple console reporting
  previously present in `features/files.py`.
"""
from typing import Dict, Any
import os
import logging

logger = logging.getLogger(__name__)

# In-memory run log collector. Other modules should call `add_log()` to
# capture verbose messages that should appear in the final `run_log.txt`.
_RUN_MESSAGES = []


def add_log(msg: str) -> None:
    """Append a single-line message to the in-memory run log.

    Keep this lightweight: modules can call it unconditionally; the
    collector will be serialized into `run_log.txt` by
    `activitymatrix_report`.
    """
    try:
        _RUN_MESSAGES.append(str(msg))
    except Exception:
        # Be defensive; logging should not break the pipeline
        try:
            logger.debug('Failed to append run log message')
        except Exception:
            pass


def add_resolution_report(report: dict) -> None:
    """Append a formatted cell line resolution report to the run log.
    
    Args:
        report: dict with keys like 'total_rows', 'resolved', 'ambiguous', 'unresolved'
    """
    try:
        add_log("-------------------------------------------------------------------")
        add_log("CELL LINE IDENTIFIER RESOLUTION REPORT")
        add_log("-------------------------------------------------------------------")
        add_log(f"Total rows processed: {report.get('total_rows', 0)}")
        add_log(f"Successfully resolved: {report.get('resolved', 0)}")
        add_log(f"Ambiguous (multiple matches): {report.get('ambiguous', 0)}")
        add_log(f"Unresolved: {report.get('unresolved', 0)}")
        
        cache_stats = report.get('cache_stats', {})
        if cache_stats.get('cache_enabled'):
            add_log(f"Cache: {cache_stats.get('cached_identifiers', 0)} cached identifiers")
    except Exception:
        logger.debug('Failed to add resolution report to log')


def get_run_logs() -> list:
    return list(_RUN_MESSAGES)


def clear_run_logs() -> None:
    try:
        _RUN_MESSAGES.clear()
    except Exception:
        pass


def activitymatrix_report(directory: str, data: Dict[str, Any], verbose: bool = False) -> str:
    os.makedirs(directory, exist_ok=True)
    report_lines = []
    report_lines.append("------------------------------------------------------")
    report_lines.append("ACTIVITY MATRIX RUN REPORT")
    report_lines.append("------------------------------------------------------")
    report_lines.append("")
    report_lines.append("VARIABLES")
    report_lines.append("")
    report_lines.append(f"Selected sources: {data.get('selected_sources')}")

    cmd = data.get('command')
    if cmd:
        report_lines.append("")
        report_lines.append(f"Command: {cmd}")
    rt = data.get('runtime_seconds')
    if rt is not None:
        report_lines.append(f"Run time (seconds): {rt:.2f}")

    node_dict = data.get('node_dict')
    if node_dict is None:
        report_lines.append("")
        report_lines.append("Node dictionary: NOT PROVIDED")
    else:
        # Include node step metadata if available
        node_report = data.get('node_report')
        if node_report:
            report_lines.append("")
            report_lines.append("NODE STEP LOG")
            report_lines.append("")
            report_lines.append(f"Input: {node_report.get('input')}")
            report_lines.append(f"HGNC file: {node_report.get('hgnc_file')}")
            report_lines.append(f"Directory output: {node_report.get('directory_output')}")
            report_lines.append(f"Used SIF parser: {node_report.get('used_sif_parser')}")
            if node_report.get('saved_path'):
                report_lines.append(f"Saved node dict: {node_report.get('saved_path')}")
            
            # Include timing summary if available
            if node_report and node_report.get('timing'):
                timing = node_report.get('timing')
                report_lines.append("")
                report_lines.append("NODE STEP TIMING SUMMARY")
                report_lines.append("")
                if timing.get('total') is not None:
                    report_lines.append(f"Total node step time (s): {timing.get('total'):.3f}")
                if timing.get('hgnc_load_time') is not None:
                    report_lines.append(f"HGNC load time (s): {timing.get('hgnc_load_time'):.3f}")
                if timing.get('lookup_build_time') is not None:
                    report_lines.append(f"Lookup build time (s): {timing.get('lookup_build_time'):.3f}")
                per = timing.get('per_node', {})
                if per:
                    report_lines.append(f"Per-node mapping (s): mean={per.get('mean',0):.4f} median={per.get('median',0):.4f} max={per.get('max',0):.4f}")
                if timing.get('alias_expansion_time') is not None:
                    report_lines.append(f"Alias/prev expansion time (s): {timing.get('alias_expansion_time'):.4f}")
            report_lines.append("")
        
        report_lines.append("")
        report_lines.append("NODE DICTIONARY")
        report_lines.append("")
        report_lines.append(f"Length of node dictionary: {len(node_dict)}")
        report_lines.append("")
        report_lines.append("Node dictionary entries:")
        for node, syms in node_dict.items():
            report_lines.append(f"{node}: {syms}")
        report_lines.append("")
        missing = data.get('missing_hgnc_nodes', [])
        report_lines.append(f"Nodes with missing HGNC_symbol entries: {len(missing)}")
        if missing:
            report_lines.append(f"Examples: {missing[:20]}")

    sidm_list = data.get('sidm_list')
    sidm_dict = data.get('sidm_dict')
    if sidm_list is not None:
        report_lines.append("")
        report_lines.append("CELL LINES / SIDM LIST")
        report_lines.append("")
        report_lines.append(f"SIDM count: {len(sidm_list)}")
        report_lines.append("")
        report_lines.append("SIDM list (ordered):")
        for sidm in sidm_list:
            report_lines.append(sidm)
        report_lines.append("")
        if sidm_dict:
            report_lines.append("Cell line names (SIDM -> name):")
            for sidm, name in sidm_dict.items():
                report_lines.append(f"{sidm}: {name}")

    report_lines.append("")
    report_lines.append("TRAINING TRAINING DATA")
    report_lines.append("")
    for key in ['raw_shape', 'prep_shape', 'norm_shape', 'node_expr_shape', 'muts_shape', 'cnv_shape', 'tfm_shape']:
        val = data.get(key)
        if val is not None:
            report_lines.append(f"{key}: {val}")

    master_fp = data.get('master_fp')
    master_shape = data.get('master_shape')
    if master_fp:
        report_lines.append("")
        report_lines.append(f"Master matrix saved: {master_fp} (shape={master_shape})")

    final_fp = data.get('final_fp')
    final_shape = data.get('final_shape')
    if final_fp:
        report_lines.append(f"Final activity saved: {final_fp} (shape={final_shape})")

    nodes_with_all_missing = data.get('nodes_with_all_missing', [])
    report_lines.append("")
    report_lines.append(f"Nodes with all-missing values in final matrix: {len(nodes_with_all_missing)}")
    if nodes_with_all_missing:
        report_lines.append(f"Examples: {nodes_with_all_missing[:20]}")

    report_path = os.path.join(directory, 'run_log.txt')
    with open(report_path, 'w', encoding='utf-8') as fh:
        # Insert any collected run messages before writing the report
        run_msgs = get_run_logs()
        if run_msgs:
            report_lines.insert(0, "VERBOSE LOG MESSAGES")
            report_lines[0:0] = [m for m in run_msgs]
        fh.write('\n'.join(report_lines))

    # Clear run messages after writing the report to avoid duplication
    clear_run_logs()

    if verbose:
        for l in report_lines:
            logger.info(l)

    return report_path


def cellfiles_report(activity_matrix, directory_cell_lines, proliferation_state=True, dna_damage_state=False):
    """Capture a short training-file creation report into the run log.

    Instead of printing to stdout, append a small human-friendly summary
    into the in-memory run-log collector used by `activitymatrix_report`.
    """
    try:
        add_log('-------------------------------------------------------------------')
        add_log('REPORT')
        add_log('Running training file writer')
        add_log('-------------------------------------------------------------------')
        add_log('Creating training files for each cell line')
        # Try to capture basic activity_matrix info
        try:
            cols = getattr(activity_matrix, 'columns', '<unknown>')
            add_log(f'Activity frame columns: {cols}')
            try:
                preview = None
                if hasattr(activity_matrix, 'head'):
                    try:
                        preview = activity_matrix.head(3)
                    except Exception:
                        preview = None
                if preview is not None:
                    add_log('Activity frame preview:\n' + str(preview))
                else:
                    add_log('<unable to display activity_matrix preview>')
            except Exception:
                add_log('<unable to display activity_matrix preview>')
        except Exception:
            add_log('<unable to inspect activity_matrix>')

        add_log(f'Adding proliferation state: {proliferation_state}')
        add_log(f'Adding DNA damage state: {dna_damage_state}')
    except Exception:
        # Never raise from reporting helpers
        logger.debug('Failed to capture cellfiles report')