import os
import pandas as pd
from pathlib import Path

from celios.features.training import extract_omics, resolve_cell_lines


def test_training_import_and_api():
    # Ensure public API exists
    assert callable(extract_omics)
    assert callable(resolve_cell_lines)


def test_extract_omics_runs_with_pre_resolved(tmp_path):
    # Prepare minimal expression CSV with symbols x SIDMs
    activity_csv = tmp_path / "activity.csv"
    df = pd.DataFrame({
        'symbol': ['GENE1', 'GENE2'],
        'SIDM1': [10.0, 0.0],
        'SIDM2': [5.0, 2.0],
    }).set_index('symbol')
    df.to_csv(activity_csv)

    # Minimal node dict mapping genes to nodes
    node_dict = {'NodeA': ['GENE1'], 'NodeB': ['GENE2']}

    # Pre-resolved cell-line mapping to avoid network resolution in tests
    pre_resolved = {
        'sidm_list': ['SIDM1', 'SIDM2'],
        'sidm_dict': {'SIDM1': 'CELL1', 'SIDM2': 'CELL2'},
        'alias_to_sidm': {},
    }

    out_dir = tmp_path / 'out'
    out_dir.mkdir()

    final = extract_omics(
        activity_file=str(activity_csv),
        cell_line_file=None,
        node_dict=node_dict,
        tf_activity_file=None,
        mutations_file=None,
        cnv_file=None,
        directory_output=str(out_dir),
        verbose=False,
        data_sources=['expression'],
        save_master=True,
        make_report=False,
        pre_resolved_cell_lines=pre_resolved,
    )

    # Check that outputs were written
    master_fp = out_dir / 'activity_master_matrix.csv'
    final_fp = out_dir / 'activity_from_master.csv'
    assert master_fp.exists(), f"Expected master file at {master_fp}"
    assert final_fp.exists(), f"Expected final file at {final_fp}"

    # Load final and do basic shape checks
    loaded_final = pd.read_csv(final_fp, index_col=0)
    assert loaded_final.shape[0] >= 1
    assert 'CELL1' in loaded_final.columns or 'SIDM1' in loaded_final.columns


def test_format_override_arg_accepts_value(tmp_path):
    # Verify the function accepts format_override without raising
    activity_csv = tmp_path / "activity2.csv"
    pd.DataFrame({'symbol': ['G'], 'SIDM1': [1.0]}).set_index('symbol').to_csv(activity_csv)

    node_dict = {'N': ['G']}
    pre_resolved = {'sidm_list': ['SIDM1'], 'sidm_dict': {'SIDM1': 'C1'}, 'alias_to_sidm': {}}

    # Should not raise
    extract_omics(
        activity_file=str(activity_csv),
        cell_line_file=None,
        node_dict=node_dict,
        directory_output=str(tmp_path / 'out2'),
        verbose=False,
        data_sources=['expression'],
        format_override='old',
        pre_resolved_cell_lines=pre_resolved,
        save_master=False,
        make_report=False,
    )
