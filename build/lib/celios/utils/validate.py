"""Validation helpers for CELIOS configuration and input files.

This module exposes a focused helper `validate_config_inputs` that
resolves paths from a pipeline config and checks presence/quick-load of
important inputs (node files, hgnc file, activity files, cell-line file,
binary tables, TF file).

It accepts either a configuration dictionary or a path to a config file
(JSON or a Python file that exposes a `config` dict). A small CLI is
provided for convenience.
"""
from pathlib import Path
from typing import Dict, Any, Union
import traceback
import json
import importlib.util
import argparse

from ..features import node as node_mod


def _resolve(p, input_path: Path):
    if p is None:
        return None
    if isinstance(p, (str,)):
        p = Path(p)
        if not p.is_absolute():
            p = (input_path / p).resolve()
    return p


def _load_config_from_path(path: Union[str, Path]) -> Dict[str, Any]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f'Config file not found: {p}')

    # Try JSON first
    try:
        with open(p, 'r', encoding='utf-8') as fh:
            return json.load(fh)
    except Exception:
        pass

    # Fallback: try executing Python file and extracting `config` variable
    if p.suffix == '.py':
        spec = importlib.util.spec_from_file_location('celios_user_config', str(p))
        module = importlib.util.module_from_spec(spec)
        try:
            spec.loader.exec_module(module)  # type: ignore
            if hasattr(module, 'config'):
                return getattr(module, 'config')
        except Exception as e:
            raise RuntimeError(f'Failed to load config from {p}: {e}') from e

    raise ValueError('Unsupported config format. Provide a JSON file or a Python file exposing `config`')


def validate_config_inputs(user_config: Union[Dict[str, Any], str, Path]) -> Dict[str, Any]:
    """Validate and return a diagnostics dict for the provided config.

    `user_config` may be either a dict (already parsed) or a path to a
    config file (JSON or Python). The returned dict contains resolved paths
    and existence booleans and attempts lightweight loads where appropriate.
    It never writes files.
    """
    diagnostics: Dict[str, Any] = {}
    if isinstance(user_config, (str, Path)):
        user_config = _load_config_from_path(user_config)

    # Import _build_config lazily to avoid circular import when core imports this module
    from ..core import _build_config
    cfg = _build_config(user_config)
    diagnostics['resolved_config'] = cfg

    base = Path(cfg.get('paths', {}).get('base', '.')).expanduser().resolve()
    input_path = Path(cfg.get('paths', {}).get('input', base / 'input')).expanduser()
    diagnostics['base'] = str(base)
    diagnostics['input_path'] = str(input_path)

    steps = cfg.get('steps', {})
    node_cfg = steps.get('Node', {}) or {}
    act_cfg = steps.get('Activity', {}) or {}

    node_input = node_cfg.get('node_input')
    hgnc = node_cfg.get('hgnc_symbols_file')
    manual = node_cfg.get('manual_symbols_file')
    node_input_p = _resolve(node_input, input_path) if isinstance(node_input, (str,)) else node_input
    hgnc_p = _resolve(hgnc, input_path) if isinstance(hgnc, (str,)) else hgnc
    manual_p = _resolve(manual, input_path) if isinstance(manual, (str,)) else manual

    diagnostics['node_input'] = str(node_input_p) if node_input_p is not None else None
    diagnostics['hgnc_file'] = str(hgnc_p) if hgnc_p is not None else None
    diagnostics['manual_symbols_file'] = str(manual_p) if manual_p is not None else None
    diagnostics['node_input_exists'] = node_input_p.exists() if isinstance(node_input_p, Path) else None
    diagnostics['hgnc_exists'] = hgnc_p.exists() if isinstance(hgnc_p, Path) else None
    diagnostics['manual_symbols_exists'] = manual_p.exists() if isinstance(manual_p, Path) else None

    # Attempt node mapping if file provided
    try:
        if isinstance(node_input_p, Path) and node_input_p.exists():
            if node_input_p.suffix.lower() == '.sif':
                # from_sif returns (node_instance, mapped_dict)
                try:
                    # Pass manual_symbols_file through if provided
                    try:
                        _, mapped = node_mod.Node.from_sif(sif_path=str(node_input_p), hgnc_symbols_file=str(hgnc_p) if hgnc_p else None, verbose=False, manual_symbols_file=str(manual_p) if manual_p else None)
                    except TypeError:
                        # In case from_sif signature differs, fall back to from_object after extracting nodes
                        mapped = node_mod.Node.from_object(node_input=str(node_input_p), hgnc_symbols_file=str(hgnc_p) if hgnc_p else None, verbose=False, manual_symbols_file=str(manual_p) if manual_p else None)
            else:
                    mapped = node_mod.Node.from_object(node_input=str(node_input_p), hgnc_symbols_file=str(hgnc_p) if hgnc_p else None, verbose=False, manual_symbols_file=str(manual_p) if manual_p else None)
            diagnostics['node_map_count'] = len(mapped) if isinstance(mapped, dict) else None
            diagnostics['node_map_missing'] = [n for n, syms in (mapped or {}).items() if not syms]
    except Exception:
        diagnostics['node_map_error'] = traceback.format_exc()

    # Resolve activity inputs
    diagnostics['activity_file'] = str(_resolve(act_cfg.get('activity_file'), input_path)) if act_cfg.get('activity_file') else None
    diagnostics['cell_line_file'] = str(_resolve(act_cfg.get('cell_line_file'), input_path)) if act_cfg.get('cell_line_file') else None
    diagnostics['tf_activity_file'] = str(_resolve(act_cfg.get('tf_activity_file'), input_path)) if act_cfg.get('tf_activity_file') else None
    diagnostics['mutations_file'] = str(_resolve(act_cfg.get('mutations_file'), input_path)) if act_cfg.get('mutations_file') else None
    diagnostics['cnv_file'] = str(_resolve(act_cfg.get('cnv_file'), input_path)) if act_cfg.get('cnv_file') else None

    for k in ['activity_file', 'cell_line_file', 'tf_activity_file', 'mutations_file', 'cnv_file']:
        v = diagnostics.get(k)
        diagnostics[f'{k}_exists'] = Path(v).exists() if v else None

    # manual_symbols file existence flag already set above for node step

    return diagnostics


    # (validate_built_config removed — keep a single entrypoint and a
    # lightweight `simple_validate_paths` for config-build-time checks.)


def simple_validate_paths(cfg: Dict[str, Any]) -> Dict[str, Any]:
    """Lightweight validation: resolve paths from a built config and check existence.

    Does not attempt to parse or map files. Returns a dict with resolved
    paths and boolean flags indicating presence. Safe to run during config
    building and very fast.
    """
    result: Dict[str, Any] = {}
    base = Path(cfg.get('paths', {}).get('base', '.')).expanduser().resolve()
    input_path = Path(cfg.get('paths', {}).get('input', base / 'input')).expanduser()
    result['base'] = str(base)
    result['input_path'] = str(input_path)

    steps = cfg.get('steps', {})
    node_cfg = steps.get('Node', {}) or {}
    act_cfg = steps.get('Activity', {}) or {}

    def _res(p):
        if p is None:
            return None
        p = Path(p)
        if not p.is_absolute():
            p = (input_path / p).resolve()
        return str(p)

    # Node files
    node_input = node_cfg.get('node_input')
    hgnc = node_cfg.get('hgnc_symbols_file')
    node_input_p = _res(node_input) if isinstance(node_input, (str,)) else None
    hgnc_p = _res(hgnc) if isinstance(hgnc, (str,)) else None
    result['node_input'] = node_input_p
    result['hgnc_file'] = hgnc_p
    result['node_input_exists'] = bool(Path(node_input_p).exists()) if node_input_p else False
    result['hgnc_exists'] = bool(Path(hgnc_p).exists()) if hgnc_p else False

    # Activity files
    result['activity_file'] = _res(act_cfg.get('activity_file')) if act_cfg.get('activity_file') else None
    result['cell_line_file'] = _res(act_cfg.get('cell_line_file')) if act_cfg.get('cell_line_file') else None
    result['tf_activity_file'] = _res(act_cfg.get('tf_activity_file')) if act_cfg.get('tf_activity_file') else None
    result['mutations_file'] = _res(act_cfg.get('mutations_file')) if act_cfg.get('mutations_file') else None
    result['cnv_file'] = _res(act_cfg.get('cnv_file')) if act_cfg.get('cnv_file') else None

    for k in ['activity_file', 'cell_line_file', 'tf_activity_file', 'mutations_file', 'cnv_file']:
        v = result.get(k)
        result[f'{k}_exists'] = bool(Path(v).exists()) if v else False

    return result


def _cli():
    parser = argparse.ArgumentParser(description='Validate CELIOS pipeline config and inputs')
    parser.add_argument('--config', '-c', help='Path to pipeline config (JSON or Python file)')
    parser.add_argument('--print-json', action='store_true', help='Print diagnostics as JSON')
    args = parser.parse_args()

    if not args.config:
        parser.error('Please provide --config <path>')

    try:
        diag = validate_config_inputs(args.config)
        if args.print_json:
            print(json.dumps(diag, indent=2))
        else:
            # Simple printable summary
            print('Validation diagnostics:')
            print(f"Node input: {diag.get('node_input')} (exists: {diag.get('node_input_exists')})")
            print(f"HGNC file: {diag.get('hgnc_file')} (exists: {diag.get('hgnc_exists')})")
            print(f"Node map count: {diag.get('node_map_count')}")
            print(f"Activity file: {diag.get('activity_file')} (exists: {diag.get('activity_file_exists')})")
    except Exception:
        print('Validation failed:')
        traceback.print_exc()


if __name__ == '__main__':
    _cli()
