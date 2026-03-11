"""CELIOS pipeline runner.

This module implements a small, opinionated pipeline runner `run_celios`
"""

from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, Optional

from .base_defaults import CONFIG_TEMPLATE
from .features import node as node_mod, training as training_mod, files as files_mod, tissue as tissue_mod
import os
import traceback
import pandas as pd


def _deep_merge(target: Dict[str, Any], src: Dict[str, Any]) -> None:
    """Recursively merge src into target (in-place).

    Only dict values are merged recursively; other values overwrite.
    """
    for k, v in src.items():
        if k in target and isinstance(target[k], dict) and isinstance(v, dict):
            _deep_merge(target[k], v)
        else:
            target[k] = deepcopy(v)


def _build_config(user_config: Dict[str, Any]) -> Dict[str, Any]:
    """Merge user config with BASE_DEFAULTS and perform light validation.

    Raises ValueError when a minimal required structure is missing.
    """
    if not isinstance(user_config, dict):
        raise ValueError("Configuration must be a dictionary")

    # Start from the template defaults and overlay any user-provided values.
    final_paths = deepcopy(CONFIG_TEMPLATE.get("paths", {}))
    user_paths = user_config.get("paths") or {}
    if not isinstance(user_paths, dict):
        raise ValueError("'paths' must be a mapping if provided in config")
    final_paths.update(user_paths)

    final_general = deepcopy(CONFIG_TEMPLATE.get("general", {}))
    user_general = user_config.get("general") or {}
    if not isinstance(user_general, dict):
        raise ValueError("'general' must be a mapping if provided in config")
    final_general.update(user_general)

    final = {
        "paths": final_paths,
        "general": final_general,
        "steps": deepcopy(CONFIG_TEMPLATE.get("steps", {})),
    }

    # Apply user overrides: accept both legacy `advance` and `steps` keys
    user_steps = user_config.get("steps", {})
    if isinstance(user_steps, dict):
        _deep_merge(final["steps"], user_steps)

    advance = user_config.get("advance", {})
    if isinstance(advance, dict):
        _deep_merge(final["steps"], advance)

    # Run lightweight validation on the built config and attach diagnostics.
    try:
        # import validator lazily to avoid circular imports
        from .utils import validate as _validate
        # Use the simple path-only validator here to keep config building fast
        diagnostics = _validate.simple_validate_paths(final)
        # Attach diagnostics under a special key for callers to inspect
        final["_validation"] = diagnostics
    except Exception:
        # Do not fail config building for validation errors; attach error info
        try:
            final["_validation_error"] = traceback.format_exc()
        except Exception:
            final["_validation_error"] = "<failed to record validation error>"

    return final
    


def get_omics(
    config: Optional[Dict[str, Any]] = None,
    activity_file: Optional[str] = None,
    cell_line_file: Optional[str] = None,
    node_dict: Optional[dict] = None,
    tf_activity_file: Optional[str] = None,
    mutations_file: Optional[str] = None,
    cnv_file: Optional[str] = None,
    directory_output: Optional[str] = None,
    data_sources: Optional[list] = None,
    verbose: bool = True,
):
    """Convenience wrapper that calls `training.extract_omics`.

    Supports two usage modes:
      - `config` provided (recommended): reads Activity step from config.
      - `config` is None: caller supplies explicit activity arguments.

    Saving outputs and creating reports is triggered only when
    `directory_output` is provided (either explicitly or via config).
    """
    # If a config is provided, derive defaults from it
    if config is not None:
        cfg = _build_config(config)
        steps = cfg["steps"]
        act_cfg = steps.get("Activity", {}) or {}

        # resolve input base so relative filenames in the config point to the
        # `paths.input` directory (keeps behaviour similar to get_nodedict)
        base = Path(cfg.get("paths", {}).get("base", ".")).expanduser().resolve()
        input_path = Path(cfg.get("paths", {}).get("input", base / "input")).expanduser()

        # collect args from config when not explicitly provided
        activity_file = activity_file or act_cfg.get("activity_file")
        cell_line_file = cell_line_file or act_cfg.get("cell_line_file")
        tf_activity_file = tf_activity_file or act_cfg.get("tf_activity_file")
        mutations_file = mutations_file or act_cfg.get("mutations_file")
        cnv_file = cnv_file or act_cfg.get("cnv_file")
        node_dict = node_dict or act_cfg.get("node_dict")
        directory_output = directory_output or act_cfg.get("directory_output") or cfg["paths"].get("output")
        # If files are given as relative paths, assume they live under input_path
        def _resolve(p):
            if p is None:
                return None
            if isinstance(p, (str,)) and not Path(p).is_absolute():
                return str((input_path / p).resolve())
            return p

        activity_file = _resolve(activity_file)
        cell_line_file = _resolve(cell_line_file)
        tf_activity_file = _resolve(tf_activity_file)
        mutations_file = _resolve(mutations_file)
        cnv_file = _resolve(cnv_file)
        # Allow node_dict to be a path string; resolve relative to input_path and
        node_dict = _resolve(node_dict) if isinstance(node_dict, (str,)) else node_dict
        data_sources = data_sources or act_cfg.get("data_sources")
    else:
        data_sources = data_sources

    # Determine whether to save outputs / make reports based on directory presence
    do_save = directory_output is not None

    return training_mod.extract_omics(
        activity_file=activity_file,
        cell_line_file=cell_line_file,
        node_dict=node_dict,
        tf_activity_file=tf_activity_file,
        mutations_file=mutations_file,
        cnv_file=cnv_file,
        directory_output=directory_output,
        verbose=bool(verbose),
        data_sources=data_sources,
        save_master=do_save,
        make_report=do_save,
    )


def _validate_stop_after(stop_after: Optional[str]) -> None:
    if stop_after is None:
        return
    valid = ["fetch", "sif", "node_dict", "activity"]
    if stop_after not in valid:
        raise ValueError(f"Invalid stop_after '{stop_after}'. Valid: {valid}")

def get_nodedict(
    config: Optional[Dict[str, Any]] = None,
    node_input: object = None,
    hgnc_symbols_file: str = None,
    directory_output: str = None,
    include_alias_prev: Optional[bool] = None,
    plan: bool = False,
    verbose: bool = True,
    manual_symbols_file: Optional[str] = None,
) -> Dict[str, Any]:
    """Run pipeline up to node dictionary creation and return artifacts.

    The function supports two usage modes:
      - `config` provided: behavior unchanged (reads `steps -> Node` config).
      - `config` is None: the caller must supply `node_input` and any
        additional parameters via keyword arguments.

    Returns a dict containing the merged configuration (when available) and
    the resulting `node_dict`.
    """

    # Prepare mode-specific variables
    if config is not None:
        cfg = _build_config(config)
        paths = cfg["paths"]
        steps = cfg["steps"]

        base = Path(paths.get("base", ".")).expanduser().resolve()
        input_path = Path(paths.get("input", base / "input")).expanduser()

        # Ensure the input and (optional) output directories exist. The
        # legacy runs path has been removed; no default `runs` directory
        # will be created by the pipeline.
        output_path_config = paths.get("output")
        output = Path(output_path_config).expanduser() if output_path_config else None

        dirs = [input_path]
        if output:
            dirs.append(output)
        for p in dirs:
            p.mkdir(parents=True, exist_ok=True)
            if verbose:
                print(f"Ensured directory: {p}")

        if plan:
            return {"celios_configuration": cfg}

        artifacts: Dict[str, Any] = {"celios_configuration": cfg}
        config_mode = True

        node_cfg = steps.get("Node", {}) or {}

        # If caller provided node_input explicitly to the function, prefer it
        # Prefer the new, shorter key names when available but accept
        if node_input is None:
            # New preferred single-entry key
            if node_cfg.get("node_input") is not None:
                node_input = node_cfg.get("node_input")

        # Resolve relative paths in config: make strings absolute under input_path
        def _resolve(p):
            if p is None:
                return None
            if isinstance(p, (str,)) and not Path(p).is_absolute():
                return str((input_path / p).resolve())
            return p

        node_input = _resolve(node_input) if isinstance(node_input, (str,)) else node_input

        include_alias_config = node_cfg.get("include_alias_previous_symbols")
        if include_alias_config is None:
            include_alias_config = node_cfg.get("include_alias_prev", False)
        include_alias = bool(include_alias_config) if include_alias_prev is None else bool(include_alias_prev)

        hgnc_file_cfg = node_cfg.get("hgnc_symbols_file") if hgnc_symbols_file is None else hgnc_symbols_file
        hgnc_file_cfg = _resolve(hgnc_file_cfg) if isinstance(hgnc_file_cfg, (str,)) else hgnc_file_cfg
        dirout_cfg = str(output) if output else node_cfg.get("directory_output")
        manual_symbols_cfg = node_cfg.get("manual_symbols_file") if manual_symbols_file is None else manual_symbols_file
        manual_symbols_cfg = _resolve(manual_symbols_cfg) if isinstance(manual_symbols_cfg, (str,)) else manual_symbols_cfg

    else:
        # explicit-args mode
        if node_input is None:
            raise ValueError("Either 'config' or 'node_input' must be provided to run_node()")

        base = Path(".").expanduser().resolve()
        input_path = base / "input"
        output = Path(directory_output).expanduser() if directory_output else None

        if output:
            output.mkdir(parents=True, exist_ok=True)
            if verbose:
                print(f"Ensured directory: {output}")

        if plan:
            return {"celios_configuration": {}}

        artifacts = {"celios_configuration": {}}
        config_mode = False

        include_alias = bool(include_alias_prev) if include_alias_prev is not None else False
        hgnc_file_cfg = hgnc_symbols_file
        dirout_cfg = str(output) if output else directory_output
        manual_symbols_cfg = manual_symbols_file

        # Node step: call classmethods on Node
    try:
        # `run_celios()` prints the STEP header; avoid duplicating it here
        # Reject file-backed inputs only in explicit-args mode: callers must
        # provide objects (dataframes, lists). When running from a config the
        # code accepts string paths and will let `Node.from_object` handle them.
        if isinstance(node_input, str) and not config_mode:
            raise ValueError(
                "String/file-backed node_input is not supported in explicit-args mode. Load your file and pass the node object "
                "(e.g., dataframe or list) or provide a config that points to the file so get_nodedict() can build it."
            )

        if include_alias and not hgnc_file_cfg:
            raise ValueError("'hgnc_symbols_file' is required when 'include_alias_previous_symbols' is True for non-file inputs")

        # Prepare a small report structure to capture node-step metadata
        node_report = {
            'input': node_input,
            'hgnc_file': hgnc_file_cfg,
            'directory_output': dirout_cfg,
            'manual_symbols_file': manual_symbols_cfg,
            'used_sif_parser': False,
            'saved_path': None,
            'node_count': None,
            'missing_nodes': [],
        }

        # If node_input is a path to a SIF file, use the specialised constructor
        if isinstance(node_input, str) and Path(node_input).suffix.lower() == '.sif':
            node_report['used_sif_parser'] = True
            # Node.from_sif returns (node_instance, mapped_dict)
            node_instance, mapped = node_mod.Node.from_sif(
                sif_path=node_input,
                hgnc_symbols_file=hgnc_file_cfg,
                directory_output=dirout_cfg,
                include_alias_prev=include_alias,
                verbose=bool(verbose),
                manual_symbols_file=manual_symbols_cfg,
            )
            # Node.from_sif saves the node dict to directory_output when provided
            if dirout_cfg and hgnc_file_cfg:
                # Filename chosen by Node is 'node_HGNC_dict.csv' when HGNC mapping applied
                node_report['saved_path'] = os.path.join(dirout_cfg, 'node_HGNC_dict.csv')
        else:
            # Generic input: delegate to from_object which accepts lists, DataFrames, or file paths
            mapped = node_mod.Node.from_object(
                node_input=node_input,
                hgnc_symbols_file=hgnc_file_cfg,
                directory_output=dirout_cfg,
                include_alias_prev=include_alias,
                verbose=bool(verbose),
                manual_symbols_file=manual_symbols_cfg,
            )
            # If mapping was produced and we have a directory_output and hgnc file,
            # assume Node saved a HGNC node dict file.
            if dirout_cfg and hgnc_file_cfg:
                node_report['saved_path'] = os.path.join(dirout_cfg, 'node_HGNC_dict.csv')

        # Fill report metadata from mapped dict
        if isinstance(mapped, dict):
            node_report['node_count'] = len(mapped)
            node_report['missing_nodes'] = [n for n, syms in mapped.items() if not syms]

        artifacts["node_dict"] = mapped
        artifacts["node_report"] = node_report
    except Exception as e:
        raise RuntimeError(f"Error in Node step: {e}") from e

    return artifacts


#///////////////////////////////////////////////////////////////////////////////////////////////////////
# MAIN CELIOS PIPELINE 
#///////////////////////////////////////////////////////////////////////////////////////////////////////

def run_celios(
    config: Dict[str, Any],
    plan: bool = False,
    stop_after: Optional[str] = None,
    get_DLPfiles: bool = False,
    proliferation_state: bool = True,
    dna_damage_state: bool = False,
    cell_directory: Optional[str] = None,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Legacy full-run wrapper that uses `run_node` and `run_activity`.

    This keeps backward compatibility while delegating the two main phases to
    separate functions.
    """
    # validate
    _validate_stop_after(stop_after)

    if plan:
        return {"celios_configuration": _build_config(config)}

    artifacts: Dict[str, Any] = {"celios_configuration": _build_config(config)}

    # If a cell_directory is explicitly provided, treat that as an intent
    # to produce DL training files. Prefer the explicit `cell_directory`
    # when writing DL outputs below.
    if cell_directory and not get_DLPfiles:
        get_DLPfiles = True

    # If the config declares a `paths.cellfiles_dir`, treat that as an
    # indication the user wants DL training files produced by default.
    cfg_paths = artifacts["celios_configuration"].get("paths", {})
    cfg_cellfiles = cfg_paths.get("cellfiles_dir")
    if cfg_cellfiles and not get_DLPfiles:
        get_DLPfiles = True
        # Only set cell_directory if not explicitly provided
        if not cell_directory:
            cell_directory = cfg_cellfiles

    # 1) Get node dictionary
    cfg = artifacts["celios_configuration"]
    
    # Check if user explicitly defined Node step (not just from template defaults)
    user_steps = config.get("steps", {})
    node_step_user_defined = "Node" in user_steps
    
    act_cfg = cfg.get("steps", {}).get("Activity", {}) or {}
    node_dict_from_activity = act_cfg.get("node_dic") or act_cfg.get("node_dict")
    
    # Run Node step if user explicitly defined it in their config
    if node_step_user_defined:
        # Standard Node step: generate node dictionary
        if verbose:
            print("STEP 1: Node dictionary")
        node_artifacts = get_nodedict(config=config, verbose=verbose)
        # get_nodedict returns an artifacts dict (it may include the mapped node
        # dictionary under the 'node_dict' key). Extract the actual mapping so the
        # Activity step receives node->list mapping rather than the wrapper dict.
        if isinstance(node_artifacts, dict) and "node_dict" in node_artifacts:
            node_dict = node_artifacts.get("node_dict")
        else:
            node_dict = node_artifacts
        artifacts["node_dict"] = node_dict

        if stop_after == "node_dict":
            if verbose:
                print("Stopping after STEP: node_dict")
            return artifacts
    elif node_dict_from_activity:
        # Skip Node step and load the provided node_dict (only if Node step not user-defined)
        if verbose:
            print("STEP 1: Skipping node dictionary generation (already provided in Activity config)")
            print(f"Using node dictionary from: {node_dict_from_activity}")
        
        # Resolve relative path if necessary
        node_dict_path = Path(node_dict_from_activity)
        if isinstance(node_dict_from_activity, str) and not node_dict_path.is_absolute():
            # Resolve relative to input_path
            base = Path(cfg.get("paths", {}).get("base", ".")).expanduser().resolve()
            input_path = Path(cfg.get("paths", {}).get("input", base / "input")).expanduser()
            node_dict_path = input_path / node_dict_from_activity
        
        # Load node dictionary using utility function
        from .utils.io import load_node_dict_from_csv
        node_dict = load_node_dict_from_csv(str(node_dict_path), verbose=verbose)
        artifacts["node_dict"] = node_dict
        
        if stop_after == "node_dict":
            if verbose:
                print("Stopping after STEP: node_dict")
            return artifacts
    else:
        raise ValueError("No node dictionary source found. Either define 'steps.Node' in config or provide 'node_dic' in 'steps.Activity'")

    # 2) Extract omics / activity matrix
    if verbose:
        print("\n\nSTEP 2: Extracting omics - activity matrix")
    activity_df = get_omics(config=config, node_dict=node_dict, verbose=verbose)
    artifacts["activity_matrix"] = activity_df

    # 3) Optionally write DL training files
    if get_DLPfiles:
        cfg = artifacts["celios_configuration"]
        tissue_dir = cfg.get("paths", {}).get("tissue_dir")
        
        if tissue_dir:
            # Tissue-aware mode: organize by tissue/cell_line
            if verbose:
                print("\n\nSTEP 3: Writing DrugLogic pipeline training files (tissue-organized)")
                print(f"Using tissue directory: {tissue_dir}")
            
            # Parse cell_line_list.csv to build tissue mapping
            cell_line_file = cfg.get("steps", {}).get("Activity", {}).get("cell_line_file")
            if not cell_line_file:
                raise ValueError("cell_line_file must be specified in Activity config for tissue-aware output")
            
            # Resolve path
            base = Path(cfg.get("paths", {}).get("base", ".")).expanduser().resolve()
            input_path = Path(cfg.get("paths", {}).get("input", base / "input")).expanduser()
            if not Path(cell_line_file).is_absolute():
                cell_line_file = str(input_path / cell_line_file)
            
            artifacts["pipeline_dir"] = tissue_mod.write_tissue_files(
                activity_df, tissue_dir, cell_line_file, proliferation_state, dna_damage_state, verbose
            )
        
        else:
            # Legacy mode: use cellfiles_dir
            if verbose:
                print("\n\nSTEP 3: Writing DrugLogic pipeline training files")
                print(f"Using cell files directory: {cell_directory or cfg_cellfiles or 'default location'}")
            # prefer explicit dir, then pipeline config cellfiles_dir, then pipeline output path, then cwd
            out_base = cell_directory or cfg.get("paths", {}).get("cellfiles_dir") or cfg.get("paths", {}).get("output") or os.getcwd()
            # Avoid doubling 'cell_lines' when `cellfiles_dir` already points to such a folder
            if os.path.basename(str(out_base)).lower() == 'cell_lines':
                dl_dir = str(out_base)
            else:
                dl_dir = os.path.join(out_base, "cell_lines")
            os.makedirs(dl_dir, exist_ok=True)
            try:
                files_mod.DL_trainingfiles(activity_df, dl_dir, proliferation_state=proliferation_state, dna_damage_state=dna_damage_state,)
                artifacts["pipeline_dir"] = dl_dir
            except Exception as e:
                # don't fail the whole run on DL file writing; record error
                artifacts["pipeline_error"] = str(e)
            if verbose:
                # Count how many cell files were written (walk subdirectories)
                try:
                    num_files = 0
                    num_dirs = 0
                    for root, dirs, files in os.walk(dl_dir):
                        # count only actual files (not directories)
                        num_files += len([f for f in files if os.path.isfile(os.path.join(root, f))])
                        # count immediate subdirectories of the base as cell-line dirs
                        if os.path.abspath(root) == os.path.abspath(dl_dir):
                            num_dirs = len(dirs)
                    if verbose:
                        print(f"Wrote {num_files} training files across {num_dirs} cell-line directories to: {dl_dir}")
                except Exception:
                    pass
    if verbose:
        print("\nCELIOS pipeline completed successfully.")

    return artifacts

# celios/core.py