"""Command-line interface for CELIOS.

Supports two modes:
- `run`: run the full pipeline via `celios.core.run_celios`
- `node-from-sif` / `node-from-object`: run the Node feature helpers

Usage examples:
  python -m celios.cli run --config config.yaml
  python -m celios.cli node-from-sif --sif examples/DNAdamage.sif --hgnc configs/hgnc_complete_set.txt
  python -m celios.cli node-from-object --input nodes.txt --hgnc configs/hgnc_complete_set.txt
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any

try:
    import yaml  # optional
except Exception:
    yaml = None


def _read_config(path: str) -> Any:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Config file not found: {p}")
    if p.suffix.lower() in (".yml", ".yaml"):
        if yaml is None:
            raise RuntimeError("PyYAML is not installed. Install pyyaml to use YAML config files.")
        return yaml.safe_load(p.read_text())
    else:
        # assume JSON
        return json.loads(p.read_text())


def _save_node_dict(node_dict: dict, out_path: Path) -> None:
    import pandas as pd

    items = [(k, ", ".join(v) if isinstance(v, (list, tuple)) else v) for k, v in node_dict.items()]
    df = pd.DataFrame(items, columns=["node_name", "symbols"])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)


def main(argv: list[str] | None = None) -> int:
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser(prog="celios", description="CELIOS pipeline and feature runner")
    sub = parser.add_subparsers(dest="command")

    # run full pipeline
    p_run = sub.add_parser("run", help="Run the full CELIOS pipeline")
    p_run.add_argument("--config", required=True, help="Path to pipeline config (JSON or YAML)")
    p_run.add_argument("--plan", action="store_true", help="Only print plan and exit")
    p_run.add_argument("--stop-after", default=None, help="Stop after step")
    p_run.add_argument("--verbose", action="store_true", help="Verbose output")

    # node-from-sif
    p_sif = sub.add_parser("node-from-sif", help="Run Node.from_sif() to extract and map nodes from a SIF file")
    p_sif.add_argument("--sif", required=True, help="Path to SIF file")
    p_sif.add_argument("--hgnc", required=True, help="Path to HGNC symbols file (tab-separated) for mapping")
    p_sif.add_argument("--out", default=None, help="Optional output CSV file to save mapped node dictionary")
    p_sif.add_argument("--include_alias_prev", action="store_true", dest="include_alias_prev",
                         help="Include alias and previous HGNC symbols when mapping (default: False)")
    p_sif.add_argument("--verbose", action="store_true")

    # node-from-object
    p_obj = sub.add_parser("node-from-object", help="Run Node.from_object() from a list, file, or DataFrame")
    p_obj.add_argument("--input", required=True, help="Input list/file (comma-separated list or path to file) or path to CSV/TSV")
    p_obj.add_argument("--hgnc", default=None, help="Optional HGNC symbols file for mapping")
    p_obj.add_argument("--out", default=None, help="Optional output CSV file to save mapped node dictionary")
    p_obj.add_argument("--include_alias_prev", action="store_true", dest="include_alias_prev",
                       help="Include alias and previous HGNC symbols when mapping (default: False)")
    p_obj.add_argument("--verbose", action="store_true")

    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return 1

    if args.command == "run":
        from .core import run_celios

        cfg = _read_config(args.config)
        artifacts = run_celios(cfg, plan=args.plan, stop_after=args.stop_after, verbose=args.verbose)
        print("Run completed. Produced artifacts keys:", list(artifacts.keys()))
        return 0

    if args.command == "node-from-sif":
        from .features.node import Node

        node, mapped = Node.from_sif(sif_path=args.sif, hgnc_symbols_file=args.hgnc, directory_output=None,
                         include_alias_prev=bool(args.include_alias_prev), verbose=args.verbose)
        print("Number of nodes mapped:", len(mapped) if mapped else 0)
        if args.out:
            out_path = Path(args.out)
            _save_node_dict(mapped, out_path)
            print(f"Saved mapped node dict to: {out_path}")
        return 0

    if args.command == "node-from-object":
        from .features.node import Node
        inp = args.input
        # If input looks like a comma-separated list, split it; otherwise pass as path
        if "," in inp and not Path(inp).exists():
            node_input = [s.strip() for s in inp.split(",") if s.strip()]
        else:
            node_input = inp

        mapped = Node.from_object(node_input=node_input, hgnc_symbols_file=args.hgnc, directory_output=None,
                      include_alias_prev=bool(args.include_alias_prev), verbose=args.verbose)
        print("Number of nodes mapped:", len(mapped) if mapped else 0)
        if args.out:
            _save_node_dict(mapped, Path(args.out))
            print(f"Saved mapped node dict to: {args.out}")
        return 0

    parser.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
