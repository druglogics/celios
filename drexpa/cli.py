"""
DREXPA Pipeline Command Line Interface
"""
import argparse
import logging
import sys
from .main import run_pipeline
from .config import get_default_config, merge_config
from .. import __version__
from .step_registry import ordered_step_names, until_choices


def setup_logging(verbose=False):
    """Configure CLI logging format and level."""
    logging.basicConfig(
        level=logging.INFO if verbose else logging.WARNING,
        format="%(asctime)s %(levelname)s %(name)s %(message)s"
    )


def main():
    """Main CLI entry point."""
    available_steps = ordered_step_names()

    parser = argparse.ArgumentParser(
        description="Run DREXPA drug panel generation pipeline. "
                   "Configure data file paths in config.py or use --config option.",
        epilog="""
Examples:
  python -m drexpa                                    # Run full pipeline with default config
  python -m drexpa --until profiles                   # Run until drug profiles
  python -m drexpa --steps profiles,panel             # Run only profiles and panel steps
  python -m drexpa --synergy-data my_data.csv         # Use custom synergy data file
  python -m drexpa --config my_config.json            # Use custom configuration file
  python -m drexpa --output-dir ./results --verbose   # Custom output directory with verbose output

Note: Modify paths in drexpa/config.py or provide a custom config file to point to your data files.
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '--version', 
        action='version', 
        version=f'DREXPA {__version__}'
    )
    
    parser.add_argument(
        '--synergy-data', 
        type=str, 
        help='Path to synergy data CSV file'
    )
    
    parser.add_argument(
        '--config', 
        type=str, 
        help='Path to custom config JSON file (optional)'
    )
    
    parser.add_argument(
        '--output-dir', 
        type=str, 
        help='Output directory (overrides config)'
    )
    
    parser.add_argument(
        '--verbose', 
        action='store_true', 
        help='Enable verbose output'
    )
    
    parser.add_argument(
        '--until', 
        type=str,
           choices=until_choices(),
        help='Run pipeline until specified step (inclusive). Runs all prerequisite steps. '
               f'Options: {", ".join(until_choices())}'
    )
    
    parser.add_argument(
        '--steps', 
        type=str,
        help='Comma-separated list of specific steps to run. Automatically includes prerequisites. '
             f'Options: {", ".join(available_steps)}. Example: --steps profiles,panel'
    )
    
    args = parser.parse_args()
    
    # Load default config
    config_dict = get_default_config()
    
    # Override with command line args
    if args.output_dir:
        config_dict['global']['output_dir'] = args.output_dir
    
    if args.verbose:
        config_dict['global']['verbose'] = True

    setup_logging(args.verbose)
    
    # Load custom config if provided
    if args.config:
        import json
        with open(args.config, 'r') as f:
            custom_config = json.load(f)
        config_dict = merge_config(config_dict, custom_config)
    
    # Determine which steps to run
    steps_to_run = None
    if args.until:
        steps_to_run = f"until_{args.until}"
    elif args.steps:
        steps_to_run = [step.strip() for step in args.steps.split(',')]
    
    # Run pipeline
    try:
        run_pipeline(config_dict, args.synergy_data, steps_to_run)
    except Exception as e:
        print(f"Error running pipeline: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()