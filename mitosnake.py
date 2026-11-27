#!/usr/bin/env python3

import argparse
import subprocess
import sys
import os
import logging
from pathlib import Path
import yaml
import tempfile
import shutil
from datetime import datetime

# Default configuration
DEFAULT_CONFIG = {
    "Resources": {
        "default": {
            "threads": 4,
            "mem_mb": 65536,
            "disk_mb": 51200
        },
        "high": {
            "threads": 8,
            "mem_mb": 131072,
            "disk_mb": 102400
        },
        "ultra": {
            "threads": 16,
            "mem_mb": 262144,  # 256GB
            "disk_mb": 204800   # 200GB
        }
    },
    "fastp": {
        "cut_front": True,
        "cut_front_window_size": 5,
        "cut_front_mean_quality": 25,
        "average_qual": 25,
        "length_required": 45
    },
    "subsample": {
        "short_cov": 2.0,
        "long_cov": 10.0,
        "genome_size": 400,  # MB
        "min_length": 10000  # Minimum length for subsampling
    },
    "filter": {
        "min_match": 2000,     # Minimum number of matching bases
        "min_coverage": 0.3    # Minimum coverage ratio requirement
    },
    # Seed segmentation parameters
    "seed_segment_length": 25000,  # Length of seed segments
    "seed_overlap": 5000,          # Overlap between seed segments
    "min_segment_length": 10000,   # Minimum segment length
    
    # Deep learning filter configuration
    "use_dl_filter": True,  # Use deep learning model by default
    
    # Iterative assembly configuration
    "use_iterative": True,  # Use iterative assembly by default
    "keep_intermediate_assemblies": False,  # Don't keep intermediate results by default
    "iterative_rounds": 2  # Default number of iterative rounds
}

__version__ = "0.9.0"

def parse_arguments():
    parser = argparse.ArgumentParser(description="MitoSnake: Pipeline for plant mitochondrial genome assembly")
    
    # Required arguments
    parser.add_argument("--module", choices=["st", "ht"], default="ht",
                      help="Module to run: st (short-read) or ht (hybrid)")
    
    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--sample-id", help="Sample identifier")
    input_group.add_argument("--sample-list", help="Path to sample list file")
    
    parser.add_argument("--read1", help="Path to forward reads")
    parser.add_argument("--read2", help="Path to reverse reads")
    parser.add_argument("--long-read", help="Path to long reads (required for hybrid mode)")
    parser.add_argument("--longread-type", choices=["clr", "hifi", "ont"],
                      help="Type of long reads: clr (PacBio CLR), hifi (PacBio HiFi), ont (Oxford Nanopore)")

    # Directory options
    parser.add_argument("--working-dir", default=os.getcwd(),
                      help="Working directory (default: current directory)")
    
    # Resource options
    parser.add_argument("--resource-mode", choices=["default", "high", "ultra"],
                      default="default", help="Resource usage mode (ultra for >3GB genomes)")
    parser.add_argument("--cores", type=int, default=4,
                      help="Number of cores (default: 4)")
    parser.add_argument("--jobs", type=int, default=1,
                      help="Number of parallel jobs (default: 1)")
    
    # Subsampling options
    parser.add_argument("--short-cov", type=float, default=2.0,
                      help="Target coverage for short reads (default: 2.0)")
    parser.add_argument("--long-cov", type=float, default=10.0,
                      help="Target coverage for long reads (default: 10.0)")
    parser.add_argument("--genome-size", type=int, default=400,
                      help="Estimated genome size in MB (default: 400)")
    parser.add_argument("--min-length", type=int, default=10000,
                      help="Minimum sequence length for subsampling (default: 10000)")
    
    # Filtering options 
    parser.add_argument("--min-match", type=int, default=2000,
                      help="Minimum number of matching bases for filtering (default: 2000)")
    parser.add_argument("--min-coverage", type=float, default=0.3,
                      help="Minimum coverage (matches/query_length) for filtering (default: 0.3)")
    
    # skip deep learning filter
    parser.add_argument("--skip-dl-filter", action="store_true",
                      help="Skip deep learning model (Tiara) for sequence filtering")
    
    # Workflow control
    parser.add_argument("--profile", help="Snakemake profile for cluster execution")
    parser.add_argument("--unlock", action="store_true",
                      help="Remove lock on working directory (recommended)")
    parser.add_argument("--clean", action="store_true",
                      help="Remove .snakemake directory (use with caution)")
    parser.add_argument("--dry-run", action="store_true",
                      help="Perform dry run")
    parser.add_argument("--verbose", action="store_true",
                      help="Verbose output")
    
    # Version argument
    parser.add_argument("--version", action="version", 
                      version=f"MitoSnake version {__version__}")
    
    # New parameter
    parser.add_argument("--conda-path", 
                       help="Path to conda installation (for example, a Miniconda3 prefix)")
    
    parser.add_argument("--keep-going", action="store_true",
                      help="Continue with independent jobs if some jobs fail")
    
    # Added new parameter
    parser.add_argument("--generate-resource-summary", action="store_true",
                      help="Generate resource usage summary (for advanced users)")
    
    parser.set_defaults(
        mito_seed=None,
        skip_iterative=False,
        keep_intermediate=False,
        iterative_rounds=2,
        seed_segment_length=25000,
        seed_overlap=5000,
        min_segment_length=10000,
    )

    args = parser.parse_args()
    
    # Validate arguments
    if args.module == "ht" and args.sample_id and not args.long_read:
        parser.error("--long-read is required for hybrid mode (--module ht)")
    
    if args.sample_id:
        if args.module in ["st", "ht"] and not (args.read1 and args.read2):
            parser.error("--read1 and --read2 are required for short-read (st) and hybrid (ht) modes")
        elif args.module == "lt-seed":
            if not args.mito_seed:
                parser.error("--mito-seed is required for seed-guided (lt-seed) mode")
            if not args.long_read:
                parser.error("--long-read is required for seed-guided (lt-seed) mode")
            if not args.longread_type:
                parser.error("--longread-type is required for seed-guided (lt-seed) mode")
    
    if args.unlock and args.clean:
        parser.error("--unlock and --clean cannot be used together")
    
    return args

def load_config(args):
    """Load and merge configuration from defaults and command line."""
    config = DEFAULT_CONFIG.copy()
    
    # Update resource mode
    if args.resource_mode:
        config["resource_mode"] = args.resource_mode
    
    # Update subsampling parameters
    if args.short_cov:
        config["subsample"]["short_cov"] = args.short_cov
    if args.long_cov:
        config["subsample"]["long_cov"] = args.long_cov
    if args.genome_size:
        config["subsample"]["genome_size"] = args.genome_size
    if args.min_length:
        config["subsample"]["min_length"] = args.min_length
    
    # Update filtering parameters
    if args.min_match:
        config["filter"]["min_match"] = args.min_match
    if args.min_coverage:
        config["filter"]["min_coverage"] = args.min_coverage
    
    # Update lt-seed mode parameters
    if args.module == "lt-seed":
        config["seed_segment_length"] = args.seed_segment_length
        config["seed_overlap"] = args.seed_overlap
        config["min_segment_length"] = args.min_segment_length
    
    # Update deep learning filter parameter
    if args.skip_dl_filter:
        config["use_dl_filter"] = False
    
    # Update iterative assembly parameters
    if args.skip_iterative:
        config["use_iterative"] = False
    if args.keep_intermediate:
        config["keep_intermediate_assemblies"] = True
    
    # Update generate_resource_summary parameter
    if args.generate_resource_summary:
        config["generate_resource_summary"] = True
    
    # Update iterative rounds parameter
    if args.iterative_rounds:
        config["iterative_rounds"] = args.iterative_rounds
    
    return config

def setup_logging(verbose):
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format='%(asctime)s - %(levelname)s - %(message)s')

def clean_snakemake_dir(working_dir):
    """Remove .snakemake directory."""
    snakemake_dir = os.path.join(working_dir, ".snakemake")
    if os.path.exists(snakemake_dir):
        try:
            shutil.rmtree(snakemake_dir)
            logging.info(f"Successfully removed {snakemake_dir}")
        except Exception as e:
            logging.error(f"Failed to remove {snakemake_dir}: {e}")
            sys.exit(1)

def create_sample_directories(working_dir, sample_id):
    """Create necessary directories for a sample."""
    directories = [
        os.path.join(working_dir, sample_id, "subsample"),
        os.path.join(working_dir, sample_id, "assembly"),
        os.path.join(working_dir, sample_id, "logs")
    ]
    for directory in directories:
        os.makedirs(directory, exist_ok=True)
        logging.debug(f"Created directory: {directory}")

def get_conda_info(custom_conda_path=None):
    """Get conda related path information
    
    Args:
        custom_conda_path: Custom conda installation path
    Returns:
        conda_base: conda base path
        envs_dir: conda environment directory
    """
    if custom_conda_path:
        conda_base = custom_conda_path
        envs_dir = os.path.join(conda_base, "envs")
        # set necessary environment variables
        os.environ["CONDA_PREFIX"] = conda_base
        os.environ["CONDA_EXE"] = os.path.join(conda_base, "bin", "conda")
        logging.info(f"Using custom conda path: {conda_base}")
        return conda_base, envs_dir
    
    try:
        # First try to get from environment variables
        conda_base = os.environ.get('CONDA_PREFIX')
        if conda_base and 'envs' in conda_base:
            conda_base = os.path.dirname(os.path.dirname(conda_base))
        
        if not conda_base:
            # If environment variables don't exist, try executing command
            conda_base_cmd = ["conda", "info", "--base"]
            conda_base = subprocess.check_output(
                conda_base_cmd,
                text=True,
                env=os.environ.copy()
            ).strip()
        
        return conda_base, os.path.join(conda_base, "envs")
    except Exception as e:
        logging.error(f"Failed to get conda base path: {e}")
        # Final fallback solution
        conda_base = os.path.dirname(os.path.dirname(sys.executable))
        logging.warning(f"Falling back to conda base path: {conda_base}")
        return conda_base, os.path.join(conda_base, "envs")

def main():
    """Main function to run the pipeline."""
    try:
        args = parse_arguments()
        setup_logging(args.verbose)
        
        logging.info("Starting MitoSnake pipeline...")
        logging.info(f"Working directory: {args.working_dir}")
        logging.info(f"Module: {args.module}")
        
        # Create working directory
        os.makedirs(args.working_dir, exist_ok=True)
        os.makedirs(os.path.join(args.working_dir, "tmp"), exist_ok=True)
        
        # Initialize configuration
        config = load_config(args)
        
        # Get conda information
        conda_base, envs_dir = get_conda_info(args.conda_path)
        
        # Update configuration
        config["conda_base"] = conda_base
        config["conda_envs_dir"] = envs_dir
        
        # add tiara_conda_env to configuration
        config.update({
            "tiara_conda_env": "mitosnake_tiara"  # ensure using correct tiara environment
        })
        
        if args.sample_id:
            logging.info(f"Sample ID: {args.sample_id}")
            create_sample_directories(args.working_dir, args.sample_id)
        elif args.sample_list:
            logging.info(f"Using sample list: {args.sample_list}")
            # Read sample list file
            sample_ids = []
            samples_info = {}
            
            with open(args.sample_list) as f:
                # Read header line
                header = next(f).strip().split('\t')
                header_map = {col: idx for idx, col in enumerate(header)}

                # Check required columns
                if args.module in ["st", "ht"]:
                    required_columns = ['sample_id', 'read1', 'read2']
                elif args.module == "lt-seed":
                    required_columns = ['sample_id', 'long_read', 'longread_type']
                
                missing_columns = [col for col in required_columns if col not in header_map]
                if missing_columns:
                    raise ValueError(f"Missing required columns in sample list file for {args.module} mode: {', '.join(missing_columns)}")


                # Read each line of sample information
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    parts = line.strip().split('\t')
                    if len(parts) < len(required_columns):
                        logging.warning(f"Skipping incomplete line: {line.strip()}")
                        continue
                    
                    sample_id = parts[header_map['sample_id']]
                    sample_info = {
                        'read1': os.path.abspath(parts[header_map['read1']]),
                        'read2': os.path.abspath(parts[header_map['read2']])
                    }
                    
                    # Get long read information (if present)
                    if 'long_read' in header_map and len(parts) > header_map['long_read']:
                        long_read = parts[header_map['long_read']]
                        if long_read and long_read.lower() != 'na' and long_read.lower() != 'none':
                            sample_info['long_read'] = os.path.abspath(long_read)
                            
                            # Get long read type (if present)
                            if 'longread_type' in header_map and len(parts) > header_map['longread_type']:
                                longread_type = parts[header_map['longread_type']]
                                if longread_type and longread_type.lower() not in ('na', 'none'):
                                    sample_info['longread_type'] = longread_type
                    
                    # Get genome size (if present)
                    if 'genome_size' in header_map and len(parts) > header_map['genome_size']:
                        try:
                            genome_size = int(parts[header_map['genome_size']])
                            sample_info['genome_size'] = genome_size
                            logging.info(f"Sample {sample_id} genome size: {genome_size} MB")
                        except (ValueError, TypeError):
                            logging.warning(f"Invalid genome size for sample {sample_id}, using default: {config['subsample']['genome_size']} MB")
                    
                    # Get mito_seed (if present for lt-seed mode)
                    if 'mito_seed' in header_map and len(parts) > header_map['mito_seed']:
                        mito_seed = parts[header_map['mito_seed']]
                        if mito_seed and mito_seed.lower() not in ('na', 'none'):
                            sample_info['mito_seed'] = os.path.abspath(mito_seed)
                    
                    sample_ids.append(sample_id)
                    samples_info[sample_id] = sample_info
                    create_sample_directories(args.working_dir, sample_id)
            
            if not sample_ids:
                raise ValueError("No valid samples found in the sample list file")
                
            # Add sample ID list and information to configuration
            config["sample_id"] = sample_ids
            config["samples_info"] = samples_info
            
            # Add sample list file path to configuration
            config["sample_list"] = os.path.abspath(args.sample_list)
        
        # Handle unlock and clean options
        if args.unlock:
            logging.info("Unlocking working directory...")
            cmd = [
                "snakemake",
                "--snakefile", os.path.join(os.path.dirname(__file__), "snakefile"),
                "--unlock",
                "--directory", args.working_dir
            ]
            try:
                subprocess.run(cmd, check=True)
                logging.info("Successfully unlocked working directory")
                sys.exit(0)
            except subprocess.CalledProcessError as e:
                logging.error(f"Failed to unlock working directory: {e}")
                sys.exit(1)

        if args.clean:
            logging.info("Cleaning .snakemake directory...")
            clean_snakemake_dir(args.working_dir)
            logging.info("Cleaned .snakemake directory")
            sys.exit(0)
        
        # Add required parameters to config
        config.update({
            "working_dir": os.path.abspath(args.working_dir),
            "module_to_run": args.module,
            "script_dir": os.path.dirname(os.path.abspath(__file__)),
            "dag_mode": args.dry_run
        })
        
        # Add mito_seed parameter to configuration
        if args.mito_seed:
            config["mito_seed"] = os.path.abspath(args.mito_seed)
            logging.info(f"Using user-provided mitochondrial seed: {config['mito_seed']}")
        
        # Add sample information
        if args.sample_id:
            config.update({
                "sample_id": [args.sample_id],  # Use list format
            })
            
            # Add different parameters based on the mode
            if args.module in ["st", "ht"]:
                # Short-read and hybrid modes require read1 and read2
                config.update({
                    "read1": os.path.abspath(args.read1),
                    "read2": os.path.abspath(args.read2)
                })
            
            if args.long_read:
                # Hybrid and long-read modes require long_read
                config.update({
                    "long_read": os.path.abspath(args.long_read),
                    "longread_type": args.longread_type
                })
            
            if args.module == "lt-seed" and args.mito_seed:
                # Seed-guided mode requires mito_seed
                config.update({
                    "mito_seed": os.path.abspath(args.mito_seed)
                })
        
        # Create configuration file
        config_file = os.path.join(args.working_dir, "config.yaml")
        with open(config_file, 'w') as f:
            yaml.dump(config, f)
        logging.info(f"Configuration saved to {config_file}")
        
        # Prepare Snakemake command
        cmd = [
            "snakemake",
            "--snakefile", os.path.join(os.path.dirname(__file__), "snakefile"),
            "--configfile", config_file,
            "--cores", str(args.cores),
            "--jobs", str(args.jobs),
            "--use-conda",
            "--conda-frontend", "mamba",
            "--conda-prefix", envs_dir
        ]
        
        if args.conda_path:
            cmd.extend([
                "--conda-base-path", conda_base
            ])
        
        # Add optional arguments
        if args.profile:
            cmd.extend(["--profile", args.profile])
        if args.dry_run:
            cmd.append("--dryrun")
        if args.verbose:
            cmd.append("--verbose")
        
        # Add keep-going argument
        if args.keep_going:
            cmd.extend(["--keep-going"])
        
        # Log the command
        logging.info(f"Running command: {' '.join(cmd)}")
        
        # Run Snakemake
        try:
            logging.info("Starting Snakemake workflow...")
            result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            logging.info("Snakemake workflow completed successfully")
            # Output Snakemake output
            print(result.stdout)
            if result.stderr:
                logging.warning("Snakemake stderr output:")
                print(result.stderr)
            
            # Generate error summary if workflow failed
            if "error" in result.stderr.lower():
                logging.error("Errors detected in workflow. Generating error summary...")
                try:
                    summary_cmd = ["grep", "-A", "3", "-B", "3", "ERROR", "**/logs/*.log"]
                    summary = subprocess.run(summary_cmd, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                    print("\n==== ERROR SUMMARY ====")
                    print(summary.stdout)
                    print("=======================\n")
                except Exception as e:
                    logging.error(f"Failed to generate error summary: {e}")
        except subprocess.CalledProcessError as e:
            logging.error(f"Snakemake failed with error code {e.returncode}")
            if e.stdout:
                logging.info("Snakemake stdout output:")
                print(e.stdout)
            if e.stderr:
                logging.error("Snakemake stderr output:")
                print(e.stderr)
            sys.exit(1)
    
    except Exception as e:
        logging.error(f"Pipeline execution failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
