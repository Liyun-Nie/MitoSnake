#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import glob
import sys

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Generate resource usage summary for Snakemake workflow")
    parser.add_argument("--sample", required=True, help="Sample ID")
    parser.add_argument("--output", required=True, help="Output file path")
    parser.add_argument("--logs-dir", help="Logs directory path, default is {sample}/logs")
    parser.add_argument("--module", choices=["st", "ht"], help="Run module, auto-detect by default")
    parser.add_argument("--verbose", action="store_true", help="Display detailed output")
    
    return parser.parse_args()

def get_benchmark_files(logs_dir, module=None):
    """Get list of benchmark files"""
    # Get all .benchmark.txt files
    benchmark_files = glob.glob(os.path.join(logs_dir, "*.benchmark.txt"))
    
    # If module is not specified, try to auto-detect
    if module is None:
        # Check if flye.benchmark.txt or polish.benchmark.txt exists
        hybrid_files = [f for f in benchmark_files if os.path.basename(f) in ["flye.benchmark.txt", "polish.benchmark.txt"]]
        module = "ht" if hybrid_files else "st"
    
    # Filter files based on module
    if module == "st":
        # Rules included in short-read mode
        rules = ["fastp", "seqkit", "spades"]
    else:
        # Rules included in hybrid mode
        rules = ["fastp", "seqkit", "spades", "tiara", "minimap2", "filter_longreads", "flye", "polish"]
    
    # Filter files
    filtered_files = [f for f in benchmark_files if os.path.basename(f).replace(".benchmark.txt", "") in rules]
    
    return filtered_files

def generate_summary(benchmark_files, output_file, verbose=False):
    """Generate resource usage summary"""
    if verbose:
        print(f"Processing {len(benchmark_files)} benchmark files...")
    
    dfs = []
    for bench_file in benchmark_files:
        if os.path.exists(bench_file):
            rule_name = os.path.basename(bench_file).replace(".benchmark.txt", "")
            try:
                df = pd.read_csv(bench_file, sep="\t")
                df['rule'] = rule_name
                dfs.append(df)
                if verbose:
                    print(f"Read {bench_file}")
            except Exception as e:
                print(f"Warning: Cannot read {bench_file}: {e}")
    
    if not dfs:
        print("Warning: No benchmark files found")
        # Create an empty summary file
        with open(output_file, 'w') as f:
            f.write("# No benchmark files found\n")
        return
    
    # Merge all dataframes
    summary_df = pd.concat(dfs)
    
    # Group by rule and calculate statistics
    summary = summary_df.groupby('rule').agg({
        's': ['mean', 'max'],
        'max_rss': ['mean', 'max'],
        'max_vms': ['mean', 'max'],
        'io_in': ['sum'],
        'io_out': ['sum'],
        'mean_load': ['mean', 'max']
    }).round(2)
    
    # Add readable time and memory units
    summary.columns = ['_'.join(col) for col in summary.columns]
    summary['runtime_hours'] = summary['s_mean'] / 3600
    summary['max_memory_gb'] = summary['max_rss_max'] / 1024 / 1024
    
    # Save summary
    summary.to_csv(output_file, sep="\t")
    
    if verbose:
        print(f"Summary saved to {output_file}")
        print("\nResource usage summary:")
        print(summary)

def main():
    """Main function"""
    args = parse_arguments()
    
    # Set logs directory
    logs_dir = args.logs_dir if args.logs_dir else os.path.join(args.sample, "logs")
    
    # Ensure logs directory exists
    if not os.path.exists(logs_dir):
        print(f"Error: Logs directory {logs_dir} does not exist")
        sys.exit(1)
    
    # Get benchmark files
    benchmark_files = get_benchmark_files(logs_dir, args.module)
    
    # Generate summary
    generate_summary(benchmark_files, args.output, args.verbose)

if __name__ == "__main__":
    main() 