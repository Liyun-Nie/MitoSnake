#!/usr/bin/env python3

import argparse
import os
import re
import sys
import glob
from collections import defaultdict
from tabulate import tabulate

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Generate error report for MitoSnake workflow")
    parser.add_argument("--log-file", help="Snakemake log file path")
    parser.add_argument("--working-dir", required=True, help="Working directory")
    parser.add_argument("--output", default="error_report.txt", help="Output file path")
    parser.add_argument("--verbose", action="store_true", help="Display detailed output")
    
    return parser.parse_args()

def parse_snakemake_log(log_file):
    """Parse Snakemake log file to extract error information"""
    if not os.path.exists(log_file):
        return None
    
    with open(log_file, 'r') as f:
        content = f.read()
    
    # Extract error information
    error_info = defaultdict(list)
    
    # Pattern to extract error in rule
    error_pattern = re.compile(r"Error in rule ([^:]+):\s+jobid: (\d+)\s+input: ([^\n]+)\s+output: ([^\n]+)\s+log: ([^\n]+)")
    
    for match in error_pattern.finditer(content):
        rule_name, job_id, inputs, outputs, logs = match.groups()
        error_info[rule_name].append({
            'job_id': job_id,
            'inputs': inputs.strip(),
            'outputs': outputs.strip(),
            'logs': logs.strip()
        })
    
    return error_info

def collect_sample_errors(working_dir):
    """Collect error information from individual sample logs"""
    sample_errors = defaultdict(list)
    
    # Find all sample directories
    sample_dirs = [d for d in os.listdir(working_dir) 
                  if os.path.isdir(os.path.join(working_dir, d)) and not d.startswith('.')]
    
    for sample in sample_dirs:
        logs_dir = os.path.join(working_dir, sample, 'logs')
        if not os.path.exists(logs_dir):
            continue
        
        # Check all log files for errors
        log_files = glob.glob(os.path.join(logs_dir, '*.log'))
        for log_file in log_files:
            try:
                rule_name = os.path.basename(log_file).replace('.log', '')
                with open(log_file, 'r') as f:
                    content = f.read()
                
                # Look for error messages
                if 'error' in content.lower() or 'exception' in content.lower():
                    # Find the last few lines that contain error information
                    lines = content.splitlines()
                    error_lines = []
                    for line in reversed(lines):
                        if 'error' in line.lower() or 'exception' in line.lower() or 'warning' in line.lower():
                            error_lines.append(line)
                        if len(error_lines) >= 3:
                            break
                    
                    sample_errors[sample].append({
                        'rule': rule_name,
                        'log_file': log_file,
                        'error_details': '\n'.join(reversed(error_lines))
                    })
            except Exception as e:
                print(f"Warning: Error parsing log file {log_file}: {e}")
    
    return sample_errors

def generate_report(error_info, sample_errors, output_file, verbose=False):
    """Generate error report"""
    with open(output_file, 'w') as f:
        f.write("===== MITOSNAKE ERROR REPORT =====\n\n")
        
        # Write workflow error summary
        f.write("== Workflow Error Summary ==\n")
        if error_info:
            error_table = []
            for rule_name, errors in error_info.items():
                for error in errors:
                    # Extract sample name from input/output paths
                    sample_match = re.search(r'([^/]+)/assembly|([^/]+)/subsample', error['inputs'])
                    sample_name = sample_match.group(1) if sample_match and sample_match.group(1) else \
                                 (sample_match.group(2) if sample_match and sample_match.group(2) else "Unknown sample")
                    
                    error_table.append([
                        sample_name,
                        rule_name,
                        error['job_id'],
                        os.path.basename(error['logs'].split(' ')[0])
                    ])
            
            if error_table:
                f.write(tabulate(error_table, headers=["Sample", "Rule", "Job ID", "Log File"], tablefmt="grid"))
                f.write("\n\n")
            else:
                f.write("No workflow errors recorded\n\n")
        else:
            f.write("No Snakemake log file found\n\n")
        
        # Write sample-specific error details
        f.write("== Sample Error Details ==\n")
        if sample_errors:
            for sample, errors in sample_errors.items():
                f.write(f"\nSample: {sample}\n")
                f.write("-" * 40 + "\n")
                
                for error in errors:
                    f.write(f"Rule: {error['rule']}\n")
                    f.write(f"Log file: {error['log_file']}\n")
                    f.write("Error message:\n")
                    f.write(error['error_details'])
                    f.write("\n" + "-" * 40 + "\n")
        else:
            f.write("No sample error logs found\n")
        
        # Add recommendation
        f.write("\n== Troubleshooting Recommendations ==\n")
        f.write("1. Check the error log files mentioned above for detailed information\n")
        f.write("2. Common causes of errors:\n")
        f.write("   - Input files are empty or missing\n")
        f.write("   - Insufficient memory or disk space\n")
        f.write("   - Missing required software or dependencies\n")
        f.write("3. If tiara_longreads rule failed, it might be due to no mitochondrial sequences found\n")
        f.write("4. If minimap2 rule failed, check if mitochondrion_contigs.fasta is valid\n")
        
        f.write("\n===== REPORT GENERATION COMPLETE =====\n")
    
    if verbose:
        with open(output_file, 'r') as f:
            print(f.read())
    
    print(f"Error report saved to: {output_file}")

def find_latest_snakemake_log(working_dir):
    """Find the latest Snakemake log file"""
    log_dir = os.path.join(working_dir, ".snakemake", "log")
    if not os.path.exists(log_dir):
        return None
    
    log_files = glob.glob(os.path.join(log_dir, "*.snakemake.log"))
    if not log_files:
        return None
    
    # Sort by modification time (newest first)
    latest_log = max(log_files, key=os.path.getmtime)
    return latest_log

def main():
    """Main function"""
    args = parse_arguments()
    
    # Find the latest Snakemake log file if not specified
    log_file = args.log_file
    if not log_file:
        log_file = find_latest_snakemake_log(args.working_dir)
        if log_file:
            print(f"Using latest Snakemake log file: {log_file}")
        else:
            print("Warning: No Snakemake log file found")
    
    # Parse Snakemake log
    error_info = parse_snakemake_log(log_file) if log_file else None
    
    # Collect errors from sample logs
    sample_errors = collect_sample_errors(args.working_dir)
    
    # Generate report
    generate_report(error_info, sample_errors, args.output, args.verbose)

if __name__ == "__main__":
    main()