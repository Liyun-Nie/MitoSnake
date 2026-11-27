# snakefile

import os
import sys
import re
import pandas as pd
from snakemake.shell import shell
from tabulate import tabulate
import logging
import datetime

# Get command line parameters from Snakemake config
working_dir = config.get("working_dir", os.getcwd())
module_to_run = config.get("module_to_run", "ht")
script_dir = config.get("script_dir", os.path.dirname(workflow.snakefile))
dag_mode = config.get("dag_mode", False)
resource_mode = config.get("resource_mode", "default")
generate_resource_summary = config.get("generate_resource_summary", False)

# Define global variables for all rules
globalvars = {
    "WORKING_DIR": working_dir,
    "SCRIPT_DIR": script_dir,
    "TMP_DIR": "tmp"
}

# Add global variables to workflow namespace
globals().update(globalvars)

# Function to get resources based on rule name
def get_resources(rule_name=None):
    """Get resource requirements for a rule.
    
    Args:
        rule_name: Name of the rule. If None, uses default resources.
    
    Returns:
        Dictionary containing resource specifications.
    """
    # Get global resource mode
    global_mode = config.get("resource_mode", "default")
    
    # Get global resource configuration
    global_resources = config.get("Resources", {}).get(global_mode, {})
    
    # Check for tool-specific resource configuration
    tool_resources = {}
    if rule_name and rule_name in config.get("tool_resources", {}):
        tool_resources = config.get("tool_resources", {}).get(rule_name, {})
    
    # Merge resource configurations, tool-specific configurations take precedence
    return {
        "threads": tool_resources.get("threads", global_resources.get("threads", 4)),
        "mem_mb": tool_resources.get("mem_mb", global_resources.get("mem_mb", 8192)),
        "disk_mb": tool_resources.get("disk_mb", global_resources.get("disk_mb", 10240)),
        "runtime": tool_resources.get("runtime", global_resources.get("runtime", 7200))
    }

# Validate module choice
if module_to_run not in ["st", "ht", "lt-seed"]:
    raise ValueError('Module must be either "st" (short-read task), "ht" (hybrid task), or "lt-seed" (long-read with seed task)')

# Create necessary directories for each sample
for sample in config.get("sample_id", []):
    for dir_path in [
        os.path.join(working_dir, sample, "subsample"),
        os.path.join(working_dir, sample, "assembly"),
        os.path.join(working_dir, sample, "logs")
    ]:
        os.makedirs(dir_path, exist_ok=True)

# Create tmp directory
os.makedirs(os.path.join(working_dir, "tmp"), exist_ok=True)

# Function to print paths
def print_paths(paths):
    print("\n[Paths]")
    path_data = [[name, path] for name, path in paths.items()]
    print(tabulate(path_data, headers=["Path Name", "Path Value"], tablefmt="grid"))

# Function to print samples info
def print_samples_info(samples_info):
    print("\n[Samples Information]")
    sample_info_table = []
    for sample_id, info in samples_info.items():
        sample_info_table.append([
            sample_id,
            info.get('read1', 'N/A'),
            info.get('read2', 'N/A'),
            info.get('long_read', 'N/A'),
            info.get('longread_type', 'N/A'),
            info.get('mito_seed', 'N/A')
        ])
    headers = ["Sample ID", "Read 1", "Read 2", "Long Read", "Long Read Type", "Mito Seed"]
    print(tabulate(sample_info_table, headers=headers, tablefmt="grid"))

# Handle single sample or sample list
samples_info = {}
if config.get("sample_id") and not config.get("sample_list"):
    # Single sample mode (not from sample list)
    sample_ids = config["sample_id"]  # This is a list
    read1 = config.get("read1")
    read2 = config.get("read2")
    long_read = config.get("long_read")
    longread_type = config.get("longread_type")
    mito_seed = config.get("mito_seed")
    
    # Validate required parameters based on module
    if module_to_run in ["st", "ht"]:
        if not (read1 and read2):
            raise ValueError("Read 1 and Read 2 files are required for short-read (st) and hybrid (ht) modes")
        if module_to_run == "ht" and not long_read:
            raise ValueError("Long read file is required for hybrid (ht) mode")
    elif module_to_run == "lt-seed":
        if not (long_read and mito_seed):
            raise ValueError("Long read file and mitochondrial seed file are required for seed-guided (lt-seed) mode")
    
    # Iterate through each sample ID in the list
    for sample_id in sample_ids:
        sample_info = {}
        
        # Add parameters based on module
        if module_to_run in ["st", "ht"]:
            sample_info.update({
                'read1': read1,
                'read2': read2
            })
        
        if long_read:
            sample_info.update({
                'long_read': long_read,
                'longread_type': longread_type
            })
        
        if module_to_run == "lt-seed" and mito_seed:
            sample_info.update({
                'mito_seed': mito_seed
            })
        
        samples_info[sample_id] = sample_info

elif config.get("sample_list"):
    # Multi-sample mode from sample list file
    sample_list_file = config.get("sample_list")
    if not sample_list_file:
        raise ValueError("Either sample_id or sample_list must be provided")
        
    try:
        sample_df = pd.read_csv(sample_list_file, sep='\t')
        
        # Determine required columns based on module
        if module_to_run in ["st", "ht"]:
            required_columns = {'sample_id', 'read1', 'read2'}
            if module_to_run == "ht":
                required_columns.add('long_read')
                required_columns.add('longread_type')
        elif module_to_run == "lt-seed":
            required_columns = {'sample_id', 'long_read', 'longread_type'}
            if 'mito_seed' not in config:
                required_columns.add('mito_seed')
        
        missing_columns = required_columns - set(sample_df.columns)
        if missing_columns:
            raise ValueError(f"Missing required columns in sample list file for {module_to_run} mode: {missing_columns}")
            
        for _, row in sample_df.iterrows():
            sample_id = row['sample_id']
            sample_info = {}
            
            # Add parameters based on the module
            if module_to_run in ["st", "ht"]:
                sample_info.update({
                    'read1': row['read1'] if os.path.isabs(row['read1']) else os.path.abspath(row['read1']),
                    'read2': row['read2'] if os.path.isabs(row['read2']) else os.path.abspath(row['read2'])
                })
            
            # Handle long read if present
            if 'long_read' in row and not pd.isna(row['long_read']):
                long_read = row['long_read'] if os.path.isabs(row['long_read']) else os.path.abspath(row['long_read'])
                sample_info['long_read'] = long_read
                
                # Get long read type if present
                if 'longread_type' in row and not pd.isna(row['longread_type']):
                    sample_info['longread_type'] = row['longread_type']
            
            # Handle mito seed if present
            if 'mito_seed' in row and not pd.isna(row['mito_seed']):
                mito_seed = row['mito_seed'] if os.path.isabs(row['mito_seed']) else os.path.abspath(row['mito_seed'])
                sample_info['mito_seed'] = mito_seed
            
            # Get genome size if present
            if 'genome_size' in row and not pd.isna(row['genome_size']):
                try:
                    genome_size = int(row['genome_size'])
                    sample_info['genome_size'] = genome_size
                except (ValueError, TypeError):
                    pass
            
            samples_info[sample_id] = sample_info
    except Exception as e:
        raise ValueError(f"Error reading sample list file: {e}")

SAMPLES = list(samples_info.keys())

if not SAMPLES:
    raise ValueError("No samples found. Please check your input parameters.")

# Function to get sample list
def get_sample_list():
    """Get the list of samples."""
    return SAMPLES

# Print information in non-DAG mode
if not dag_mode:
    paths = {
        "WORKING_DIR": working_dir,
        "SCRIPT_DIR": script_dir
    }
    print_paths(paths)
    print_samples_info(samples_info)

# Function to get genome size for a sample
def get_genome_size(wildcards):
    """Get the genome size for a sample."""
    if config.get("samples_info") and wildcards.sample in config["samples_info"]:
        if "genome_size" in config["samples_info"][wildcards.sample]:
            return config["samples_info"][wildcards.sample]["genome_size"]
    return config["subsample"]["genome_size"]

# Determine which module to run and set output files
def get_final_outputs():
    """Get final outputs based on module."""
    outputs = []
    
    for sample in get_sample_list():
        outputs.append(f"{sample}/logs/consolidated_workflow.log")
        
        if config["module_to_run"] == "st":
            outputs.extend([
                f"{sample}/assembly/contigs.fasta"
            ])
        elif config["module_to_run"] == "ht":
            outputs.extend([
                f"{sample}/assembly/hybrid_assembly.fasta"
            ])
        elif config["module_to_run"] == "lt-seed":
            iterative_rounds = config.get("iterative_rounds", 2)
            outputs.extend([
                f"{sample}/assembly/final_assembly.fasta",
                f"{sample}/assembly/final_assembly.gfa",
                f"{sample}/assembly/all_rounds_assemblies.fasta",
                f"{sample}/assembly/flye_out_round{iterative_rounds}/assembly.fasta",
                f"{sample}/assembly/flye_out_round{iterative_rounds}/assembly_graph.gfa"
            ])
    
    return outputs

# Include all rules - this must be done after defining get_genome_size function
include: os.path.join(script_dir, "rules/__init__.smk")

# Define the rule all
rule all:
    """
    Main target rule.
    """
    input:
        get_final_outputs()

# Add validate_wildcards function
def validate_wildcards(wildcards):
    """Validate wildcard values."""
    if not isinstance(wildcards.sample, str):
        raise ValueError(f"Invalid sample wildcard: {wildcards.sample}")
    if wildcards.sample not in SAMPLES:
        raise ValueError(f"Unknown sample: {wildcards.sample}")

# Add wildcard_constraints
wildcard_constraints:
    sample = "|".join([re.escape(sample) for sample in SAMPLES])

# Add debugging information
logging.debug(f"Samples: {SAMPLES}")

onstart:
    shell("echo 'Workflow started at {}'".format(datetime.datetime.now()))

onsuccess:
    shell("echo 'Workflow completed successfully at {}'".format(datetime.datetime.now()))

globals()["samples_info"] = samples_info
