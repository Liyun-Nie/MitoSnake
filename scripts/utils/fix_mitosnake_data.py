#!/usr/bin/env python3
# fix_mitosnake_data.py - Fix problematic FASTQ data for mitosnake workflow

import os
import sys
import argparse
import subprocess
import logging
import multiprocessing
import shutil

def setup_logger(log_file=None):
    """Set up logging"""
    logger = logging.getLogger('fix_mitosnake_data')
    logger.setLevel(logging.INFO)
    
    # Console handler
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(console)
    
    # File handler if specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)
    
    return logger

def run_command(cmd, logger):
    """Run command and log output"""
    try:
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return True, result.stdout
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {' '.join(cmd)}")
        logger.error(f"Error: {e.stderr}")
        return False, e.stderr

def fix_fastq_data(input_file, output_file, input_format="fastq", threads=1, logger=None):
    """Fix sequence file using seqkit sana"""
    cmd = ['seqkit', 'sana', '-i', input_format, '-j', str(threads), input_file, '-o', output_file]
    logger.info(f"Fixing file with sana: {input_file} with {threads} threads")
    success, output = run_command(cmd, logger)
    
    if not success:
        logger.error(f"Failed to fix {input_file}")
        return False
    
    if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
        logger.error(f"Output file is empty or missing: {output_file}")
        return False
    
    return True

def pair_fastq_files(input_file1, input_file2, output_dir, threads=1, logger=None):
    """Pair FASTQ files using seqkit pair"""
    # Create a temporary pairing directory to store paired output
    pair_output_dir = os.path.join(output_dir, "pair_temp")
    os.makedirs(pair_output_dir, exist_ok=True)
    
    # seqkit pair expects read1, read2 and out-dir
    cmd = ['seqkit', 'pair', 
           '--read1', input_file1, 
           '--read2', input_file2, 
           '--out-dir', pair_output_dir, 
           '-j', str(threads)]
    
    logger.info(f"Pairing files: {input_file1} and {input_file2} with {threads} threads")
    success, output = run_command(cmd, logger)
    
    if not success:
        logger.error(f"Failed to pair {input_file1} and {input_file2}")
        return False
    
    # Check if paired files were created
    paired_read1 = os.path.join(pair_output_dir, os.path.basename(input_file1))
    paired_read2 = os.path.join(pair_output_dir, os.path.basename(input_file2))
    
    if not os.path.exists(paired_read1) or os.path.getsize(paired_read1) == 0 or \
       not os.path.exists(paired_read2) or os.path.getsize(paired_read2) == 0:
        logger.error(f"Paired output files are empty or missing: {paired_read1}, {paired_read2}")
        return False
    
    # Move files from temporary directory to final location
    # The output directory should already exist, as it was created for sana output
    final_read1 = os.path.join(output_dir, os.path.basename(input_file1).replace("sana_", ""))
    final_read2 = os.path.join(output_dir, os.path.basename(input_file2).replace("sana_", ""))
    
    try:
        shutil.move(paired_read1, final_read1)
        shutil.move(paired_read2, final_read2)
        # Clean up temporary directory
        shutil.rmtree(pair_output_dir)
    except Exception as e:
        logger.error(f"Failed to move paired files: {e}")
        return False
    
    if not os.path.exists(final_read1) or os.path.getsize(final_read1) == 0 or \
       not os.path.exists(final_read2) or os.path.getsize(final_read2) == 0:
        logger.error(f"Final paired files are empty or missing: {final_read1}, {final_read2}")
        return False
    
    return True

def parse_sample_list(sample_list_file, logger):
    """Parse mitosnake sample list to extract file paths"""
    sample_data = {}
    
    logger.info(f"Reading sample list: {sample_list_file}")
    try:
        with open(sample_list_file, 'r') as f:
            # Skip header
            header = f.readline().strip()
            
            # Parse each line
            for line in f:
                if line.strip() and not line.startswith('#'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:  # At least has sample_id, read1, read2
                        sample_id = parts[0]
                        sample_data[sample_id] = {
                            'read1': parts[1] if parts[1] != "NA" else None,
                            'read2': parts[2] if parts[2] != "NA" else None
                        }
                        
                        # Add long read data if available
                        if len(parts) >= 4:
                            sample_data[sample_id]['long_read'] = parts[3] if parts[3] != "NA" else None
                            
                            if len(parts) >= 5:
                                sample_data[sample_id]['longread_type'] = parts[4]
                                
                                if len(parts) >= 6:
                                    sample_data[sample_id]['genome_size'] = parts[5]
    except Exception as e:
        logger.error(f"Error parsing sample list: {e}")
        return {}
    
    logger.info(f"Found {len(sample_data)} samples in the list")
    return sample_data

def create_fixed_sample_list(original_list, sample_replacements, output_file, logger):
    """Create a new sample list with fixed data paths"""
    try:
        with open(original_list, 'r') as f:
            content = f.readlines()
        
        # Process each line
        new_content = []
        for line in content:
            new_line = line
            
            # Apply replacements
            for sample_id, replacements in sample_replacements.items():
                for old_path, new_path in replacements.items():
                    new_line = new_line.replace(old_path, new_path)
            
            new_content.append(new_line)
        
        # Write new file
        with open(output_file, 'w') as f:
            f.writelines(new_content)
        
        logger.info(f"Created fixed sample list: {output_file}")
        return True
    except Exception as e:
        logger.error(f"Error creating fixed sample list: {e}")
        return False

def process_sample(sample_id, sample_info, output_dir, threads, logger, keep_intermediate=False):
    """Process single sample data"""
    replacements = {}
    
    # Process short reads
    if sample_info.get('read1') and sample_info.get('read2'):
        # Paired-end data needs both sana and pair processing
        read1 = sample_info['read1']
        read2 = sample_info['read2']
        
        if os.path.exists(read1) and os.path.exists(read2):
            # Create output directories in original file locations
            read1_dir = os.path.dirname(read1)
            read2_dir = os.path.dirname(read2)
            
            out_dir1 = os.path.join(read1_dir, output_dir)
            out_dir2 = os.path.join(read2_dir, output_dir)
            
            os.makedirs(out_dir1, exist_ok=True)
            if out_dir1 != out_dir2:
                os.makedirs(out_dir2, exist_ok=True)
            
            # Create intermediate file paths directly in the output directory
            sana_out_read1 = os.path.join(out_dir1, "sana_" + os.path.basename(read1))
            sana_out_read2 = os.path.join(out_dir2, "sana_" + os.path.basename(read2))
            
            # Output file paths will be created by pair_fastq_files
            out_read1 = os.path.join(out_dir1, os.path.basename(read1))
            out_read2 = os.path.join(out_dir2, os.path.basename(read2))
            
            # Step 1: Fix with sana
            r1_fixed = fix_fastq_data(read1, sana_out_read1, "fastq", threads, logger)
            r2_fixed = fix_fastq_data(read2, sana_out_read2, "fastq", threads, logger)
            
            if r1_fixed and r2_fixed:
                # Step 2: Apply pairing using the output directory
                paired_success = pair_fastq_files(
                    sana_out_read1, sana_out_read2,
                    out_dir1,  # Use the first output directory for pairing
                    threads, logger
                )
                
                # Clean up intermediate files if not keeping them
                if not keep_intermediate:
                    try:
                        if os.path.exists(sana_out_read1):
                            os.remove(sana_out_read1)
                        if os.path.exists(sana_out_read2):
                            os.remove(sana_out_read2)
                    except Exception as e:
                        logger.warning(f"Failed to clean up intermediate files: {e}")
                
                if paired_success:
                    logger.info(f"Successfully fixed paired-end reads for {sample_id}")
                    replacements[read1] = out_read1
                    replacements[read2] = out_read2
                else:
                    logger.error(f"Failed to pair reads for {sample_id}")
            else:
                logger.error(f"Failed to fix short reads for {sample_id}")
                # Clean up any created intermediate files
                if not keep_intermediate:
                    try:
                        if os.path.exists(sana_out_read1):
                            os.remove(sana_out_read1)
                        if os.path.exists(sana_out_read2):
                            os.remove(sana_out_read2)
                    except:
                        pass
        else:
            logger.warning(f"Short read files not found for {sample_id}")
    elif sample_info.get('read1') and not sample_info.get('read2'):
        # Single-end data only needs sana processing
        read1 = sample_info['read1']
        
        if os.path.exists(read1):
            # Create output directory
            read1_dir = os.path.dirname(read1)
            out_dir1 = os.path.join(read1_dir, output_dir)
            os.makedirs(out_dir1, exist_ok=True)
            
            # Create output file path
            out_read1 = os.path.join(out_dir1, os.path.basename(read1))
            
            # Only use sana to fix
            r1_fixed = fix_fastq_data(read1, out_read1, "fastq", threads, logger)
            
            if r1_fixed:
                logger.info(f"Successfully fixed single-end read for {sample_id}")
                replacements[read1] = out_read1
            else:
                logger.error(f"Failed to fix single-end read for {sample_id}")
        else:
            logger.warning(f"Single-end read file not found for {sample_id}")
    
    # Process long reads (if present) - only needs sana processing
    if sample_info.get('long_read'):
        long_read = sample_info['long_read']
        
        if os.path.exists(long_read):
            # Create output directory
            long_read_dir = os.path.dirname(long_read)
            out_dir_long = os.path.join(long_read_dir, output_dir)
            os.makedirs(out_dir_long, exist_ok=True)
            
            # Create output file path
            out_long_read = os.path.join(out_dir_long, os.path.basename(long_read))
            
            # Long reads only use sana
            long_fixed = fix_fastq_data(long_read, out_long_read, "fastq", threads, logger)
            
            if long_fixed:
                logger.info(f"Successfully fixed long read for {sample_id}")
                replacements[long_read] = out_long_read
            else:
                logger.error(f"Failed to fix long read for {sample_id}")
        else:
            logger.warning(f"Long read file not found for {sample_id}")
    
    return replacements

def main():
    parser = argparse.ArgumentParser(description="Fix problematic FASTQ data for mitosnake workflow")
    
    # Required arguments
    parser.add_argument("--sample-list", required=True, help="Original mitosnake sample list file")
    
    # Optional arguments
    parser.add_argument("--output-dir", default="fixed_data", 
                       help="Output directory name for fixed files (will be created under each original data directory)")
    parser.add_argument("--output-list", help="Output sample list file (default: sample_list_fixed.txt)")
    parser.add_argument("--log-file", help="Log file path")
    parser.add_argument("--threads", type=int, default=multiprocessing.cpu_count(), 
                       help=f"Number of threads to use (default: {multiprocessing.cpu_count()})")
    parser.add_argument("--keep-intermediate", action="store_true",
                       help="Keep intermediate files after processing")
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logger(args.log_file)
    
    # Log thread count
    logger.info(f"Using {args.threads} threads for processing")
    if args.keep_intermediate:
        logger.info("Intermediate files will be kept after processing")
    
    # Parse sample list
    sample_data = parse_sample_list(args.sample_list, logger)
    if not sample_data:
        logger.error("Failed to parse sample list or no samples found")
        sys.exit(1)
    
    # Track replacements for updating sample list
    all_replacements = {}
    
    # Process each sample
    for sample_id, sample_info in sample_data.items():
        logger.info(f"Processing sample: {sample_id}")
        
        replacements = process_sample(
            sample_id, 
            sample_info, 
            args.output_dir, 
            args.threads, 
            logger,
            args.keep_intermediate
        )
        
        if replacements:
            all_replacements[sample_id] = replacements
    
    # Create fixed sample list
    if all_replacements:
        output_list = args.output_list or f"{os.path.splitext(args.sample_list)[0]}_fixed.txt"
        if create_fixed_sample_list(args.sample_list, all_replacements, output_list, logger):
            logger.info(f"Successfully created fixed sample list: {output_list}")
            logger.info(f"Run mitosnake with: --sample-list {output_list}")
        else:
            logger.error("Failed to create fixed sample list")
    else:
        logger.warning("No replacements made, sample list unchanged")

if __name__ == "__main__":
    main()