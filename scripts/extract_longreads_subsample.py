import subprocess
import sys
import logging
import os
from pathlib import Path

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', stream=sys.stderr)

# Arguments passed from Snakemake
longreads = sys.argv[1]
output = sys.argv[2]
threads = int(sys.argv[3])
whole_genome_size = int(sys.argv[4])
sampling_cov = float(sys.argv[5])

# Fixed minimum length for filtering
MIN_LENGTH = 10000  # 10kb as minimum sequence length threshold

try:
    logging.info(f"Input reads: {longreads}")
    logging.info(f"Output FASTA: {output}")
    logging.info(f"Threads: {threads}")
    logging.info(f"Genome size: {whole_genome_size} MB")
    logging.info(f"Target coverage: {sampling_cov}x")
    logging.info(f"Minimum sequence length filter: {MIN_LENGTH} bp")

    # Determine cat command based on input file extension
    cat_cmd = "zcat" if longreads.endswith(".gz") else "cat"
    logging.info(f"Using '{cat_cmd}' for input file.")

    sum_len = None
    try:
        # Step 1: Calculate stats for sequences after length filtering
        cmd_stats_filtered = f"{cat_cmd} {longreads} | seqtk seq -A | seqkit seq -m {MIN_LENGTH} | seqkit stats -j {threads} -a -T -"
        logging.info(f"Running command to calculate filtered sequences stats: {cmd_stats_filtered}")
        
        process = subprocess.Popen(
            cmd_stats_filtered, 
            shell=True, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE,
            universal_newlines=True
        )
        stats_filtered_output, stderr_filtered_output = process.communicate()
        
        if process.returncode != 0:
            logging.warning(f"Filtered stats command returned non-zero exit status {process.returncode}")
            logging.warning(f"Standard error output: {stderr_filtered_output}")
            raise ValueError(f"Failed to calculate stats for filtered sequences")
            
        logging.info(f"seqkit filtered stats output:\n{stats_filtered_output}")
        
        # Parse filtered stats output
        filtered_lines = stats_filtered_output.strip().split('\n')
        if len(filtered_lines) < 2:
            if len(filtered_lines) > 0 and "file" in filtered_lines[0] and "sum_len" in filtered_lines[0]:
                logging.warning(f"No sequences found with length >= {MIN_LENGTH} bp!")
                sum_len = 0
            else:
                raise ValueError("Unexpected output from seqkit stats - not enough lines")
        else:
            headers = filtered_lines[0].split('\t')
            if 'sum_len' not in headers or 'num_seqs' not in headers:
                raise ValueError(f"Required columns not found in headers: {headers}")

            sum_len_idx = headers.index('sum_len')
            num_seqs_idx = headers.index('num_seqs')
            sum_len = int(filtered_lines[1].split('\t')[sum_len_idx])
            num_seqs = int(filtered_lines[1].split('\t')[num_seqs_idx])
            
            logging.info(f"Filtered sequences (>= {MIN_LENGTH} bp): {num_seqs} sequences, {sum_len} total bases")
            
            if num_seqs == 0:
                logging.warning(f"No sequences found with length >= {MIN_LENGTH} bp! Output will be empty.")
                sum_len = 0

    except Exception as e:
        # Catch any errors during stats calculation
        logging.error(f"Error calculating filtered sequence stats: {e}")
        # We won't continue if we can't get proper stats
        raise

    # Ensure sum_len is not zero to avoid division by zero error
    if sum_len == 0:
        logging.warning("Total filtered length calculated as 0. Cannot calculate sampling proportion. Output will be empty.")
        p = 0.0
        actual_depth = 0.0
    else:
        # Calculate the proportion (p) to extract and the actual depth
        p = (whole_genome_size * sampling_cov * 1000000) / sum_len
        actual_depth = sampling_cov if p < 1.0 else sum_len / (whole_genome_size * 1000000)

    logging.info(f"Calculated sampling proportion: {p:.6f} ({p*100:.2f}%)")
    logging.info(f"Actual extraction depth: {actual_depth:.2f}x")

    # Adjust output filename if actual depth is different from target
    output_base = os.path.dirname(output)
    original_filename = os.path.basename(output)
    
    # Update output filename if actual depth differs significantly from target
    if abs(actual_depth - sampling_cov) > 0.1 and p < 1.0:
        # Extract base name pattern and create new name with actual depth
        prefix = original_filename.split('_subsample_')[0]
        new_filename = f"{prefix}_subsample_{actual_depth:.1f}x.fa"
        new_output = os.path.join(output_base, new_filename)
        logging.info(f"Actual depth differs from target. Updating output filename to: {new_filename}")
        output = new_output

    # Ensure p is between 0.0 and 1.0
    if p < 0.0:
        p = 0.0
    elif p > 1.0:
        logging.warning(f"Sampling proportion {p:.6f} exceeds 1.0, using 1.0 (100%)")
        p = 1.0

    # Determine if full copy/conversion or sampling is required
    if p >= 1.0:
        # Convert all filtered sequences to FASTA
        cmd_copy = f"{cat_cmd} {longreads} | seqtk seq -A | seqkit seq -m {MIN_LENGTH} > {output}"
        logging.info(f"Proportion >= 1.0. Converting full input to FASTA with length filtering >= {MIN_LENGTH}bp.")
        logging.info(f"Running command: {cmd_copy}")
        
        process = subprocess.run(
            cmd_copy, 
            shell=True, 
            capture_output=True,
            text=True,
            check=False
        )
        
        if process.returncode != 0:
            logging.error(f"Command failed with exit code {process.returncode}")
            if process.stderr:
                logging.error(f"Error output: {process.stderr}")
            raise subprocess.CalledProcessError(process.returncode, cmd_copy, 
                                               output=process.stdout, stderr=process.stderr)
        
    elif p > 0.0:
        # Sample from filtered sequences
        cmd_sample = f"{cat_cmd} {longreads} | seqtk seq -A | seqkit seq -m {MIN_LENGTH} | seqtk sample -s 100 - {p} > {output}"
        logging.info(f"Proportion < 1.0. Sampling with length filtering >= {MIN_LENGTH}bp.")
        logging.info(f"Running command: {cmd_sample}")
        
        process = subprocess.run(
            cmd_sample, 
            shell=True, 
            capture_output=True,
            text=True,
            check=False
        )
        
        if process.returncode != 0:
            logging.error(f"Command failed with exit code {process.returncode}")
            if process.stderr:
                logging.error(f"Error output: {process.stderr}")
            raise subprocess.CalledProcessError(process.returncode, cmd_sample, 
                                               output=process.stdout, stderr=process.stderr)
        
    else: # p == 0.0
        logging.info("Sampling proportion is 0. Creating an empty output file.")
        Path(output).touch() # Create empty file

    # Verify output file exists
    if not os.path.exists(output):
        logging.error("Output file was not created.")
        Path(output).touch() # Ensure file exists anyway
        sys.exit(1)
    elif os.path.getsize(output) == 0 and p > 0.0:
        logging.warning("Output file is empty despite sampling proportion > 0. Input might have been effectively empty after filtering.")
    else:
        logging.info(f"Subsampling/conversion complete. Output saved to {output}")

except subprocess.CalledProcessError as e:
    logging.error(f"Error occurred during subprocess execution: {str(e)}")
    logging.error(f"Command: {e.cmd}")
    if hasattr(e, 'stderr') and e.stderr:
        logging.error(f"Error output: {e.stderr}")
    Path(output).touch()
    sys.exit(1)
except Exception as e:
    logging.error(f"General error occurred: {str(e)}")
    Path(output).touch()
    sys.exit(1)