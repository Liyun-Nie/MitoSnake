import re
import subprocess
import sys
import logging
import os
from pathlib import Path

# Setup logging to file
def setup_logging(log_file):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    
    # File handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(file_handler)
    
    # Also write to stderr
    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(stderr_handler)
    
    return logger

# Arguments from Snakemake
fq1 = sys.argv[1]
fq2 = sys.argv[2]
output1 = sys.argv[3]
output2 = sys.argv[4]
threads = int(sys.argv[5])
whole_genome_size = int(sys.argv[6])
sampling_cov = float(sys.argv[7])
log_stats = sys.argv[8]
log_subsample = sys.argv[9]

# Setup logging
logger = setup_logging(log_subsample)

try:
    logger.info(f"Input R1: {fq1}")
    logger.info(f"Input R2: {fq2}")
    logger.info(f"Output R1: {output1}")
    logger.info(f"Output R2: {output2}")
    logger.info(f"Threads: {threads}")
    logger.info(f"Genome size: {whole_genome_size} MB")
    logger.info(f"Target coverage: {sampling_cov}x")
    
    # Try using seqkit stats first
    try:
        # Step 1: Run seqkit stats to get the total length (sum_len) of the input files
        cmd_stats = f"seqkit stats -j {threads} {fq1} {fq2} -a -T"
        logger.info(f"Running command: {cmd_stats}")
        stats_output = subprocess.check_output(cmd_stats, shell=True, universal_newlines=True)
        
        # Write stats output to separate log
        with open(log_stats, 'w') as log_file:
            log_file.write(stats_output)
        
        logger.info(f"seqkit stats output:\n{stats_output}")
        
        # Step 2: Extract the total length (sum_len) for both R1 and R2 sequences
        lines = stats_output.strip().split('\n')
        if len(lines) < 3:
            raise ValueError("Unexpected output from seqkit stats - not enough lines")
        
        # Get headers and find index of sum_len
        headers = lines[0].split('\t')
        if 'sum_len' not in headers:
            raise ValueError(f"sum_len column not found in headers: {headers}")
        
        sum_len_idx = headers.index('sum_len')
        sum_len_r1 = int(lines[1].split('\t')[sum_len_idx])
        sum_len_r2 = int(lines[2].split('\t')[sum_len_idx])
        
    except Exception as e:
        # Fallback to seqtk fqchk for each file separately
        logger.warning(f"seqkit stats failed: {e}")
        logger.info("Falling back to seqtk fqchk...")
        
        # Function to get bases from seqtk fqchk
        def get_bases_from_fqchk(fastq_file):
            # Determine if input is gzipped
            cat_cmd = "zcat" if fastq_file.endswith(".gz") else "cat"
            
            cmd_fqchk = f"{cat_cmd} {fastq_file} | seqtk fqchk -"
            logger.info(f"Running command: {cmd_fqchk}")
            fqchk_output = subprocess.check_output(cmd_fqchk, shell=True, universal_newlines=True)
            
            # Example output: ALL     724087  13.7 28.8    0.0 0.5 0.3 0.6 0.0 0.2    16  38  23 4.5e+09    0    0
            lines = fqchk_output.strip().split('\n')
            first_line = lines[0].split()
            if len(first_line) < 14:
                raise ValueError(f"Unexpected format from seqtk fqchk: {first_line}")
            
            # The 14th column typically contains the sum of all bases
            bases = float(first_line[13])
            if 'e+' in first_line[13]:  # Handle scientific notation
                bases = int(float(first_line[13]))
            else:
                bases = int(first_line[13])
            
            return bases
        
        # Get bases from both files
        sum_len_r1 = get_bases_from_fqchk(fq1)
        sum_len_r2 = get_bases_from_fqchk(fq2)
        
        # Write stats to the log files
        with open(log_stats, 'w') as log_file:
            log_file.write(f"file\tformat\ttype\tnum_seqs\tsum_len\n")
            log_file.write(f"{fq1}\tfastq\tDNA\tN/A\t{sum_len_r1}\n")
            log_file.write(f"{fq2}\tfastq\tDNA\tN/A\t{sum_len_r2}\n")
    
    logger.info(f"Total bases R1: {sum_len_r1}")
    logger.info(f"Total bases R2: {sum_len_r2}")
    
    # Step 3: Calculate the proportion p for each sequence
    p_r1 = (whole_genome_size * sampling_cov * 1000000) / sum_len_r1
    p_r2 = (whole_genome_size * sampling_cov * 1000000) / sum_len_r2
    
    # Use the minimum of both proportions to ensure balanced pairs
    p = min(p_r1, p_r2)
    
    # Ensure p is between 0.001 and 1.0
    if p < 0.001:
        logger.warning(f"Sampling proportion {p:.6f} is too low, using minimum 0.001")
        p = 0.001
    elif p > 1.0:
        logger.warning(f"Sampling proportion {p:.6f} exceeds 1.0, using 1.0 (100%)")
        p = 1.0
    
    logger.info(f"Calculated sampling proportion: {p:.6f} ({p*100:.2f}%)")
    logger.info(f"p_r1: {p_r1:.6f}, p_r2: {p_r2:.6f}, using p: {p:.6f}")
    
    # Step 4: Run seqtk sample to subsample the sequences
    # Step 5: Compress the output with gzip
    cmd_sample1 = f"seqtk sample -s 100 {fq1} {p} | gzip > {output1}"
    cmd_sample2 = f"seqtk sample -s 100 {fq2} {p} | gzip > {output2}"
    
    logger.info(f"Running command: {cmd_sample1}")
    subprocess.run(cmd_sample1, shell=True, check=True)
    
    logger.info(f"Running command: {cmd_sample2}")
    subprocess.run(cmd_sample2, shell=True, check=True)
    
    # Verify output files exist and are not empty
    for output_file in [output1, output2]:
        if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
            logger.warning(f"Output file {output_file} is empty or missing. Creating empty file.")
            Path(output_file).touch()  # Ensure file exists
    
    logger.info("Subsampling complete.")

except subprocess.CalledProcessError as e:
    logger.error(f"Error occurred during subprocess execution: {str(e)}")
    # Create empty output files to avoid workflow failure
    Path(output1).touch()
    Path(output2).touch()
    sys.exit(1)
except Exception as e:
    logger.error(f"General error occurred: {str(e)}")
    # Create empty output files to avoid workflow failure
    Path(output1).touch()
    Path(output2).touch()
    sys.exit(1)