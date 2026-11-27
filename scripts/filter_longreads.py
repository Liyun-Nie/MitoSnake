#!/usr/bin/env python3
"""
Filter long reads based on mapping quality for mitochondrial genome assembly.

This script processes PAF mapping files to select reads that meet quality criteria
for mitochondrial genome assembly.

Author: MITO-SNAKE Development Team
Version: 0.3.1
"""

import argparse
import sys
from collections import defaultdict
import os
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', stream=sys.stderr)

def parse_args():
    parser = argparse.ArgumentParser(description='Filter long reads based on mapping results')
    parser.add_argument('--paf', required=True, help='Input PAF file')
    parser.add_argument('--reads', required=True, help='Input reads file (FASTQ/FASTA)')
    parser.add_argument('--output', required=True, help='Output filtered reads file')
    parser.add_argument('--min-match', type=int, default=2000,
                      help='Minimum number of matching bases in the alignment')
    parser.add_argument('--min-coverage', type=float, default=0.3,
                      help='Minimum coverage ratio (matching bases / query length)')
    parser.add_argument('--fallback-count', type=int, default=10,
                      help='Number of longest reads to include when no reads pass filters')
    return parser.parse_args()

def process_paf(paf_file, min_match, min_coverage):
    """Process PAF file and return set of read IDs that pass filters."""
    passed_reads = set()
    alignment_stats = {
        'total': 0,
        'passed': 0,
        'failed_match': 0,
        'failed_coverage': 0
    }
    
    # Check if PAF file is empty
    if os.path.getsize(paf_file) == 0:
        logging.warning("PAF file is empty. No alignments found.")
        logging.warning("Possible reasons: 1) Reference seed sequence doesn't match long reads; 2) Long read quality is too low; 3) Alignment parameters are too strict")
        logging.info("Suggestion: Try different reference seeds or decrease filtering thresholds (min_match and min_coverage)")
        return passed_reads, alignment_stats
    
    with open(paf_file) as f:
        for line in f:
            alignment_stats['total'] += 1
            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue
                
            # Parse required fields
            query_name = fields[0]
            query_length = int(fields[1])
            matching_bases = int(fields[9])
            
            # Calculate coverage ratio
            coverage_ratio = matching_bases / query_length
            
            # Record filtering statistics
            if matching_bases < min_match:
                alignment_stats['failed_match'] += 1
            elif coverage_ratio < min_coverage:
                alignment_stats['failed_coverage'] += 1
            else:
                alignment_stats['passed'] += 1
                passed_reads.add(query_name)
    
    # Log detailed statistics
    if alignment_stats['total'] > 0:
        logging.info(f"Total alignments: {alignment_stats['total']}")
        logging.info(f"Alignments passing filters: {alignment_stats['passed']} ({alignment_stats['passed']/alignment_stats['total']*100:.1f}%)")
        logging.info(f"Alignments failing on match count: {alignment_stats['failed_match']} ({alignment_stats['failed_match']/alignment_stats['total']*100:.1f}%)")
        logging.info(f"Alignments failing on coverage ratio: {alignment_stats['failed_coverage']} ({alignment_stats['failed_coverage']/alignment_stats['total']*100:.1f}%)")

    return passed_reads, alignment_stats

def detect_file_format(filename):
    """Detect if file is FASTA or FASTQ format."""
    with open(filename) as f:
        first_char = f.read(1)
        if first_char == '>':
            return 'fasta'
        elif first_char == '@':
            return 'fastq'
        else:
            raise ValueError(f"Unable to detect format of {filename}. File must start with '>' (FASTA) or '@' (FASTQ).")

def filter_reads(reads_file, output_file, passed_reads, fallback_count=10):
    """Write filtered reads to output file."""
    from Bio import SeqIO
    import warnings
    from Bio import BiopythonWarning
    
    # Suppress Biopython warnings
    warnings.simplefilter('ignore', BiopythonWarning)
    
    # Detect input format
    input_format = detect_file_format(reads_file)
    
    # Determine output format (always output as FASTA)
    output_format = 'fasta'
    
    # Count filtered reads
    total_reads = 0
    passed_count = 0
    
    # Read all records into memory to select longest reads if none pass filtering
    all_records = list(SeqIO.parse(reads_file, input_format))
    total_reads = len(all_records)
    
    with open(output_file, 'w') as out:
        # First try to write reads that passed filtering
        for record in all_records:
            if record.id in passed_reads:
                passed_count += 1
                SeqIO.write(record, out, output_format)
        
        # If no reads passed filtering, select the longest reads as fallback
        if passed_count == 0:
            logging.warning(f"No reads passed filtering criteria. Using the {fallback_count} longest reads as fallback.")
            logging.warning("This is a fallback mechanism, results may not be optimal. Consider adjusting filtering parameters and rerunning.")
            
            # Sort by length and select longest reads
            all_records.sort(key=lambda x: len(x.seq), reverse=True)
            top_records = all_records[:fallback_count]
            
            for record in top_records:
                passed_count += 1
                SeqIO.write(record, out, output_format)
            
            logging.info(f"Selected {passed_count} longest reads to continue the pipeline.")
    
    # Print summary
    logging.info(f"Total reads: {total_reads}")
    logging.info(f"Reads passing filters: {passed_count}")
    logging.info(f"Reads filtered out: {total_reads - passed_count}")
    
    return passed_count

def main():
    args = parse_args()
    
    # Validate input files
    if not os.path.exists(args.paf):
        logging.error(f"Error: PAF file not found: {args.paf}")
        sys.exit(1)
    if not os.path.exists(args.reads):
        logging.error(f"Error: Reads file not found: {args.reads}")
        sys.exit(1)
    
    try:
        # Process PAF file and get passed read IDs
        logging.info("Processing PAF file...")
        passed_reads, stats = process_paf(args.paf, args.min_match, args.min_coverage)
        
        # Filter and write reads
        logging.info("Filtering reads...")
        passed_count = filter_reads(args.reads, args.output, passed_reads, args.fallback_count)
        
        if passed_count > 0:
            logging.info("Filtering completed successfully.")
        else:
            logging.error("No reads retained after filtering, which should not happen (due to fallback mechanism).")
            sys.exit(1)
            
    except Exception as e:
        logging.error(f"Error: {str(e)}")
        sys.exit(1)

if __name__ == '__main__':
    main()

