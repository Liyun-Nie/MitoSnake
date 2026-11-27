#!/usr/bin/env python3
"""
Extract sequences from FASTA/FASTQ files based on Tiara classification results.
"""

import argparse
import os
import sys
from Bio import SeqIO

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Extract sequences from FASTA/FASTQ based on Tiara classification"
    )
    
    parser.add_argument(
        '--input',
        required=True,
        help='Input FASTA/FASTQ file'
    )
    
    parser.add_argument(
        '--classes',
        required=True,
        help='Tiara classification file'
    )
    
    parser.add_argument(
        '--output',
        required=True,
        help='Output FASTA file'
    )
    
    parser.add_argument(
        '--target',
        default='mitochondrion',
        choices=['mitochondrion', 'plastid', 'both'],
        help='Target class to extract (default: mitochondrion)'
    )
    
    return parser.parse_args()

def read_tiara_classes(class_file):
    """Read Tiara classification file and return dict of sequence IDs -> classes"""
    classes = {}
    with open(class_file) as f:
        # Skip header
        next(f, None)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                seq_id = parts[0]
                seq_class = parts[1].lower()
                classes[seq_id] = seq_class
    return classes

def main():
    """Main function"""
    args = parse_arguments()
    
    # Validate input files
    if not os.path.exists(args.input):
        sys.stderr.write(f"Input file not found: {args.input}\n")
        sys.exit(1)
    
    if not os.path.exists(args.classes):
        sys.stderr.write(f"Classes file not found: {args.classes}\n")
        sys.exit(1)
    
    # Read classifications
    classifications = read_tiara_classes(args.classes)
    
    # Determine target classes
    target_classes = []
    if args.target == 'mitochondrion':
        target_classes = ['mitochondrion', 'mit']
    elif args.target == 'plastid':
        target_classes = ['plastid', 'pla']
    elif args.target == 'both':
        target_classes = ['mitochondrion', 'mit', 'plastid', 'pla']
    
    # Determine input format
    fmt = 'fasta'
    if args.input.endswith(('.fq', '.fastq')):
        fmt = 'fastq'
    
    # Extract and write sequences
    with open(args.output, 'w') as out:
        count = 0
        for record in SeqIO.parse(args.input, fmt):
            seq_id = record.id
            if seq_id in classifications:
                seq_class = classifications[seq_id]
                for target in target_classes:
                    if target in seq_class:
                        SeqIO.write(record, out, 'fasta')
                        count += 1
                        break
    
    sys.stderr.write(f"Extracted {count} {args.target} sequences to {args.output}\n")

if __name__ == "__main__":
    main()
