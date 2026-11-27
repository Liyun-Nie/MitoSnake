#!/usr/bin/env python3
"""
Iterative assembly pipeline control script.
Supports multi-round iterative assembly, using results from previous round as seed for next round.
"""

import argparse
import os
import sys
import subprocess
import shutil
import logging
import tempfile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def setup_logging(verbose=False):
    """Setup logging configuration"""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stderr)
        ]
    )

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Iterative assembly pipeline for plant mitochondrial genomes"
    )
    
    parser.add_argument(
        '--seed',
        required=True,
        help='Initial seed file (FASTA)'
    )
    
    parser.add_argument(
        '--long-reads',
        required=True,
        help='Long reads file (FASTA/FASTQ)'
    )
    
    parser.add_argument(
        '--output-dir',
        required=True,
        help='Base output directory'
    )
    
    parser.add_argument(
        '--genome-size',
        type=str,
        default='1m',
        help='Estimated mitochondrial genome size (default: 1m)'
    )
    
    parser.add_argument(
        '--longread-type',
        choices=['ont', 'clr', 'hifi'],
        default='ont',
        help='Type of long reads (default: ont)'
    )
    
    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help='Number of threads (default: 4)'
    )
    
    parser.add_argument(
        '--min-match',
        type=int,
        default=2000,
        help='Minimum number of matching bases for filtering (default: 2000)'
    )
    
    parser.add_argument(
        '--min-coverage',
        type=float,
        default=0.3,
        help='Minimum coverage for filtering (default: 0.3)'
    )
    
    parser.add_argument(
        '--script-dir',
        required=True,
        help='Path to script directory'
    )
    
    parser.add_argument(
        '--rounds',
        type=int,
        default=2,
        help='Number of iterative assembly rounds (default: 2)'
    )
    
    parser.add_argument(
        '--use-dl-filter',
        action='store_true',
        help='Use deep learning filter'
    )
    
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )
    
    return parser.parse_args()

def run_round(seed_file, long_reads, output_dir, round_num, genome_size, longread_type, 
             threads, min_match, min_coverage, script_dir, use_dl_filter):
    """Run single round of assembly"""
    round_dir = os.path.join(output_dir, f"flye_out_round{round_num}")
    os.makedirs(round_dir, exist_ok=True)
    
    # Step 1: Direct mapping of long reads to seed sequences
    mapping_dir = os.path.join(round_dir, "mapping")
    os.makedirs(mapping_dir, exist_ok=True)
    
    paf_file = os.path.join(mapping_dir, "mapped_reads.paf")
    
    # Choose mapping parameters based on long read type
    map_type = "map-ont" if longread_type == "ont" else "map-pb"
    
    # Execute minimap2 mapping
    cmd = [
        'minimap2',
        '-x', map_type,
        '-t', str(threads),
        seed_file,
        long_reads,
        '-o', paf_file
    ]
    
    logging.info(f"Round {round_num}: Running mapping: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    # Step 2: Filter long reads
    filtered_reads = os.path.join(round_dir, "filtered_longreads.fa")
    
    # use strict filter parameters for subsequent rounds
    if round_num > 1:
        # use strict filter parameters for subsequent rounds
        adjusted_min_coverage = max(0.5, min_coverage)  
    else:
        adjusted_min_coverage = min_coverage
    
    # use adjusted parameters when calling filter_longreads.py
    filter_cmd = [
        'python',
        os.path.join(script_dir, 'filter_longreads.py'),
        '--paf', paf_file,
        '--reads', long_reads,
        '--output', filtered_reads,
        '--min-match', str(min_match),
        '--min-coverage', str(adjusted_min_coverage)  # use adjusted parameters
    ]
    
    logging.info(f"Round {round_num}: Filtering reads: {' '.join(filter_cmd)}")
    subprocess.run(filter_cmd, check=True)
    
    # Step 3: Optional Tiara classification
    if use_dl_filter:
        # define tiara output file      
        classified_reads = os.path.join(round_dir, "mitochondrion_filtered_longreads.fa")
        plastid_output = os.path.join(round_dir, "plastid_filtered_longreads.fa")
        tiara_out_file = os.path.join(round_dir, "tiara_out.txt")
        tiara_log_out = os.path.join(round_dir, "log_tiara_out.txt")
        
        # directly call tiara command
        tiara_cmd = [
            "conda", "run", "-n", "mitosnake_tiara", 
            "tiara",
            "-i", filtered_reads,
            "-o", tiara_out_file,
            "--tf", "mit", "pla",
            "-t", str(threads)
        ]
        
        logging.info(f"Round {round_num}: Running Tiara: {' '.join(tiara_cmd)}")
        
        try:
            # create log file
            tiara_log = os.path.join(round_dir, "tiara_log.txt")
            
            # execute tiara command
            with open(tiara_log, 'w') as log_file:
                with open(tiara_log_out, 'w') as out_file:
                    process = subprocess.Popen(
                        tiara_cmd,
                        stdout=out_file,
                        stderr=log_file
                    )
                    process.wait()
                    
                    if process.returncode != 0:
                        raise Exception(f"Tiara exited with code {process.returncode}")
            
            # check if tiara has generated mitochondrial sequences
            if not os.path.exists(classified_reads) or os.path.getsize(classified_reads) == 0:
                logging.warning(f"Round {round_num}: No mitochondrial sequences found or file is empty")
                logging.warning(f"Round {round_num}: Using all filtered reads")
                shutil.copy2(filtered_reads, classified_reads)
        except Exception as e:
            logging.warning(f"Round {round_num}: Tiara failed: {str(e)}")
            logging.warning(f"Round {round_num}: Using all filtered reads")
            shutil.copy2(filtered_reads, classified_reads)
    else:
        # skip tiara classification
        classified_reads = os.path.join(round_dir, "mitochondrion_filtered_longreads.fa")
        shutil.copy2(filtered_reads, classified_reads)
    
    # Step 4: Run assembly using Flye directly
    assembly_dir = os.path.join(round_dir, "assembly")
    os.makedirs(assembly_dir, exist_ok=True)
    
    # Get Flye read type parameter
    read_type_param = '--nano-raw' if longread_type == 'ont' else ('--pacbio-hifi' if longread_type == 'hifi' else '--pacbio-raw')
    
    # Use fixed mitochondrial genome size 1m instead of the passed parameter
    mito_genome_size = '1m'
    logging.info(f"Round {round_num}: Using mitochondrial genome size: {mito_genome_size} for Flye assembly")
    
    flye_cmd = [
        'flye',
        read_type_param, classified_reads,
        '--out-dir', assembly_dir,
        '--genome-size', mito_genome_size,  # Use fixed value 1m
        '--threads', str(threads)
    ]
    
    logging.info(f"Round {round_num}: Running Flye assembly: {' '.join(flye_cmd)}")
    try:
        subprocess.run(flye_cmd, check=True)
    except:
        logging.error(f"Round {round_num}: Assembly failed")
        # Create empty files to ensure workflow can continue
        open(os.path.join(assembly_dir, "assembly.fasta"), 'w').close()
        open(os.path.join(assembly_dir, "assembly_graph.gfa"), 'w').close()
    
    # Return assembly result file paths
    return {
        "assembly": os.path.join(assembly_dir, "assembly.fasta"),
        "graph": os.path.join(assembly_dir, "assembly_graph.gfa")
    }

def create_seed_from_assembly(assembly_file, output_file, round_num):
    """Create seed file for next round from assembly results"""
    if not os.path.exists(assembly_file) or os.path.getsize(assembly_file) == 0:
        logging.warning(f"Assembly file {assembly_file} is empty or doesn't exist")
        # Return None to indicate seed creation failed
        return None
    
    # Read assembly results and rename sequences
    records = []
    for i, record in enumerate(SeqIO.parse(assembly_file, "fasta")):
        new_id = f"round{round_num}_contig_{i+1}"
        new_desc = f"from_round_{round_num}_assembly original_id={record.id}"
        new_record = SeqRecord(
            seq=record.seq,
            id=new_id,
            description=new_desc
        )
        records.append(new_record)
    
    # If no sequences found, return None
    if not records:
        logging.warning(f"No sequences found in {assembly_file}")
        return None
    
    # Write new seed file
    with open(output_file, 'w') as f:
        SeqIO.write(records, f, "fasta")
    
    logging.info(f"Created seed file for next round: {output_file}")
    return output_file

# improve finalize_results function to ensure correct file and symlink creation
def finalize_results(rounds_results, output_dir):
    # get last round results
    last_round = max(rounds_results.keys())
    last_assembly = rounds_results[last_round]["assembly"]
    last_graph = rounds_results[last_round]["graph"]
    
    # copy to final output
    final_assembly = os.path.join(output_dir, "final_assembly.fasta")
    final_graph = os.path.join(output_dir, "final_assembly.gfa")
    
    # check if files exist
    if os.path.exists(last_assembly) and os.path.getsize(last_assembly) > 0:
        shutil.copy2(last_assembly, final_assembly)
    else:
        open(final_assembly, 'w').close()
        logging.warning("Created empty final assembly file")
    
    if os.path.exists(last_graph) and os.path.getsize(last_graph) > 0:
        shutil.copy2(last_graph, final_graph)
    else:
        open(final_graph, 'w').close()
        logging.warning("Created empty final graph file")
    
    # create symlink or copy file to ensure correct files in each round directory
    for round_num in rounds_results:
        round_dir = os.path.join(output_dir, f"flye_out_round{round_num}")
        assembly_dir = os.path.join(round_dir, "assembly")
        
        if os.path.exists(assembly_dir):
            # source file path
            src_assembly = os.path.join(assembly_dir, "assembly.fasta")
            src_graph = os.path.join(assembly_dir, "assembly_graph.gfa")
            
            # destination file path
            dest_assembly = os.path.join(round_dir, "assembly.fasta")
            dest_graph = os.path.join(round_dir, "assembly_graph.gfa")
            
            # create file or link
            for src, dest in [(src_assembly, dest_assembly), (src_graph, dest_graph)]:
                if os.path.exists(src) and not os.path.exists(dest):
                    logging.info(f"Creating file link from {src} to {dest}")
                    try:
                        # try creating hard link
                        os.link(src, dest)
                    except OSError:
                        try:
                            # if hard link fails, try creating symlink
                            os.symlink(src, dest)
                        except OSError:
                            # if symlink fails, copy file
                            shutil.copy2(src, dest)
                            logging.info(f"Copied file from {src} to {dest}")
    
    # create file containing all rounds results 
    all_rounds_file = os.path.join(output_dir, "all_rounds_assemblies.fasta")
    with open(all_rounds_file, 'w') as out:
        for round_num in sorted(rounds_results.keys()):
            assembly_file = rounds_results[round_num]["assembly"]
            if os.path.exists(assembly_file) and os.path.getsize(assembly_file) > 0:
                out.write(f"# Round {round_num} assembly\n")
                with open(assembly_file) as f:
                    out.write(f.read())
                out.write("\n\n")
    
    logging.info(f"Finalized results in {output_dir}")
    return final_assembly, final_graph, all_rounds_file

def main():
    """Main function"""
    args = parse_arguments()
    setup_logging(args.verbose)
    
    try:
        # Validate input files
        if not os.path.exists(args.seed):
            logging.error(f"Seed file not found: {args.seed}")
            sys.exit(1)
        
        if not os.path.exists(args.long_reads):
            logging.error(f"Long reads file not found: {args.long_reads}")
            sys.exit(1)
        
        # Create base output directory
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Store results from each round
        rounds_results = {}
        
        # Set initial seed
        current_seed = args.seed
        
        # Run specified number of iterative assembly rounds
        for round_num in range(1, args.rounds + 1):
            logging.info(f"Starting assembly round {round_num}")
            
            # Run current round
            result = run_round(
                seed_file=current_seed,
                long_reads=args.long_reads,
                output_dir=args.output_dir,
                round_num=round_num,
                genome_size=args.genome_size,
                longread_type=args.longread_type,
                threads=args.threads,
                min_match=args.min_match,
                min_coverage=args.min_coverage,
                script_dir=args.script_dir,
                use_dl_filter=args.use_dl_filter
            )
            
            # Save results
            rounds_results[round_num] = result
            
            # If not the last round, prepare seed for next round
            if round_num < args.rounds:
                next_seed = os.path.join(args.output_dir, f"seed_for_round{round_num+1}.fasta")
                seed_created = create_seed_from_assembly(result["assembly"], next_seed, round_num)
                
                # If seed creation failed, use current seed to continue
                if seed_created:
                    current_seed = next_seed
                else:
                    logging.warning(f"Could not create seed from round {round_num}, using previous seed")
        
        # Finalize results
        final_assembly, final_graph, all_rounds = finalize_results(rounds_results, args.output_dir)
        
        logging.info("Iterative assembly pipeline completed successfully")
        logging.info(f"Final assembly: {final_assembly}")
        logging.info(f"Final graph: {final_graph}")
        logging.info(f"All rounds assembly: {all_rounds}")
        
    except Exception as e:
        logging.error(f"Error in iterative pipeline: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()