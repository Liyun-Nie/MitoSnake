#!/usr/bin/env python3
"""
Path validation tool for mitochondrial genome assembly.
This script extracts connections from assembly graph, creates contig pairs,
and evaluates the support for each connection using long reads.
"""

import argparse
import os
import sys
import subprocess
import logging
from collections import defaultdict
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random
import re

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Validate connections in mitochondrial genome assembly using long reads"
    )
    
    parser.add_argument(
        '--gfa',
        required=True,
        help='Assembly graph in GFA format'
    )
    
    parser.add_argument(
        '--fasta',
        required=True,
        help='Assembly contigs in FASTA format'
    )
    
    parser.add_argument(
        '--long-reads',
        required=True,
        help='Long reads in FASTA/FASTQ format'
    )
    
    parser.add_argument(
        '--output-dir',
        required=True,
        help='Output directory'
    )
    
    parser.add_argument(
        '--min-support',
        type=int,
        default=3,
        help='Minimum number of reads supporting a connection (default: 3)'
    )
    
    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help='Number of threads for mapping (default: 4)'
    )
    
    parser.add_argument(
        '--overlap',
        type=int,
        default=500,
        help='Overlap size to include from each contig (default: 500 bp)'
    )
    
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose output'
    )
    
    parser.add_argument(
        '--sampling-rate',
        type=float,
        default=0.05,
        help='Fraction of reads to use for validation (default: 0.05 = 5%)'
    )
    
    return parser.parse_args()

def parse_gfa(gfa_file):
    """Parse GFA file to extract nodes and links"""
    nodes = {}  # node_id -> sequence
    links = []  # list of (from_node, from_orient, to_node, to_orient, cigar)
    
    with open(gfa_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if not parts:
                continue
                
            # Parse segment lines (nodes)
            if parts[0] == 'S':
                node_id = parts[1]
                sequence = parts[2]
                nodes[node_id] = sequence
                
            # Parse link lines (edges)
            elif parts[0] == 'L':
                from_node = parts[1]
                from_orient = parts[2]
                to_node = parts[3]
                to_orient = parts[4]
                cigar = parts[5]
                links.append((from_node, from_orient, to_node, to_orient, cigar))
    
    logging.info(f"Parsed {len(nodes)} nodes and {len(links)} links from GFA file")
    return nodes, links

def parse_fasta(fasta_file):
    """Parse FASTA file to get node sequences"""
    contigs = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        contigs[record.id] = str(record.seq)
    
    logging.info(f"Parsed {len(contigs)} contigs from FASTA file")
    return contigs

def create_contig_pairs(nodes, links, contigs, output_dir, overlap=500):
    """Create FASTA files for each link containing the pair of connected contigs"""
    os.makedirs(output_dir, exist_ok=True)
    pair_files = {}
    
    for i, (from_node, from_orient, to_node, to_orient, cigar) in enumerate(links):
        # Get the sequences
        if from_node in contigs:
            from_seq = contigs[from_node]
        elif from_node in nodes:
            from_seq = nodes[from_node]
        else:
            logging.warning(f"Node {from_node} not found in contigs or GFA, skipping")
            continue
            
        if to_node in contigs:
            to_seq = contigs[to_node]
        elif to_node in nodes:
            to_seq = nodes[to_node]
        else:
            logging.warning(f"Node {to_node} not found in contigs or GFA, skipping")
            continue
        
        # Handle orientation
        if from_orient == '-':
            from_seq = str(Seq(from_seq).reverse_complement())
        if to_orient == '-':
            to_seq = str(Seq(to_seq).reverse_complement())
        
        # Extract overlap region based on CIGAR (simplification, using fixed overlap)
        overlap_size = min(overlap, len(from_seq), len(to_seq))
        
        # Create pair sequence
        from_part = from_seq[-overlap_size:] if len(from_seq) > overlap_size else from_seq
        to_part = to_seq[:overlap_size] if len(to_seq) > overlap_size else to_seq
        pair_seq = from_seq[:-overlap_size] + from_part + to_part + to_seq[overlap_size:]
        
        # Create pair record
        pair_id = f"pair_{i+1}_{from_node}_{from_orient}_to_{to_node}_{to_orient}"
        pair_record = SeqRecord(
            Seq(pair_seq),
            id=pair_id,
            description=f"Connection from {from_node}({from_orient}) to {to_node}({to_orient})"
        )
        
        # Write to file
        pair_file = os.path.join(output_dir, f"{pair_id}.fasta")
        with open(pair_file, 'w') as f:
            SeqIO.write(pair_record, f, "fasta")
        
        # Store file path
        pair_files[(from_node, from_orient, to_node, to_orient)] = {
            'file': pair_file,
            'id': pair_id,
            'junction_start': len(from_seq) - overlap_size,
            'junction_end': len(from_seq) + overlap_size
        }
        
    logging.info(f"Created {len(pair_files)} contig pair files")
    return pair_files

def subsample_reads(input_reads, output_reads, sampling_rate=0.1):
    """对输入读长进行采样，只保留一部分用于验证"""
    records = list(SeqIO.parse(input_reads, "fasta" if input_reads.endswith(('.fa', '.fasta')) else "fastq"))
    sampled_records = random.sample(records, int(len(records) * sampling_rate))
    
    SeqIO.write(sampled_records, output_reads, "fasta")
    return output_reads

def map_reads_to_pairs_lightweight(pair_files, long_reads, output_dir, threads=4):
    """轻量级版本的读长比对，只统计跨越连接的读长数量，不保存完整BAM"""
    os.makedirs(os.path.join(output_dir, "results"), exist_ok=True)
    junction_support = {}
    
    for pair_key, pair_info in pair_files.items():
        pair_file = pair_info['file']
        pair_id = pair_info['id']
        junction_start = pair_info['junction_start']
        junction_end = pair_info['junction_end']
        
        # 定义连接区域
        junction_region = (junction_start, junction_end)
        
        # 使用minimap2直接输出SAM格式并处理
        read_format = "map-ont" if long_reads.endswith((".fq", ".fastq", ".fq.gz", ".fastq.gz")) else "map-pb"
        
        # 修改minimap2命令，使用更快速的设置
        minimap2_cmd = [
            "minimap2",
            "-x", read_format,
            "-m", "40",  # 最小匹配长度
            "--secondary=no",  # 不输出次优比对
            "-t", str(threads),
            pair_file,
            long_reads
        ]
        
        try:
            # 运行minimap2并直接处理输出
            proc = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            
            # 统计跨越连接的读长
            supporting_reads = set()
            
            # 逐行读取SAM输出
            for line in proc.stdout:
                if line.startswith('@'):  # 跳过SAM头
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) < 11 or parts[2] == '*':  # 未比对上的读长
                    continue
                
                read_id = parts[0]
                pos = int(parts[3])
                cigar = parts[5]
                
                # 计算读长覆盖区域
                read_end = pos
                expanded_cigar = re.findall(r'(\d+)([MIDNSHPX=])', cigar)
                for length, op in expanded_cigar:
                    length = int(length)
                    if op in 'M=X':  # 这些操作会消耗参考序列
                        read_end += length
                
                # 检查读长是否跨越连接区域
                if pos <= junction_start and read_end >= junction_end:
                    supporting_reads.add(read_id)
            
            # 等待进程完成
            proc.stdout.close()
            proc.wait()
            
            # 存储结果
            junction_support[pair_key] = {
                'count': len(supporting_reads),
                'reads': list(supporting_reads)
            }
            
            logging.info(f"Junction {pair_id}: {len(supporting_reads)} supporting reads")
            
        except Exception as e:
            logging.error(f"Error mapping reads to {pair_id}: {e}")
            junction_support[pair_key] = {
                'count': 0,
                'reads': []
            }
    
    # 将结果写入文件
    results_file = os.path.join(output_dir, "results", "junction_support.tsv")
    with open(results_file, 'w') as f:
        f.write("From_Node\tFrom_Orient\tTo_Node\tTo_Orient\tSupport_Count\tSupporting_Reads\n")
        
        for (from_node, from_orient, to_node, to_orient), support in junction_support.items():
            f.write(f"{from_node}\t{from_orient}\t{to_node}\t{to_orient}\t{support['count']}\t{','.join(support['reads'])}\n")
    
    return junction_support

def analyze_junction_support(pair_files, mapping_results, output_dir, min_support=3):
    """Analyze read support for each junction"""
    os.makedirs(os.path.join(output_dir, "results"), exist_ok=True)
    junction_support = {}
    
    for pair_key, bam_file in mapping_results.items():
        pair_info = pair_files[pair_key]
        pair_id = pair_info['id']
        junction_start = pair_info['junction_start']
        junction_end = pair_info['junction_end']
        
        # Define junction region with a buffer
        region = f"{pair_id}:{junction_start-50}-{junction_end+50}"
        
        # Count reads spanning the junction
        count_cmd = [
            "samtools", "view", 
            "-c",  # Count
            bam_file, 
            region
        ]
        
        try:
            result = subprocess.run(count_cmd, capture_output=True, text=True, check=True)
            support_count = int(result.stdout.strip())
            
            # Get supporting read IDs
            view_cmd = [
                "samtools", "view",
                bam_file,
                region
            ]
            
            result = subprocess.run(view_cmd, capture_output=True, text=True, check=True)
            supporting_reads = [line.split("\t")[0] for line in result.stdout.strip().split("\n") if line]
            
            # Store results
            junction_support[pair_key] = {
                'count': support_count,
                'reads': supporting_reads if supporting_reads[0] != '' else [],
                'is_supported': support_count >= min_support
            }
            
        except Exception as e:
            logging.error(f"Error analyzing junction support for {pair_id}: {e}")
            junction_support[pair_key] = {
                'count': 0,
                'reads': [],
                'is_supported': False,
                'error': str(e)
            }
    
    # Write results to file
    results_file = os.path.join(output_dir, "results", "junction_support.tsv")
    with open(results_file, 'w') as f:
        f.write("From_Node\tFrom_Orient\tTo_Node\tTo_Orient\tSupport_Count\tIs_Supported\tSupporting_Reads\n")
        
        for (from_node, from_orient, to_node, to_orient), support in junction_support.items():
            f.write(f"{from_node}\t{from_orient}\t{to_node}\t{to_orient}\t{support['count']}\t{support['is_supported']}\t{','.join(support['reads'])}\n")
    
    # Write summary
    summary_file = os.path.join(output_dir, "results", "summary.txt")
    with open(summary_file, 'w') as f:
        total_junctions = len(junction_support)
        supported_junctions = sum(1 for support in junction_support.values() if support['is_supported'])
        
        f.write(f"Total junctions analyzed: {total_junctions}\n")
        f.write(f"Supported junctions (>= {min_support} reads): {supported_junctions} ({supported_junctions/total_junctions*100:.1f}%)\n")
        f.write(f"Unsupported junctions: {total_junctions - supported_junctions} ({(total_junctions - supported_junctions)/total_junctions*100:.1f}%)\n\n")
        
        f.write("Supported Junctions:\n")
        for (from_node, from_orient, to_node, to_orient), support in junction_support.items():
            if support['is_supported']:
                f.write(f"  {from_node}({from_orient}) -> {to_node}({to_orient}): {support['count']} reads\n")
        
        f.write("\nUnsupported Junctions:\n")
        for (from_node, from_orient, to_node, to_orient), support in junction_support.items():
            if not support['is_supported']:
                f.write(f"  {from_node}({from_orient}) -> {to_node}({to_orient}): {support['count']} reads\n")
    
    logging.info(f"Analysis complete. Results written to {results_file} and {summary_file}")
    return junction_support

def main():
    args = parse_args()
    
    try:
        # Create output directory
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Step 1: Parse input files
        logging.info("Parsing GFA file...")
        nodes, links = parse_gfa(args.gfa)
        
        logging.info("Parsing FASTA file...")
        contigs = parse_fasta(args.fasta)
        
        # Step 2: Create contig pairs
        logging.info("Creating contig pairs...")
        pair_files = create_contig_pairs(
            nodes, links, contigs, 
            os.path.join(args.output_dir, "pairs"),
            overlap=args.overlap
        )
        
        # Step 3: Map reads to contig pairs
        logging.info("Mapping reads to contig pairs...")
        sampled_reads = os.path.join(args.output_dir, "sampled_reads.fa")
        subsample_reads(args.long_reads, sampled_reads, args.sampling_rate)
        mapping_results = map_reads_to_pairs_lightweight(
            pair_files, sampled_reads, 
            args.output_dir, threads=args.threads
        )
        
        # Step 4: Analyze junction support
        logging.info("Analyzing junction support...")
        junction_support = analyze_junction_support(
            pair_files, mapping_results, 
            args.output_dir, min_support=args.min_support
        )
        
        logging.info("Path validation complete!")
        
    except Exception as e:
        logging.error(f"Error in path validation: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
