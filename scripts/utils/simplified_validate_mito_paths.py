#!/usr/bin/env python3
"""
轻量级路径验证工具 - 只计数支持读长
"""
import argparse
import os
import sys
import subprocess
import logging
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_args():
    parser = argparse.ArgumentParser(description="轻量级路径验证工具")
    parser.add_argument('--gfa', required=True, help='GFA格式的组装图')
    parser.add_argument('--fasta', required=True, help='FASTA格式的组装序列')
    parser.add_argument('--long-reads', required=True, help='长读长数据')
    parser.add_argument('--output-dir', required=True, help='输出目录')
    parser.add_argument('--threads', type=int, default=4, help='线程数')
    parser.add_argument('--overlap', type=int, default=500, help='连接区域长度')
    parser.add_argument('--sampling-rate', type=float, default=0.05, help='读长采样率')
    return parser.parse_args()

def parse_gfa(gfa_file):
    nodes = {}
    links = []
    with open(gfa_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if not parts:
                continue
            if parts[0] == 'S':
                node_id = parts[1]
                sequence = parts[2]
                nodes[node_id] = sequence
            elif parts[0] == 'L':
                links.append((parts[1], parts[2], parts[3], parts[4], parts[5]))
    return nodes, links

def parse_fasta(fasta_file):
    contigs = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        contigs[record.id] = str(record.seq)
    return contigs

def subsample_reads(input_reads, output_reads, sampling_rate=0.05):
    """对读长进行采样"""
    import random
    
    # 确定输入格式
    fmt = "fasta"
    if input_reads.endswith(('.fq', '.fastq', '.fq.gz', '.fastq.gz')):
        fmt = "fastq"
    
    # 读取所有记录
    records = list(SeqIO.parse(input_reads, fmt))
    total = len(records)
    
    # 采样
    sample_size = int(total * sampling_rate)
    sampled = random.sample(records, sample_size)
    
    # 写入输出文件
    SeqIO.write(sampled, output_reads, "fasta")
    logging.info(f"采样了 {sample_size}/{total} ({sampling_rate*100:.1f}%) 条读长")
    
    return output_reads

def validate_connections(nodes, links, contigs, long_reads, output_dir, threads=4, overlap=500):
    """验证所有连接"""
    tmp_dir = os.path.join(output_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)
    
    results = {}
    
    # 对每个连接进行验证
    for i, (from_node, from_orient, to_node, to_orient, cigar) in enumerate(links):
        # 获取序列
        if from_node in contigs:
            from_seq = contigs[from_node]
        elif from_node in nodes:
            from_seq = nodes[from_node]
        else:
            logging.warning(f"节点 {from_node} 未找到，跳过")
            continue
            
        if to_node in contigs:
            to_seq = contigs[to_node]
        elif to_node in nodes:
            to_seq = nodes[to_node]
        else:
            logging.warning(f"节点 {to_node} 未找到，跳过")
            continue
        
        # 处理方向
        if from_orient == '-':
            from_seq = str(Seq(from_seq).reverse_complement())
        if to_orient == '-':
            to_seq = str(Seq(to_seq).reverse_complement())
        
        # 提取重叠区域
        overlap_size = min(overlap, len(from_seq), len(to_seq))
        
        # 创建连接序列
        from_part = from_seq[-overlap_size:] if len(from_seq) > overlap_size else from_seq
        to_part = to_seq[:overlap_size] if len(to_seq) > overlap_size else to_seq
        
        # 使用更短的序列来减少比对时间
        junction_seq = from_seq[-2000:] if len(from_seq) > 2000 else from_seq
        junction_seq += to_seq[:2000] if len(to_seq) > 2000 else to_seq
        
        # 创建连接ID
        junction_id = f"junction_{i+1}_{from_node}_{from_orient}_to_{to_node}_{to_orient}"
        
        # 写入临时文件
        junction_file = os.path.join(tmp_dir, f"{junction_id}.fa")
        with open(junction_file, 'w') as f:
            f.write(f">{junction_id}\n{junction_seq}\n")
        
        # 计算连接位置
        junction_pos = len(from_seq)
        
        # 使用minimap2快速比对并统计支持度
        cmd = [
            "minimap2",
            "-c",  # 输出PAF格式
            "-x", "map-pb",
            "--secondary=no",
            "-t", str(threads),
            junction_file,
            long_reads
        ]
        
        logging.info(f"验证连接 {junction_id}")
        
        try:
            # 运行minimap2
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            
            # 统计跨越连接的读长
            supporting_reads = set()
            
            # 从PAF输出中统计
            for line in proc.stdout:
                parts = line.strip().split('\t')
                if len(parts) < 10:
                    continue
                
                # PAF格式: query_name, query_length, query_start, query_end, strand, 
                #         target_name, target_length, target_start, target_end, matches, block_length, mapq
                
                read_id = parts[0]
                target_start = int(parts[7])
                target_end = int(parts[8])
                
                # 检查读长是否跨越连接
                junction_start = junction_pos - overlap_size
                junction_end = junction_pos + overlap_size
                
                if target_start <= junction_start and target_end >= junction_end:
                    supporting_reads.add(read_id)
            
            # 等待进程完成
            proc.stdout.close()
            proc.wait()
            
            # 存储结果
            results[(from_node, from_orient, to_node, to_orient)] = {
                'count': len(supporting_reads),
                'reads': list(supporting_reads)
            }
            
            logging.info(f"连接 {junction_id}: {len(supporting_reads)} 条支持读长")
            
        except Exception as e:
            logging.error(f"验证连接 {junction_id} 出错: {e}")
            results[(from_node, from_orient, to_node, to_orient)] = {
                'count': 0,
                'reads': []
            }
    
    # 写入结果
    results_file = os.path.join(output_dir, "connection_support.tsv")
    with open(results_file, 'w') as f:
        f.write("From_Node\tFrom_Orient\tTo_Node\tTo_Orient\tSupport_Count\tSupporting_Reads\n")
        
        for (from_node, from_orient, to_node, to_orient), support in results.items():
            f.write(f"{from_node}\t{from_orient}\t{to_node}\t{to_orient}\t{support['count']}\t{','.join(support['reads'][:10])}\n")
    
    # 写入摘要
    summary_file = os.path.join(output_dir, "summary.txt")
    with open(summary_file, 'w') as f:
        total = len(results)
        supported = sum(1 for support in results.values() if support['count'] >= 3)
        
        f.write(f"总连接数: {total}\n")
        f.write(f"有支持的连接数 (≥3 reads): {supported} ({supported/total*100:.1f}%)\n\n")
        
        f.write("==== 有支持的连接 ====\n")
        for (from_node, from_orient, to_node, to_orient), support in sorted(
                results.items(), key=lambda x: x[1]['count'], reverse=True):
            if support['count'] >= 3:
                f.write(f"{from_node}({from_orient}) -> {to_node}({to_orient}): {support['count']} reads\n")
        
        f.write("\n==== 无支持的连接 ====\n")
        for (from_node, from_orient, to_node, to_orient), support in sorted(
                results.items(), key=lambda x: x[1]['count']):
            if support['count'] < 3:
                f.write(f"{from_node}({from_orient}) -> {to_node}({to_orient}): {support['count']} reads\n")
    
    return results

def generate_html_report(results, output_dir, gfa_file):
    """生成HTML格式的可视化报告"""
    # 提取节点和连接信息
    nodes = set()
    for from_node, _, to_node, _, _ in results.keys():
        nodes.add(from_node)
        nodes.add(to_node)
    
    # 创建节点-连接图数据
    nodes_data = []
    links_data = []
    
    for node in sorted(nodes):
        nodes_data.append({'id': node, 'group': 1})
    
    for (from_node, from_orient, to_node, to_orient), support in results.items():
        links_data.append({
            'source': from_node,
            'target': to_node,
            'value': support['count'],
            'supported': support['count'] >= 3,
            'from_orient': from_orient,
            'to_orient': to_orient,
            'count': support['count']
        })
    
    # 创建HTML报告
    html_file = os.path.join(output_dir, "connection_report.html")
    with open(html_file, 'w') as f:
        f.write(f'''
        <!DOCTYPE html>
        <html>
        <head>
            <meta charset="utf-8">
            <title>连接支持度报告</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                .container {{ display: flex; }}
                .graph {{ width: 800px; height: 600px; border: 1px solid #ccc; }}
                .info {{ margin-left: 20px; }}
                table {{ border-collapse: collapse; }}
                th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                th {{ background-color: #f2f2f2; }}
                tr.supported {{ background-color: #d4edda; }}
                tr.unsupported {{ background-color: #f8d7da; }}
            </style>
            <script src="https://d3js.org/d3.v5.min.js"></script>
        </head>
        <body>
            <h1>线粒体基因组连接支持度报告</h1>
            <div class="container">
                <div class="graph" id="graph"></div>
                <div class="info">
                    <h2>连接支持度摘要</h2>
                    <p>总连接数: {len(results)}</p>
                    <p>有支持的连接数 (≥3 reads): {sum(1 for s in results.values() if s['count'] >= 3)}</p>
                    <p>无支持的连接数: {sum(1 for s in results.values() if s['count'] < 3)}</p>
                    
                    <h2>所有连接</h2>
                    <table id="connections">
                        <tr>
                            <th>从节点</th>
                            <th>方向</th>
                            <th>到节点</th>
                            <th>方向</th>
                            <th>支持读长数</th>
                            <th>状态</th>
                        </tr>
                        {''.join(f"""
                        <tr class="{'supported' if support['count'] >= 3 else 'unsupported'}">
                            <td>{from_node}</td>
                            <td>{from_orient}</td>
                            <td>{to_node}</td>
                            <td>{to_orient}</td>
                            <td>{support['count']}</td>
                            <td>{'支持' if support['count'] >= 3 else '不支持'}</td>
                        </tr>
                        """ for (from_node, from_orient, to_node, to_orient), support in
                        sorted(results.items(), key=lambda x: x[1]['count'], reverse=True))}
                    </table>
                </div>
            </div>
            
            <script>
            // 使用D3.js创建交互式图形
            const nodes = {str(nodes_data).replace("'", '"')};
            const links = {str(links_data).replace("'", '"')};
            
            const width = 800;
            const height = 600;
            
            const svg = d3.select("#graph")
                .append("svg")
                .attr("width", width)
                .attr("height", height);
            
            // 创建力导向图
            const simulation = d3.forceSimulation(nodes)
                .force("link", d3.forceLink(links).id(d => d.id).distance(100))
                .force("charge", d3.forceManyBody().strength(-300))
                .force("center", d3.forceCenter(width / 2, height / 2));
            
            // 绘制连接
            const link = svg.append("g")
                .selectAll("line")
                .data(links)
                .enter().append("line")
                .attr("stroke-width", d => Math.sqrt(d.value))
                .attr("stroke", d => d.supported ? "#28a745" : "#dc3545");
            
            // 绘制节点
            const node = svg.append("g")
                .selectAll("circle")
                .data(nodes)
                .enter().append("circle")
                .attr("r", 10)
                .attr("fill", "#1f77b4")
                .call(d3.drag()
                    .on("start", dragstarted)
                    .on("drag", dragged)
                    .on("end", dragended));
            
            // 添加节点标签
            const text = svg.append("g")
                .selectAll("text")
                .data(nodes)
                .enter().append("text")
                .text(d => d.id)
                .attr("font-size", 10)
                .attr("dx", 12)
                .attr("dy", 4);
            
            // 更新位置
            simulation.on("tick", () => {
                link
                    .attr("x1", d => d.source.x)
                    .attr("y1", d => d.source.y)
                    .attr("x2", d => d.target.x)
                    .attr("y2", d => d.target.y);
                
                node
                    .attr("cx", d => d.x)
                    .attr("cy", d => d.y);
                
                text
                    .attr("x", d => d.x)
                    .attr("y", d => d.y);
            });
            
            // 拖拽函数
            function dragstarted(d) {
                if (!d3.event.active) simulation.alphaTarget(0.3).restart();
                d.fx = d.x;
                d.fy = d.y;
            }
            
            function dragged(d) {
                d.fx = d3.event.x;
                d.fy = d3.event.y;
            }
            
            function dragended(d) {
                if (!d3.event.active) simulation.alphaTarget(0);
                d.fx = null;
                d.fy = null;
            }
            </script>
        </body>
        </html>
        ''')
    
    logging.info(f"HTML报告已生成: {html_file}")

def main():
    args = parse_args()
    
    try:
        # 创建输出目录
        os.makedirs(args.output_dir, exist_ok=True)
        
        # 解析输入文件
        logging.info("解析GFA文件...")
        nodes, links = parse_gfa(args.gfa)
        
        logging.info("解析FASTA文件...")
        contigs = parse_fasta(args.fasta)
        
        # 对长读长进行采样
        sampled_reads = os.path.join(args.output_dir, "sampled_reads.fa")
        logging.info(f"对长读长进行采样 ({args.sampling_rate*100:.1f}%)...")
        subsample_reads(args.long_reads, sampled_reads, args.sampling_rate)
        
        # 验证连接
        logging.info("验证所有连接...")
        results = validate_connections(
            nodes, links, contigs, sampled_reads, 
            args.output_dir, args.threads, args.overlap
        )
        
        logging.info(f"验证完成! 结果已写入 {args.output_dir}")
        
        # 生成HTML报告
        generate_html_report(results, args.output_dir, args.gfa)
        
    except Exception as e:
        logging.error(f"错误: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()