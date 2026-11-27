rule summarize_resources:
    """
    Summarize resource usage from all rules.
    """
    input:
        # Use conditional statement to determine input files based on run mode
        lambda wildcards: (
            # Benchmark files used in lt-seed mode
            [
                f"{wildcards.sample}/logs/subsample_longreads.benchmark.txt",
                f"{wildcards.sample}/logs/iterative_assembly.benchmark.txt"
            ] if config.get("module_to_run") == "lt-seed" else
            
            # Benchmark files used in ht mode - polish.benchmark.txt removed
            [
                f"{wildcards.sample}/logs/fastp.benchmark.txt",
                f"{wildcards.sample}/logs/seqkit.benchmark.txt", 
                f"{wildcards.sample}/logs/spades.benchmark.txt", 
                f"{wildcards.sample}/logs/tiara.benchmark.txt",
                f"{wildcards.sample}/logs/minimap2.benchmark.txt",
                f"{wildcards.sample}/logs/filter_longreads.benchmark.txt", 
                f"{wildcards.sample}/logs/tiara_longreads.benchmark.txt",
                f"{wildcards.sample}/logs/flye.benchmark.txt",
                f"{wildcards.sample}/logs/finalize.benchmark.txt"
                # Reference to polish.benchmark.txt removed
            ] if config.get("module_to_run") == "ht" else
            
            # Benchmark files used in st mode
            [
                f"{wildcards.sample}/logs/fastp.benchmark.txt",
                f"{wildcards.sample}/logs/seqkit.benchmark.txt", 
                f"{wildcards.sample}/logs/spades.benchmark.txt", 
                f"{wildcards.sample}/logs/tiara.benchmark.txt"
            ]
        )
    output:
        "{sample}/logs/resource_summary.tsv"
    run:
        import pandas as pd
        import os
        
        dfs = []
        for bench_file in input:
            if os.path.exists(bench_file):
                rule_name = os.path.basename(bench_file).replace(".benchmark.txt", "")
                try:
                    df = pd.read_csv(bench_file, sep="\t")
                    df['rule'] = rule_name
                    dfs.append(df)
                except Exception as e:
                    print(f"Warning: Could not read {bench_file}: {e}")
        
        if dfs:
            summary_df = pd.concat(dfs)
            summary = summary_df.groupby('rule').agg({
                's': ['mean', 'max'],
                'max_rss': ['mean', 'max'],
                'max_vms': ['mean', 'max'],
                'io_in': ['sum'],
                'io_out': ['sum'],
                'mean_load': ['mean', 'max']
            }).round(2)
            
            # Add readable time and memory units
            summary.columns = ['_'.join(col) for col in summary.columns]
            summary['runtime_hours'] = summary['s_mean'] / 3600
            summary['max_memory_gb'] = summary['max_rss_max'] / 1024 / 1024
            
            summary.to_csv(output, sep="\t") 
        else:
            # If no benchmark files are found, create an empty summary file
            with open(output, 'w') as f:
                f.write("No benchmark files found\n") 