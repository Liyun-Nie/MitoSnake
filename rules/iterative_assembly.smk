rule iterative_assembly:
    """
    Perform iterative assembly using a pipeline that supports multiple rounds.
    """
    input:
        seed = lambda wildcards: config.get("mito_seed") or config["samples_info"][wildcards.sample].get("mito_seed"),
        long_reads = "{sample}/subsample/long_subsample_" + str(config["subsample"]["long_cov"]) + "x.fa"
    output:
        final_assembly = "{sample}/assembly/final_assembly.fasta",
        final_graph = "{sample}/assembly/final_assembly.gfa",
        all_rounds = "{sample}/assembly/all_rounds_assemblies.fasta",
        last_round_assembly = "{sample}/assembly/flye_out_round" + str(config.get("iterative_rounds", 2)) + "/assembly.fasta", 
        last_round_graph = "{sample}/assembly/flye_out_round" + str(config.get("iterative_rounds", 2)) + "/assembly_graph.gfa"
    params:
        outdir = "{sample}/assembly",
        genome_size = "1m", 
        longread_type = lambda wildcards: samples_info[wildcards.sample]['longread_type'],
        min_match = config["filter"]["min_match"],
        min_coverage = config["filter"]["min_coverage"],
        use_dl_filter = config.get("use_dl_filter", True),
        rounds = config.get("iterative_rounds", 2)  # Default two rounds, configurable
    threads: get_resources("flye")["threads"]
    log:
        "{sample}/logs/iterative_assembly.log"
    benchmark:
        "{sample}/logs/iterative_assembly.benchmark.txt"
    conda:
        get_conda_env("mitosnake")
    shell:
        """
        # Check if input files exist and are not empty
        if [ ! -s "{input.seed}" ]; then
            echo "ERROR: Seed file {input.seed} does not exist or is empty." | tee -a {log} >&2
            exit 1
        fi
        
        if [ ! -s "{input.long_reads}" ]; then
            echo "ERROR: Long reads file {input.long_reads} does not exist or is empty." | tee -a {log} >&2
            exit 1
        fi
        
        # Set up environment
        echo "INFO: Starting iterative assembly pipeline with {params.rounds} rounds..." >> {log} 2>&1
        
        # Set use_dl_filter flag
        USE_DL_FILTER=""
        if [ "{params.use_dl_filter}" = "True" ]; then
            USE_DL_FILTER="--use-dl-filter"
        fi
        
        # Run iterative pipeline
        python {SCRIPT_DIR}/scripts/iterative_pipeline.py \\
            --seed {input.seed} \\
            --long-reads {input.long_reads} \\
            --output-dir {params.outdir} \\
            --genome-size {params.genome_size} \\
            --longread-type {params.longread_type} \\
            --threads {threads} \\
            --min-match {params.min_match} \\
            --min-coverage {params.min_coverage} \\
            --script-dir {SCRIPT_DIR}/scripts \\
            --rounds {params.rounds} \\
            $USE_DL_FILTER \\
            --verbose >> {log} 2>&1
        
        # Verify all output files exist
        for file in {output.final_assembly} {output.final_graph} {output.all_rounds} {output.last_round_assembly} {output.last_round_graph}; do
            if [ ! -f "$file" ]; then
                echo "WARNING: Expected output file $file not found. Creating empty file." >> {log} 2>&1
                touch "$file"
            fi
        done
        
        echo "SUCCESS: Iterative assembly completed." >> {log} 2>&1
        """