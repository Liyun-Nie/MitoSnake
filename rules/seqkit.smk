rule seqkit:
    """
    Rule to subsample short reads to desired coverage.
    """
    input:
        r1 = "{sample}/subsample/R1_clean.fq.gz",
        r2 = "{sample}/subsample/R2_clean.fq.gz"
    output:
        r1 = "{sample}/subsample/R1_subsample_" + str(config["subsample"]["short_cov"]) + "x.fq.gz",
        r2 = "{sample}/subsample/R2_subsample_" + str(config["subsample"]["short_cov"]) + "x.fq.gz"
    params:
        genome_size = lambda wildcards: get_genome_size(wildcards),
        coverage = config["subsample"]["short_cov"]
    log:
        stats = "{sample}/logs/seqkit_stats.log",
        subsample = "{sample}/logs/shortreads_subsample.log"
    benchmark:
        "{sample}/logs/seqkit.benchmark.txt"
    threads: get_resources("seqkit")["threads"]
    conda:
        get_conda_env("mitosnake")
    shell:
        """
        # Check if input files exist and are not empty
        for file in {input.r1} {input.r2}; do
            if [ ! -s "$file" ]; then
                echo "ERROR: Input file $file does not exist or is empty. Cannot proceed." | tee -a {log.subsample} >&2
                exit 1
            fi
        done
        
        # Check if previous output exists but is empty, delete if so
        for file in {output.r1} {output.r2}; do
            if [ -f "$file" ] && [ ! -s "$file" ]; then
                echo "INFO: Removing empty file: $file" >> {log.subsample} 2>&1
                rm -f "$file"
            fi
        done

        echo "INFO: Starting short read subsampling..." >> {log.subsample} 2>&1
        python {SCRIPT_DIR}/scripts/seqkit_subsample.py \
        {input.r1} {input.r2} \
        {output.r1} {output.r2} \
        {threads} {params.genome_size} {params.coverage} \
        {log.stats} {log.subsample}
        
        # Verify output files are not empty
        for file in {output.r1} {output.r2}; do
            if [ ! -s "$file" ]; then
                echo "ERROR: Output file $file is empty after subsampling." | tee -a {log.subsample} >&2
                exit 1
            fi
        done
        
        echo "SUCCESS: Short read subsampling completed." >> {log.subsample} 2>&1
        """
