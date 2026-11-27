rule subsample_longreads:
    """
    Rule to subsample long reads to desired coverage before mapping.
    This allows reusing the subsampled reads across the workflow.
    """
    input:
        reads = lambda wildcards: samples_info[wildcards.sample]['long_read']
    output:
        subsampled = "{sample}/subsample/long_subsample_" + str(config["subsample"]["long_cov"]) + "x.fa"
    params:
        genome_size = lambda wildcards: get_genome_size(wildcards),
        coverage = config["subsample"]["long_cov"]
    log:
        "{sample}/logs/longreads_subsample.log"
    benchmark:
        "{sample}/logs/subsample_longreads.benchmark.txt"
    threads: get_resources("minimap2")["threads"]
    conda:
        get_conda_env("mitosnake")
    shell:
        """
        # Check if input file exists and is not empty
        if [ ! -s "{input.reads}" ]; then
            echo "ERROR: Input file {input.reads} does not exist or is empty." | tee -a {log} >&2
            exit 1
        fi
        
        # Check if previous output exists but is empty, delete if so
        if [ -f "{output.subsampled}" ] && [ ! -s "{output.subsampled}" ]; then
            echo "INFO: Removing empty file: {output.subsampled}" >> {log} 2>&1
            rm -f "{output.subsampled}"
        fi

        # Subsample long reads
        echo "INFO: Starting long reads subsampling..." >> {log} 2>&1
        python {SCRIPT_DIR}/scripts/extract_longreads_subsample.py \\
            {input.reads} \\
            {output.subsampled} \\
            {threads} {params.genome_size} {params.coverage} \\
            >> {log} 2>&1
        
        # Verify output file is not empty
        if [ ! -s {output.subsampled} ]; then
            echo "WARNING: Output file {output.subsampled} is empty. No reads passed filtering criteria." | tee -a {log} >&2
        else
            echo "SUCCESS: Long read subsampling completed successfully." >> {log} 2>&1
        fi
        """
