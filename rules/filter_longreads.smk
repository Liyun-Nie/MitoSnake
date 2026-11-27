rule filter_longreads:
    """
    Filter long reads based on mapping quality and coverage.
    """
    input:
        paf = "{sample}/assembly/mapped_reads.paf",
        reads = "{sample}/subsample/long_subsample_" + str(config["subsample"]["long_cov"]) + "x.fa"
    output:
        filtered = "{sample}/assembly/filtered_longreads.fa"
    params:
        min_match = config["filter"]["min_match"],
        min_coverage = config["filter"]["min_coverage"]
    log:
        "{sample}/logs/filter_longreads.log"
    benchmark:
        "{sample}/logs/filter_longreads.benchmark.txt"
    conda:
       get_conda_env("mitosnake")
    shell:
        """
        # Check if input files exist and are not empty
        for file in {input.paf} {input.reads}; do
            if [ ! -s "$file" ]; then
                echo "ERROR: Input file $file does not exist or is empty. Cannot proceed." | tee -a {log} >&2
                exit 1
            fi
        done
        
        # Check if previous output exists but is empty, delete if so
        if [ -f "{output.filtered}" ] && [ ! -s "{output.filtered}" ]; then
            echo "INFO: Removing empty file: {output.filtered}" >> {log} 2>&1
            rm -f "{output.filtered}"
        fi

        echo "INFO: Starting long read filtering..." >> {log} 2>&1
        python {SCRIPT_DIR}/scripts/filter_longreads.py \
        --paf {input.paf} \
        --reads {input.reads} \
        --output {output.filtered} \
        --min-match {params.min_match} \
        --min-coverage {params.min_coverage} \
        > {log} 2>&1
        
        # Verify output file is not empty
        if [ ! -s {output.filtered} ]; then
            echo "WARNING: Output file {output.filtered} is empty. No reads passed filtering criteria." | tee -a {log} >&2
        else
            echo "SUCCESS: Long read filtering completed successfully." >> {log} 2>&1
        fi
        """

