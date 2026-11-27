rule tiara_longreads:
    """
    Rule to classify filtered long reads using Tiara.
    """
    input:
        reads = "{sample}/assembly/filtered_longreads.fa"
    output:
        mit_reads = "{sample}/assembly/mitochondrion_filtered_longreads.fa",
        out = "{sample}/assembly/tiara_long_out.txt",
        log_out = "{sample}/assembly/tiara_long_log_out.txt"
    log:
        "{sample}/logs/tiara_longreads.log"
    benchmark:
        "{sample}/logs/tiara_longreads.benchmark.txt"
    threads: get_resources("default")["threads"]
    conda:
        get_conda_env("mitosnake_tiara")
    shell:
        """
        # Check if input file exists and is not empty
        if [ ! -s "{input.reads}" ]; then
            echo "ERROR: Input file {input.reads} does not exist or is empty. Cannot proceed." | tee -a {log} >&2
            exit 1
        fi
        
        # Check if previous output exists but is empty, delete if so
        for file in {output.mit_reads} {output.out} {output.log_out}; do
            if [ -f "$file" ] && [ ! -s "$file" ]; then
                echo "INFO: Removing empty file: $file" >> {log} 2>&1
                rm -f "$file"
            fi
        done

        echo "INFO: Starting Tiara classification for long reads..." >> {log} 2>&1
        # Run Tiara on filtered long reads
        tiara -i {input.reads} \
        -o {output.out} \
        --tf mit pla \
        -t {threads} \
        > {output.log_out} 2> {log}
        
        # Verify output files are not empty
        for file in {output.mit_reads} {output.out} {output.log_out}; do
            if [ ! -s "$file" ]; then
                echo "ERROR: Output file $file is empty after classification." | tee -a {log} >&2
                exit 1
            fi
        done
        
        echo "SUCCESS: Tiara long reads classification completed." >> {log} 2>&1
        """ 