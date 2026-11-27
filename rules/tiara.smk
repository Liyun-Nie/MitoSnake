rule tiara:
    """
    Rule to classify contigs using Tiara.
    """
    input:
        contigs = "{sample}/assembly/contigs.fasta"
    output:
        mit_contigs = "{sample}/assembly/mitochondrion_contigs.fasta",
        out = "{sample}/assembly/tiara_out.txt",
        log_out = "{sample}/assembly/tiara_log_out.txt"
    log:
        "{sample}/logs/tiara.log"
    benchmark:
        "{sample}/logs/tiara.benchmark.txt"
    threads: get_resources("default")["threads"]
    conda:
         get_conda_env("mitosnake_tiara")
    shell:
        """
        # Check if input file exists and is not empty
        if [ ! -s "{input.contigs}" ]; then
            echo "ERROR: Input file {input.contigs} does not exist or is empty. Cannot proceed." | tee -a {log} >&2
            exit 1
        fi
        
        # Check if previous output exists but is empty, delete if so
        for file in {output.mit_contigs} {output.out} {output.log_out}; do
            if [ -f "$file" ] && [ ! -s "$file" ]; then
                echo "INFO: Removing empty file: $file" >> {log} 2>&1
                rm -f "$file"
            fi
        done

        echo "INFO: Starting Tiara classification..." >> {log} 2>&1
        tiara -i {input.contigs} \
        -o {output.out} \
        --tf mit pla \
        -t {threads} \
        > {output.log_out} 2> {log}
        
        # Verify output files are not empty
        for file in {output.mit_contigs} {output.out} {output.log_out}; do
            if [ ! -s "$file" ]; then
                echo "ERROR: Output file $file is empty after classification." | tee -a {log} >&2
                exit 1
            fi
        done
        
        echo "SUCCESS: Tiara classification completed." >> {log} 2>&1
        """
