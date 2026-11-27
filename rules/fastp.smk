rule fastp:
    """
    Rule to perform quality control using fastp.
    """
    input:
        fq1 = lambda wildcards: samples_info[wildcards.sample]['read1'],
        fq2 = lambda wildcards: samples_info[wildcards.sample]['read2']
    output:
        fq1 = "{sample}/subsample/R1_clean.fq.gz",
        fq2 = "{sample}/subsample/R2_clean.fq.gz",
        html = "{sample}/logs/fastp.html",
        json = "{sample}/logs/fastp.json"
    log:
        "{sample}/logs/fastp.log"
    benchmark:
        "{sample}/logs/fastp.benchmark.txt"
    threads: get_resources("fastp")["threads"]
    resources:
        mem_mb = get_resources("fastp")["mem_mb"],
        disk_mb = get_resources("fastp")["disk_mb"]
    params:
        cut_front = config["fastp"]["cut_front"],
        cut_front_window_size = config["fastp"]["cut_front_window_size"],
        cut_front_mean_quality = config["fastp"]["cut_front_mean_quality"],
        average_qual = config["fastp"]["average_qual"],
        length_required = config["fastp"]["length_required"]
    conda:
        get_conda_env("mitosnake")
    shell:
        """
        # Check if input files exist and are not empty
        for file in {input.fq1} {input.fq2}; do
            if [ ! -s "$file" ]; then
                echo "ERROR: Input file $file does not exist or is empty. Cannot proceed." | tee -a {log} >&2
                exit 1
            fi
        done
        
        # Check if previous output files exist but are empty, delete if so
        for file in {output.fq1} {output.fq2} {output.html} {output.json}; do
            if [ -f "$file" ] && [ ! -s "$file" ]; then
                echo "INFO: Removing empty file: $file" >> {log} 2>&1
                rm -f "$file"
            fi
        done
        
        echo "INFO: Starting fastp quality control..." >> {log} 2>&1
        set -euo pipefail
        fastp \
        -i {input.fq1} \
        -I {input.fq2} \
        -o {output.fq1} \
        -O {output.fq2} \
        --cut_front {params.cut_front} \
        --cut_front_window_size {params.cut_front_window_size} \
        --cut_front_mean_quality {params.cut_front_mean_quality} \
        --average_qual {params.average_qual} \
        --length_required {params.length_required} \
        --html {output.html} \
        --json {output.json} \
        --thread {threads} \
        2>> {log}
        
        # Verify output files are not empty
        for file in {output.fq1} {output.fq2} {output.html} {output.json}; do
            if [ ! -s "$file" ]; then
                echo "ERROR: Output file $file is empty after processing." | tee -a {log} >&2
                exit 1
            fi
        done
        
        echo "SUCCESS: Fastp quality control completed." >> {log} 2>&1
        """

