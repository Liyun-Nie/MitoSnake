rule spades:
    """
    Rule to assemble short reads using SPAdes.
    """
    input:
        r1 = "{sample}/subsample/R1_subsample_" + str(config["subsample"]["short_cov"]) + "x.fq.gz",
        r2 = "{sample}/subsample/R2_subsample_" + str(config["subsample"]["short_cov"]) + "x.fq.gz"
    output:
        contigs = "{sample}/assembly/contigs.fasta",
        scaffolds = "{sample}/assembly/scaffolds.fasta"
    log:
        "{sample}/logs/spades.log"
    benchmark:
        "{sample}/logs/spades.benchmark.txt"
    params:
        outdir = "{sample}/assembly/spades_out",
        memory = lambda wildcards, resources: resources.mem_mb // 1024
    threads: get_resources("spades")["threads"]
    resources:
        mem_mb = get_resources("spades")["mem_mb"],
        disk_mb = get_resources("spades")["disk_mb"]
    conda:
        get_conda_env("mitosnake")
    shell:
        """
        # Check if input files exist and are not empty
        for file in {input.r1} {input.r2}; do
            if [ ! -s "$file" ]; then
                echo "ERROR: Input file $file does not exist or is empty. Cannot proceed." | tee -a {log} >&2
                exit 1
            fi
        done
        
        # Check if previous output exists but is empty, delete if so
        for file in {output.contigs} {output.scaffolds}; do
            if [ -f "$file" ] && [ ! -s "$file" ]; then
                echo "INFO: Removing empty file: $file" >> {log} 2>&1
                rm -f "$file"
            fi
        done

        echo "INFO: Starting SPAdes assembly..." >> {log} 2>&1
        spades.py \
        -1 {input.r1} \
        -2 {input.r2} \
        -o {params.outdir} \
        -t {threads} \
        -m {params.memory} \
        --careful \
        --phred-offset 33 \
        > {log} 2>&1
        
        echo "INFO: Copying SPAdes output files..." >> {log} 2>&1
        cp {params.outdir}/contigs.fasta {output.contigs}
        cp {params.outdir}/scaffolds.fasta {output.scaffolds}
        
        # Verify output files are not empty
        for file in {output.contigs} {output.scaffolds}; do
            if [ ! -s "$file" ]; then
                echo "ERROR: Output file $file is empty after assembly." | tee -a {log} >&2
                exit 1
            fi
        done
        
        echo "SUCCESS: SPAdes assembly completed." >> {log} 2>&1
        """
