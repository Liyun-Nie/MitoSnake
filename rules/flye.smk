rule flye:
    """
    Assemble filtered mitochondrial long reads using Flye.
    """
    input:
        reads = "{sample}/assembly/mitochondrion_filtered_longreads.fa"
    output:
        assembly = "{sample}/assembly/flye_out/assembly.fasta",
        gfa = "{sample}/assembly/flye_out/assembly_graph.gfa"
    params:
        outdir = "{sample}/assembly/flye_out",
        
        mode = lambda wildcards: (
            "--pacbio-hifi" if samples_info[wildcards.sample]['longread_type'] == 'hifi' else
            "--pacbio-raw" if samples_info[wildcards.sample]['longread_type'] == 'clr' else
            "--nano-raw"
        )
    log:
        "{sample}/logs/flye.log"
    benchmark:
        "{sample}/logs/flye.benchmark.txt"
    threads: get_resources("flye")["threads"]
    conda:
        get_conda_env("mitosnake")
    shell:
        """
        # Check if input file exists and is not empty
        if [ ! -s "{input.reads}" ]; then
            echo "ERROR: Input file {input.reads} does not exist or is empty. Cannot proceed with assembly." | tee -a {log} >&2
            exit 1
        fi
        
        # Check if previous output exists but is empty, delete if so
        for file in {output.assembly} {output.gfa}; do
            if [ -f "$file" ] && [ ! -s "$file" ]; then
                echo "INFO: Removing empty file: $file" >> {log} 2>&1
                rm -f "$file"
            fi
        done

        # Also clear the output directory if it exists but contains no valid output
        if [ -d "{params.outdir}" ]; then
            if [ ! -s "{params.outdir}/assembly.fasta" ]; then
                echo "INFO: Removing potentially corrupt output directory: {params.outdir}" >> {log} 2>&1
                rm -rf "{params.outdir}"
            fi
        fi
        
        echo "INFO: Starting Flye assembly with mode: {params.mode}..." >> {log} 2>&1
        
        flye {params.mode} \
        {input.reads} \
        --genome-size 1m \
        --threads {threads} \
        --out-dir {params.outdir} \
        2> >(tee -a {log} >&2)
        
        # Verify output files exist and are not empty
        if [ ! -s {output.assembly} ]; then
            echo "ERROR: Output file {output.assembly} is empty. Assembly has failed." | tee -a {log} >&2
            exit 1
        else
            echo "SUCCESS: Flye assembly completed successfully." >> {log} 2>&1
        fi
        """