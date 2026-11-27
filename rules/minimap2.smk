rule minimap2:
    """
    Subsample long reads and map to mitochondrial candidate contigs.
    """
    input:
        ref = "{sample}/assembly/mitochondrion_contigs.fasta",
        reads = lambda wildcards: samples_info[wildcards.sample]['long_read']
    output:
        subsampled = "{sample}/subsample/long_subsample_" + str(config["subsample"]["long_cov"]) + "x.fa",
        paf = "{sample}/assembly/mapped_reads.paf"
    log:
        subsample = "{sample}/logs/longreads_subsample.log",
        minimap = "{sample}/logs/minimap2.log"
    benchmark:
        "{sample}/logs/minimap2.benchmark.txt"
    params:
        genome_size = lambda wildcards: get_genome_size(wildcards),
        coverage = config["subsample"]["long_cov"],
        preset = lambda wildcards: "map-pb" if samples_info[wildcards.sample]['longread_type'] in ['clr', 'hifi'] else "map-ont"
    threads: get_resources("minimap2")["threads"]
    conda:
        get_conda_env("mitosnake")
    shell:
        """
        # Check if input files exist and are not empty
        for file in {input.ref} {input.reads}; do
            if [ ! -s "$file" ]; then
                echo "ERROR: Input file $file does not exist or is empty. Cannot proceed." | tee -a {log.subsample} >&2
                exit 1
            fi
        done
        
        # Check if previous output exists but is empty, delete if so
        for file in {output.subsampled} {output.paf}; do
            if [ -f "$file" ] && [ ! -s "$file" ]; then
                echo "INFO: Removing empty file: $file" >> {log.subsample} 2>&1
                rm -f "$file"
            fi
        done

        # Subsample long reads
        echo "INFO: Starting long reads subsampling..." >> {log.subsample} 2>&1
        python {SCRIPT_DIR}/scripts/extract_longreads_subsample.py \
        {input.reads} \
        {output.subsampled} \
        {threads} {params.genome_size} {params.coverage} \
        >> {log.subsample} 2>&1
        
        # Check if subsampled file is empty
        if [ ! -s {output.subsampled} ]; then
            echo "ERROR: Subsampled file is empty. Cannot proceed with mapping." | tee -a {log.minimap} >&2
            exit 1
        fi

        # Map reads and output PAF format
        echo "INFO: Starting minimap2 mapping..." >> {log.minimap} 2>&1
        minimap2 -x {params.preset} \
        -t {threads} \
        {input.ref} \
        {output.subsampled} \
        > {output.paf} 2>> {log.minimap}
        
        # Verify PAF file is not empty
        if [ ! -s {output.paf} ]; then
            echo "WARNING: PAF file is empty. No mappings were found." | tee -a {log.minimap} >&2
        else
            echo "SUCCESS: Minimap2 mapping completed." >> {log.minimap} 2>&1
        fi
        """