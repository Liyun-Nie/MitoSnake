rule bypass_tiara_longreads:
    """
    Rule to bypass Tiara classification for long reads.
    This rule is used when deep learning filtering is disabled.
    """
    input:
        reads = "{sample}/assembly/filtered_longreads.fa"
    output:
        mit_reads = "{sample}/assembly/mitochondrion_filtered_longreads.fa",
        out = "{sample}/assembly/tiara_long_out.txt",
        log_out = "{sample}/assembly/tiara_long_log_out.txt"
    log:
        "{sample}/logs/bypass_tiara_longreads.log"
    benchmark:
        "{sample}/logs/bypass_tiara_longreads.benchmark.txt"
    shell:
        """
        # Check if input file exists and is not empty
        if [ ! -s "{input.reads}" ]; then
            echo "ERROR: Input file {input.reads} does not exist or is empty. Cannot proceed." | tee -a {log} >&2
            exit 1
        fi
        
        echo "INFO: Deep learning filtering disabled. Copying all filtered long reads without Tiara classification..." >> {log} 2>&1
        
        # Copy all filtered reads to mitochondrion_filtered_longreads.fa
        cp {input.reads} {output.mit_reads}
        
        # Create dummy tiara output files
        echo "# Deep learning filtering disabled - dummy file" > {output.out}
        echo "# All long reads are treated as mitochondrial" >> {output.out}
        echo "# Bypass Tiara for long reads" > {output.log_out}
        
        echo "SUCCESS: Bypassed Tiara classification for long reads, all filtered reads treated as mitochondrial." >> {log} 2>&1
        """ 