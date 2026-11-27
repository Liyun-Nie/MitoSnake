rule bypass_tiara:
    """
    Rule to bypass Tiara classification and use all contigs.
    This rule is used when deep learning filtering is disabled.
    """
    input:
        contigs = "{sample}/assembly/contigs.fasta"
    output:
        mit_contigs = "{sample}/assembly/mitochondrion_contigs.fasta",
        out = "{sample}/assembly/tiara_out.txt",
        log_out = "{sample}/assembly/tiara_log_out.txt"
    log:
        "{sample}/logs/bypass_tiara.log"
    benchmark:
        "{sample}/logs/bypass_tiara.benchmark.txt"
    shell:
        """
        # Check if input file exists and is not empty
        if [ ! -s "{input.contigs}" ]; then
            echo "ERROR: Input file {input.contigs} does not exist or is empty. Cannot proceed." | tee -a {log} >&2
            exit 1
        fi
        
        echo "INFO: Deep learning filtering disabled. Copying all contigs without Tiara classification..." >> {log} 2>&1
        
        # Copy all contigs to mitochondrion_contigs.fasta
        cp {input.contigs} {output.mit_contigs}
        
        # Create dummy tiara output files
        echo "# Deep learning filtering disabled - dummy file" > {output.out}
        echo "# All contigs are treated as mitochondrial" >> {output.out}
        echo "# Bypass Tiara" > {output.log_out}
        
        echo "SUCCESS: Bypassed Tiara classification, all contigs treated as mitochondrial." >> {log} 2>&1
        """ 