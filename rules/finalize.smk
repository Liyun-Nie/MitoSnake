rule finalize:
    """
    Create final output files from Flye assembly for both hybrid and seed-guided mode.
    This rule processes both assembly.fasta and assembly_graph.gfa files.
    """
    input:
        assembly = "{sample}/assembly/flye_out/assembly.fasta",
        gfa = "{sample}/assembly/flye_out/assembly_graph.gfa"
    output:
        final_assembly = "{sample}/assembly/final_assembly.fasta" if config.get("module_to_run") == "lt-seed" else "{sample}/assembly/hybrid_assembly.fasta",
        final_gfa = "{sample}/assembly/final_assembly.gfa" if config.get("module_to_run") == "lt-seed" else "{sample}/assembly/hybrid_assembly.gfa"
    log:
        "{sample}/logs/finalize.log"
    benchmark:
        "{sample}/logs/finalize.benchmark.txt"
    params:
        mode = config.get("module_to_run")
    shell:
        """
        # Check if input files exist and are not empty
        if [ ! -s "{input.assembly}" ]; then
            echo "ERROR: Input file {input.assembly} does not exist or is empty. Cannot proceed." | tee -a {log} >&2
            exit 1
        fi
        
        if [ ! -s "{input.gfa}" ]; then
            echo "WARNING: Input GFA file {input.gfa} does not exist or is empty." >> {log} 2>&1
        fi
        
        # Check if previous outputs exist but are empty, delete if so
        for output_file in {output.final_assembly} {output.final_gfa}; do
            if [ -f "$output_file" ] && [ ! -s "$output_file" ]; then
                echo "INFO: Removing empty file: $output_file" >> {log} 2>&1
                rm -f "$output_file"
            fi
        done
        
        # Copy the assembly files to the final output
        if [ "{params.mode}" == "lt-seed" ]; then
            echo "INFO: Finalizing long-read seed-guided assembly..." >> {log} 2>&1
        else
            echo "INFO: Finalizing hybrid assembly..." >> {log} 2>&1
        fi
        
        # Copy assembly FASTA
        cp {input.assembly} {output.final_assembly} 2>> {log}
        
        # Copy assembly GFA if exists
        if [ -s "{input.gfa}" ]; then
            cp {input.gfa} {output.final_gfa} 2>> {log}
        else
            touch {output.final_gfa}
            echo "WARNING: Empty GFA file created as placeholder." >> {log} 2>&1
        fi
        
        # Verify the copied files are not empty
        if [ ! -s {output.final_assembly} ]; then
            echo "ERROR: Finalized assembly FASTA file is empty!" | tee -a {log} >&2
            exit 1
        fi
        
        if [ "{params.mode}" == "lt-seed" ]; then
            echo "SUCCESS: Long-read seed-guided assembly finalized." >> {log} 2>&1
        else
            echo "SUCCESS: Hybrid assembly finalized." >> {log} 2>&1
        fi
        """