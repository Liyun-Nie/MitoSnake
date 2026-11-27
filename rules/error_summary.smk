rule error_summary:
    output:
        summary = "error_summary.txt"
    shell:
        r"""
        echo "ERROR SUMMARY REPORT" > {output.summary}
        echo "===================" >> {output.summary}
        echo "Generated at: $(date)" >> {output.summary}
        echo "" >> {output.summary}
        
        # Dynamically find existing log files
        log_files=""
        for sample in {SAMPLES}; do
            for rule in fastp seqkit spades tiara minimap2 filter_longreads tiara_longreads flye polish; do
                if [ -f "$sample/logs/$rule.log" ]; then
                    log_files="$log_files $sample/logs/$rule.log"
                fi
            done
        done
        
        # Check for errors and warnings in discovered log files
        for log in $log_files; do
            if [ -f "$log" ]; then
                if grep -q "ERROR\|WARNING" "$log" 2>/dev/null; then
                    echo "Issues found in $log:" >> {output.summary}
                    grep "ERROR\|WARNING" "$log" >> {output.summary}
                    echo "-------------------" >> {output.summary}
                fi
            fi
        done
        
        # Add summary information if no errors were found
        if [ $(grep -c "ERROR\|WARNING\|missing\|empty" {output.summary}) -eq 0 ]; then
            echo "No errors or warnings found in log files." >> {output.summary}
        fi
        
        echo "" >> {output.summary}
        echo "For full logs, check individual log files in each sample's logs directory." >> {output.summary}
        
        # Output to standard output
        cat {output.summary}
        """