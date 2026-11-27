rule consolidated_log:
    """
    Consolidate all log files into a single chronological log file.
    """
    input:
        # Get all relevant log files based on mode
        lambda wildcards: (
            # lt-seed mode already has iterative_assembly.log, no need for consolidation
            [
                f"{wildcards.sample}/logs/iterative_assembly.log"
            ] if config.get("module_to_run") == "lt-seed" else
            
            # All log files for ht mode, arranged in execution order
            [
                f"{wildcards.sample}/logs/fastp.log",
                f"{wildcards.sample}/logs/seqkit_stats.log", 
                f"{wildcards.sample}/logs/shortreads_subsample.log",
                f"{wildcards.sample}/logs/spades.log", 
                f"{wildcards.sample}/logs/tiara.log",
                f"{wildcards.sample}/logs/longreads_subsample.log",
                f"{wildcards.sample}/logs/minimap2.log",
                f"{wildcards.sample}/logs/filter_longreads.log", 
                f"{wildcards.sample}/logs/tiara_longreads.log",
                f"{wildcards.sample}/logs/flye.log",
                f"{wildcards.sample}/logs/finalize.log"
            ] if config.get("module_to_run") == "ht" else
            
            # All log files for st mode, arranged in execution order
            [
                f"{wildcards.sample}/logs/fastp.log",
                f"{wildcards.sample}/logs/seqkit_stats.log", 
                f"{wildcards.sample}/logs/shortreads_subsample.log",
                f"{wildcards.sample}/logs/spades.log", 
                f"{wildcards.sample}/logs/tiara.log"
            ]
        )
    output:
        log = "{sample}/logs/consolidated_workflow.log"
    run:
        import os
        from datetime import datetime
        
        # 整合所有日志文件
        with open(output.log, "w") as outf:
            outf.write("============================================\n")
            outf.write(f"MITO-SNAKE Pipeline - {config['module_to_run'].upper()} Mode\n")
            outf.write(f"Sample: {wildcards.sample}\n")
            outf.write(f"Started at: {datetime.now().strftime('%a %b %d %H:%M:%S %Z %Y')}\n")
            outf.write("============================================\n\n")
            
            # 添加各个模块的日志
            for log_file in input:
                if os.path.exists(log_file) and os.path.getsize(log_file) > 0:
                    # 获取步骤名称
                    step_name = os.path.basename(log_file).replace('.log', '').replace('_', ' ').upper()
                    
                    # 添加步骤标题
                    outf.write("\n--------------------------------------------\n")
                    outf.write(f"STEP: {step_name}\n")
                    outf.write("--------------------------------------------\n")
                    
                    # 添加日志内容
                    with open(log_file, "r") as inf:
                        outf.write(inf.read())
                    outf.write("\n")
            
            # 添加完成标记
            outf.write("\n============================================\n")
            outf.write(f"Pipeline completed at: {datetime.now().strftime('%a %b %d %H:%M:%S %Z %Y')}\n")
            outf.write("============================================\n") 