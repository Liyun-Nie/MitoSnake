# Add get_conda_env function
def get_conda_env(env_name):
    """Get the conda environment name.
    
    Args:
        env_name: Name of the conda environment.
    
    Returns:
        Name of the conda environment.
    """
    # Always return environment name, not path
    return env_name

# Get run mode
run_mode = config.get("module_to_run", "ht")

# Check if custom mitochondrial seed exists
has_mito_seed = "mito_seed" in config and config["mito_seed"]

# Check if deep learning filtering is enabled
use_dl_filter = config.get("use_dl_filter", True)

# Load rules common to all modes
include: "consolidated_log.smk"  # Add new log consolidation rule

# Load different rules based on mode
if run_mode == "st":
    # Short-read assembly mode
    include: "fastp.smk"
    include: "seqkit.smk"
    include: "spades.smk"
    
    # Use tiara or bypass based on user choice
    if use_dl_filter:
        include: "tiara.smk"
    else:
        include: "bypass_tiara.smk"

elif run_mode == "ht":
    # Hybrid assembly mode
    include: "fastp.smk"
    include: "seqkit.smk"
    include: "spades.smk"
    
    # Use tiara or bypass based on user choice
    if use_dl_filter:
        include: "tiara.smk"
    else:
        include: "bypass_tiara.smk"
        
    include: "minimap2.smk"
    include: "filter_longreads.smk"
    
    # Use tiara_longreads or bypass based on user choice
    if use_dl_filter:
        include: "tiara_longreads.smk"
    else:
        include: "bypass_tiara_longreads.smk"
        
    include: "flye.smk"
    include: "finalize.smk"

elif run_mode == "lt-seed":
    # lt-seed mode workflow
    include: "subsample_longreads.smk"  # Step 1: Subsample long reads
    include: "iterative_assembly.smk"   # Use iterative assembly pipeline with all steps intergrated

# Summary rules needed for all modes
include: "summary.smk"

# summary error info 
include: "error_summary.smk"