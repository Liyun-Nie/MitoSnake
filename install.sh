#!/bin/bash

# Enable error tracking
set -e

# Script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Check conda installation
if ! command -v conda &> /dev/null; then
    echo "Conda not found. Please install Miniconda3 first."
    exit 1
fi

# Check if mamba is installed, install if not
if ! command -v mamba &> /dev/null; then
    echo "Mamba not found. Installing mamba..."
    conda install -c conda-forge mamba -y
fi

# Get conda base path
CONDA_BASE=$(conda info --base)
echo "Conda base path: $CONDA_BASE"

# More reliable check if environment exists
check_env_exists() {
    local env_name=$1
    if [ -d "$CONDA_BASE/envs/$env_name" ]; then
        return 0  # Environment exists
    else
        return 1  # Environment doesn't exist
    fi
}

# Parse command line arguments
REPLACE_TIARA_MODELS=true
for arg in "$@"; do
    case $arg in
        --no-replace-tiara)
            REPLACE_TIARA_MODELS=false
            shift
            ;;
    esac
done

# Create main environment
echo "Creating main mitosnake environment..."
if ! check_env_exists "mitosnake"; then
    echo "Environment mitosnake does not exist, creating..."
    mamba env create -f envs/mitosnake.yml
else
    echo "Environment mitosnake already exists, updating..."
    mamba env update -f envs/mitosnake.yml
fi

# Create tool-specific environments
echo "Creating tool-specific environments..."
for env_file in envs/mitosnake_*.yml; do
    if [[ -f "$env_file" ]]; then
        env_name=$(basename "$env_file" .yml)
        if ! check_env_exists "$env_name"; then
            echo "Creating environment: $env_name"
            mamba env create -f "$env_file"
        else
            echo "Environment $env_name already exists, updating..."
            mamba env update -f "$env_file"
        fi
    fi
done

echo "All environments created/updated successfully!"

# Check for optimized Tiara models
TIARA_MODELS_DIR="$SCRIPT_DIR/resources/tiara_models"
FLAG_FILE="$SCRIPT_DIR/.tiara_models_replaced"

if [[ "$REPLACE_TIARA_MODELS" = true && -d "$TIARA_MODELS_DIR" && ! -f "$FLAG_FILE" ]]; then
    echo "Found optimized Tiara models, proceeding with replacement..."
    
    # Check if there are sufficient permissions to modify conda environment
    TIARA_ENV_PATH="$CONDA_BASE/envs/mitosnake_tiara"
    if [[ ! -w "$TIARA_ENV_PATH" ]]; then
        echo "Warning: Insufficient permissions to modify Tiara environment ($TIARA_ENV_PATH)"
        echo "Please run this script with administrator privileges, or ensure you have permissions to modify conda environments."
        read -p "Continue attempting to replace models? (y/n) " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "Skipping model replacement."
            exit 0
        fi
    fi
    
    # Activate mitosnake environment to run Python script
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate mitosnake
    
    # Run model replacement script
    REPLACEMENT_LOG="$SCRIPT_DIR/logs/tiara_model_replacement.log"
    mkdir -p "$(dirname "$REPLACEMENT_LOG")"
    echo "Starting Tiara model replacement ($(date '+%Y-%m-%d %H:%M:%S'))" > "$REPLACEMENT_LOG"
    python "$SCRIPT_DIR/scripts/replace_tiara_models.py" --replace "$TIARA_MODELS_DIR" --conda-env "mitosnake_tiara" >> "$REPLACEMENT_LOG" 2>&1
    
    # Check execution result
    if [ $? -eq 0 ]; then
        echo "Tiara model replacement successful!"
        # Create flag file indicating models have been replaced
        echo "Tiara models replaced on $(date '+%Y-%m-%d %H:%M:%S')" > "$FLAG_FILE"
    else
        echo "Warning: Tiara model replacement failed, will use default models."
    fi
    
    # Exit environment
    conda deactivate
else
    if [[ -f "$FLAG_FILE" ]]; then
        echo "Tiara models have already been replaced with optimized versions."
    elif [[ ! -d "$TIARA_MODELS_DIR" ]]; then
        echo "Optimized Tiara models directory not found, will use default models."
    fi
fi

echo "You can now use MITO-SNAKE with: python mitosnake.py --help" 