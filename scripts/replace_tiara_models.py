#!/usr/bin/env python3
"""
Replace Tiara model files in conda environment
This script backups original models and replaces them with optimized ones
"""

import os
import sys
import shutil
import subprocess
import argparse
from datetime import datetime

def find_tiara_models_dir(conda_env=None):
    """
    Find Tiara models directory
    
    Args:
        conda_env: conda environment name or path, if None, use current environment
    
    Returns:
        Path to Tiara models directory
    """
    # If conda environment is specified, activate it
    cmd = []
    if conda_env:
        # Determine conda path
        conda_path = os.environ.get('CONDA_EXE', 'conda')
        if os.path.exists(conda_env):
            # If conda_env is a path
            cmd = [f"source $(dirname $(dirname {conda_path}))/etc/profile.d/conda.sh && conda activate {conda_env} &&"]
        else:
            # If conda_env is an environment name
            cmd = [f"source $(dirname $(dirname {conda_path}))/etc/profile.d/conda.sh && conda activate {conda_env} &&"]
    
    # Run Python command to find Tiara models directory
    python_cmd = 'python -c "import pkg_resources; print(pkg_resources.resource_filename(\'tiara\', \'models\'))"'
    full_cmd = " ".join(cmd + [python_cmd])
    
    try:
        if conda_env:
            # Use shell=True to ensure conda activate works properly
            result = subprocess.run(full_cmd, shell=True, capture_output=True, text=True)
        else:
            # Run Python command directly
            result = subprocess.run(python_cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            raise RuntimeError(f"Failed to find Tiara models directory: {result.stderr}")
        
        models_dir = result.stdout.strip()
        if not os.path.exists(models_dir):
            raise RuntimeError(f"Found models directory does not exist: {models_dir}")
        
        return models_dir
    
    except Exception as e:
        print(f"Error: {e}")
        print("Unable to automatically find Tiara models directory. Please specify the --tiara-dir parameter manually.")
        sys.exit(1)

def check_permissions(directory):
    """
    Check if we have permission to modify the specified directory
    
    Args:
        directory: Directory to check
    
    Returns:
        Boolean indicating whether we have write permission
    """
    if not os.access(directory, os.W_OK):
        return False
    
    # Test creating a temporary file
    try:
        test_file = os.path.join(directory, ".permission_test")
        with open(test_file, 'w') as f:
            f.write("test")
        os.remove(test_file)
        return True
    except:
        return False

def replace_models(source_dir, tiara_models_dir, backup=True):
    """
    Replace Tiara model files
    
    Args:
        source_dir: Source directory containing optimized models
        tiara_models_dir: Tiara models directory
        backup: Whether to backup original models
    """
    if not os.path.exists(source_dir):
        raise ValueError(f"Source directory does not exist: {source_dir}")
    
    if not os.path.exists(tiara_models_dir):
        raise ValueError(f"Tiara models directory does not exist: {tiara_models_dir}")
    
    # Ensure target directory structure exists
    tfidf_target = os.path.join(tiara_models_dir, "tfidf-models")
    nnet_target = os.path.join(tiara_models_dir, "nnet-models")
    
    if not os.path.exists(tfidf_target) or not os.path.exists(nnet_target):
        raise ValueError(f"Invalid Tiara models directory structure: {tiara_models_dir}")
    
    # Check permissions
    if not check_permissions(tiara_models_dir):
        raise PermissionError(f"Insufficient permissions to modify Tiara models directory: {tiara_models_dir}")
    
    # Backup original models
    backup_dir = None
    if backup:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_dir = os.path.join(os.path.dirname(tiara_models_dir), f"tiara_models_backup_{timestamp}")
        print(f"Backing up original models to: {backup_dir}")
        shutil.copytree(tiara_models_dir, backup_dir)
    
    try:
        # Replace TF-IDF models
        print("Replacing TF-IDF models...")
        for k in [4, 5, 6]:
            src_dir = os.path.join(source_dir, "tfidf", f"k{k}-first-stage")
            dst_dir = os.path.join(tfidf_target, f"k{k}-first-stage")
            
            if os.path.exists(src_dir):
                if os.path.exists(dst_dir):
                    shutil.rmtree(dst_dir)
                print(f"  Copying first stage k={k}: {src_dir} -> {dst_dir}")
                shutil.copytree(src_dir, dst_dir)
        
        for k in [4, 5, 6, 7]:
            src_dir = os.path.join(source_dir, "tfidf", f"k{k}-second-stage")
            dst_dir = os.path.join(tfidf_target, f"k{k}-second-stage")
            
            if os.path.exists(src_dir):
                if os.path.exists(dst_dir) and k == 7:
                    # For k=7, might need to create new directory
                    shutil.rmtree(dst_dir)
                elif os.path.exists(dst_dir):
                    shutil.rmtree(dst_dir)
                
                print(f"  Copying second stage k={k}: {src_dir} -> {dst_dir}")
                shutil.copytree(src_dir, dst_dir)
        
        # Replace neural network models
        print("Replacing neural network models...")
        for k in [4, 5, 6]:
            src_file = os.path.join(source_dir, "nnet", f"first_nnet_kmer_{k}.pkl")
            dst_file = os.path.join(nnet_target, f"first_nnet_kmer_{k}.pkl")
            
            if os.path.exists(src_file):
                print(f"  Copying first stage k={k}: {src_file} -> {dst_file}")
                shutil.copy2(src_file, dst_file)
        
        for k in [4, 5, 6, 7]:
            src_file = os.path.join(source_dir, "nnet", f"second_nnet_kmer_{k}.pkl")
            dst_file = os.path.join(nnet_target, f"second_nnet_kmer_{k}.pkl")
            
            if os.path.exists(src_file):
                print(f"  Copying second stage k={k}: {src_file} -> {dst_file}")
                shutil.copy2(src_file, dst_file)
        
        print(f"Model replacement complete! Tiara will now use optimized models.")
        if backup:
            print(f"To restore original models, run: python {__file__} --restore {backup_dir} --tiara-dir {tiara_models_dir}")
    except Exception as e:
        print(f"Error during replacement: {e}")
        if backup:
            print("Rolling back to original models...")
            restore_backup(backup_dir, tiara_models_dir)
        raise

def restore_backup(backup_dir, tiara_models_dir):
    """
    Restore backed up model files
    
    Args:
        backup_dir: Backup directory
        tiara_models_dir: Tiara models directory
    """
    if not os.path.exists(backup_dir):
        raise ValueError(f"Backup directory does not exist: {backup_dir}")
    
    if not os.path.exists(tiara_models_dir):
        raise ValueError(f"Tiara models directory does not exist: {tiara_models_dir}")
    
    print(f"Restoring backed up model files: {backup_dir} -> {tiara_models_dir}")
    
    # Clear current models directory
    for item in os.listdir(tiara_models_dir):
        path = os.path.join(tiara_models_dir, item)
        if os.path.isdir(path):
            shutil.rmtree(path)
        else:
            os.remove(path)
    
    # Copy backed up files
    for item in os.listdir(backup_dir):
        src = os.path.join(backup_dir, item)
        dst = os.path.join(tiara_models_dir, item)
        
        if os.path.isdir(src):
            shutil.copytree(src, dst)
        else:
            shutil.copy2(src, dst)
    
    print("Model restoration complete! Tiara will now use original models.")

def main():
    parser = argparse.ArgumentParser(description="Replace or restore Tiara model files")
    
    action_group = parser.add_mutually_exclusive_group(required=True)
    action_group.add_argument("--replace", metavar="SOURCE_DIR", help="Replace Tiara models with optimized models from specified directory")
    action_group.add_argument("--restore", metavar="BACKUP_DIR", help="Restore Tiara models from backup directory")
    
    # Get default values from environment variables
    default_conda_env = os.environ.get('TIARA_CONDA_ENV')
    default_tiara_dir = os.environ.get('TIARA_MODELS_DIR')
    
    parser.add_argument("--conda-env", default=default_conda_env, help="Specify conda environment name or path")
    parser.add_argument("--tiara-dir", default=default_tiara_dir, help="Specify Tiara models directory path")
    parser.add_argument("--no-backup", action="store_true", help="Do not create backup when replacing")
    parser.add_argument("--log", help="Log file path")
    
    args = parser.parse_args()
    
    # Set up logging
    if args.log:
        log_dir = os.path.dirname(args.log)
        if log_dir and not os.path.exists(log_dir):
            os.makedirs(log_dir)
        
        # Redirect standard output and error
        log_file = open(args.log, 'w')
        sys.stdout = log_file
        sys.stderr = log_file
    
    # Determine Tiara models directory
    tiara_models_dir = args.tiara_dir
    if tiara_models_dir is None:
        tiara_models_dir = find_tiara_models_dir(args.conda_env)
    
    print(f"Tiara models directory: {tiara_models_dir}")
    
    if args.replace:
        replace_models(args.replace, tiara_models_dir, not args.no_backup)
    elif args.restore:
        restore_backup(args.restore, tiara_models_dir)
    
    # Close log file
    if args.log:
        log_file.close()
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__

if __name__ == "__main__":
    main() 