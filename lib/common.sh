#!/usr/bin/env bash
# Common functions for workflow scripts
# Author: Mike Cuoco

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Print colored output
print_info() { echo -e "${BLUE}[INFO]${NC} $1"; }
print_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
print_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
print_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Validate required files
validate_setup() {
    local snakefile="$1"
    local issues=0
    
    if [[ ! -f "$snakefile" ]]; then
        print_error "Snakefile not found at: $snakefile"
        ((issues++))
    fi
    
    if ! command -v snakemake &> /dev/null; then
        print_error "snakemake command not found. Please install snakemake."
        ((issues++))
    fi
    
    if [[ $issues -gt 0 ]]; then
        print_error "Setup validation failed. Please fix the above issues."
        exit 1
    fi
}

# Build snakemake command with custom config parameters
build_snakemake_cmd() {
    local snakefile="$1"
    local output_dir="$2" 
    local dry_run="$3"
    local cores="$4"
    shift 4  # Remove first 4 arguments, rest are config key=value pairs
    
    local cmd=("snakemake -c${cores} all --snakefile $snakefile --rerun-incomplete")

    # Pass config values directly via command line
    cmd+=("--config")
    cmd+=("outdir=$(realpath "$output_dir")")
    
    # Add any additional config parameters
    for config_param in "$@"; do
        cmd+=("$config_param")
    done
    
    if [[ "${dry_run}" == "true" ]]; then
        cmd+=("--dry-run")
    fi
    
    echo "${cmd[@]}"
}

# Generic workflow runner
run_workflow() {
    local script_name="$1"
    local snakefile="$2" 
    local parse_args_func="$3"
    local validate_args_func="$4"
    local build_config_func="$5"
    shift 5
    
    # Parse arguments using script-specific function
    eval "$parse_args_func \"\$@\""
    
    # Validate arguments using script-specific function
    eval "$validate_args_func"
    
    # Create output directory
    mkdir -p "$output_dir"
    
    # Build config parameters using script-specific function
    local config_result
    config_result=$(eval "$build_config_func")
    
    # Extract cores (first item) and config params (rest)
    local cores
    cores=$(echo "$config_result" | cut -d' ' -f1)
    local config_params
    config_params=$(echo "$config_result" | cut -d' ' -f2-)
    
    # Build and run snakemake command
    local cmd
    cmd=$(build_snakemake_cmd "$snakefile" "$output_dir" "$dry_run" "$cores" $config_params)
    
    print_info "Running $script_name workflow..."
    print_info "Output directory: $output_dir"
    print_info "Using $cores core(s)"
    print_info "Command: $cmd"
    
    eval "$cmd"
}

# Activate conda environment
activate_conda_env() {
    local script_dir="$1"
    local conda_env_path="${script_dir}/.conda"

    # attempt to switch proper environment
    if ! command -v 'snakemake' &>/dev/null; then
        print_info "Attempting to switch to $conda_env_path environment"
        eval "$(conda shell.bash hook)"
        conda activate $conda_env_path || {
            print_error "Failed to activate conda environment: $conda_env_path"
            exit 1
        }
        print_success "Successfully activated conda environment"
    fi
} 

# Generic main function
run_main() {
    local script_name="$1"
    local script_dir="$2"
    local snakefile="$3"
    local show_help_func="$4"
    local parse_args_func="$5"
    local validate_args_func="$6"
    local build_config_func="$7"
    shift 7
    
    activate_conda_env "$script_dir"
    
    # Validate setup first
    validate_setup "$snakefile"
    
    if [[ $# -eq 0 ]]; then
        eval "$show_help_func"
        exit 1
    fi
    
    # Run the workflow with all arguments
    run_workflow "$script_name" "$snakefile" "$parse_args_func" "$validate_args_func" "$build_config_func" "$@"
}