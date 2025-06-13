#!/usr/bin/env bash
# Author: Mike Cuoco
# Created on: 6/11/25, 3:09 PM
#
# idat2vcf - Bash wrapper for Illumina IDAT to VCF conversion snakemake workflow
# Usage: ./idat2vcf [COMMAND] [OPTIONS]

set -euo pipefail

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_DIR="${SCRIPT_DIR}/workflow"
SNAKEFILE="${WORKFLOW_DIR}/idat2vcf.smk"

# Source common functions
source "${SCRIPT_DIR}/lib/common.sh"

# Show help
show_help() {
    cat << EOF
idat2vcf - Convert Illumina IDAT files to VCF format
           Only supports GSA-24v3-0-A2 and InfinitumCoreExome-24v1-4_A1 arrays.

USAGE:
    bash idat2vcf.sh [OPTIONS]

OPTIONS:
    -i, --input DIR     Input directory containing data of single UCSD IGM SNP Array (required)
    -o, --output DIR    Output directory to store results (required)
    -c, --cores NUM     Number of cores to use (default: 1)
    -n, --dry-run       Show what would be executed without running
    -h, --help          Show this help message

EXAMPLES:
    # Dry run to see what would be executed
    idat2vcf -i /path/to/idat/files -o results -n

    # Basic run
    idat2vcf -i /path/to/idat/files -o results

    # Use 4 cores
    idat2vcf -i /path/to/idat/files -o results -c 4

EOF
}

# Script-specific argument parsing
parse_args() {
    input_dir=""
    output_dir=""
    cores="1"  # Default to 1 core
    dry_run="false"
    
    # Parse command arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -i|--input)
                input_dir="$2"
                shift 2
                ;;
            -o|--output)
                output_dir="$2"
                shift 2
                ;;
            -c|--cores)
                cores="$2"
                shift 2
                ;;
            -n|--dry-run)
                dry_run="true"
                shift
                ;;
            -h|--help)
                show_help
                exit 0
                ;;
            *)
                print_error "Unknown option: $1"
                exit 1
                ;;
        esac
    done
}

# Script-specific argument validation
validate_args() {
    if [[ -z "$input_dir" ]]; then
        print_error "Input directory is required. Use -i/--input"
        exit 1
    fi
    
    if [[ -z "$output_dir" ]]; then
        print_error "Output directory is required. Use -o/--output"
        exit 1
    fi
    
    if [[ ! -d "$input_dir" ]]; then
        print_error "Input directory does not exist: $input_dir"
        exit 1
    fi
    
    # Validate cores is a positive integer
    if ! [[ "$cores" =~ ^[1-9][0-9]*$ ]]; then
        print_error "Cores must be a positive integer. Got: $cores"
        exit 1
    fi
    
    print_info "Using $cores core(s)"
    print_info "Input directory: $input_dir"
}

# Script-specific config building
build_config() {
    echo "$cores" "\"igm_dir=$(realpath "$input_dir")\""
}

# Main command dispatcher
main() {
    run_main "idat2vcf" "$SCRIPT_DIR" "$SNAKEFILE" "show_help" "parse_args" "validate_args" "build_config" "$@"
}

# Run main function
main "$@"