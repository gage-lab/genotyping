#!/usr/bin/env bash
# Author: Mike Cuoco
# Created on: 6/11/25, 3:09 PM
#
# impute - Bash wrapper for genotype imputation snakemake workflow
# Usage: ./impute [OPTIONS]

set -euo pipefail

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_DIR="${SCRIPT_DIR}/workflow"
SNAKEFILE="${WORKFLOW_DIR}/impute.smk"

# Source common functions
source "${SCRIPT_DIR}/lib/common.sh"

# Show help
show_help() {
    cat << EOF
impute - Run genotype imputation workflow with Michigan Imputation Server. 
         Uses hg38 reference genome and 1000 Genomes Phase 3 WGS reference panel.

USAGE:
    bash impute.sh [OPTIONS]

OPTIONS:
    -i, --input FILE    VCF file to impute (can be used multiple times, required)
    -o, --output DIR    Output directory (required)
    -t, --token TOKEN   Imputation Server API token (required)
    -c, --cores NUM     Number of cores to use (default: 1)
    -n, --dry-run       Show what would be executed without running
    -h, --help          Show this help message

EXAMPLES:
    # Dry run to see what would be executed
    impute -i file1.vcf.gz -i file2.vcf.gz -o results -t your_api_token -n

    # Basic run with multiple files
    impute -i file1.vcf.gz -i file2.vcf.gz -o results -t your_api_token

    # Run with 4 cores
    impute -i file1.vcf.gz -i file2.vcf.gz -o results -t your_api_token -c 4

EOF
}





# Script-specific argument parsing
parse_args() {
    vcf_files=()
    output_dir=""
    token=""
    cores="1"  # Default to 1 core
    dry_run="false"
    
    # Parse command arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -i|--input)
                vcf_files+=("$2")
                shift 2
                ;;
            -o|--output)
                output_dir="$2"
                shift 2
                ;;
            -t|--token)
                token="$2"
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
    if [[ ${#vcf_files[@]} -eq 0 ]]; then
        print_error "VCF files are required. Use -i/--input"
        exit 1
    fi
    
    if [[ -z "$output_dir" ]]; then
        print_error "Output directory is required. Use -o/--output"
        exit 1
    fi
    
    if [[ -z "$token" ]]; then
        print_error "ImputationBot API token is required. Use -t/--token"
        exit 1
    fi
    
    # Validate VCF files exist
    for vcf_file in "${vcf_files[@]}"; do
        if [[ ! -f "$vcf_file" ]]; then
            print_error "VCF file does not exist: $vcf_file"
            exit 1
        fi
    done
    
    # Validate cores is a positive integer
    if ! [[ "$cores" =~ ^[1-9][0-9]*$ ]]; then
        print_error "Cores must be a positive integer. Got: $cores"
        exit 1
    fi
    
    print_info "Using $cores core(s)"
    print_info "VCF files: ${vcf_files[*]}"
}

# Script-specific config building
build_config() {
    # Join array elements with spaces for snakemake config
    local vcf_list="${vcf_files[*]}"
    echo "$cores" "\"vcfs=\\\"$vcf_list\\\"\" \"token=\\\"$token\\\"\""
}

# Main command dispatcher
main() {
    run_main "imputation" "$SCRIPT_DIR" "$SNAKEFILE" "show_help" "parse_args" "validate_args" "build_config" "$@"
}

# Run main function
main "$@" 