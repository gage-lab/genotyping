# Snakemake workflows for SNP array data processing

## Setup


```bash
# Install conda environment:
conda env create -f conda.yaml -p ./.conda

# Download resources
bash download_resources.sh
```

## Workflow 1: IDAT to VCF

Converts IDAT files to VCF format, if SNP array is based on hg19, lifts over to hg38.

```bash
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
```

## Workflow 2: SNP imputation 

Merges VCF files from multiple sample sets, submits to Michigan Imputation Server for imputation, waits for results and saves single vcf file with imputed genotypes. Requires API token from Michigan Imputation Server. See [this tutorial](https://genepi.github.io/michigan-imputationserver/tutorials/api/).

```bash
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
```

## TODO

- [ ] add script to fetch resources
- [ ] set BCFTOOLS tempdir
- [ ] add TopMed functionality
- [ ] add QC checks
   - [ ] check for sample duplicates
   - [ ] report callrates in idat2vcf
   - [ ] report imputation rates in impute

