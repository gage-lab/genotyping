# Snakemake workflows for SNP array data processing

## Setup

Install conda environment:

```bash
conda env create -f conda.yaml -p ./.conda
```

## Workflow 1: IDAT to VCF

```bash
idat2vcf - Convert Illumina IDAT files to VCF format

USAGE:
    idat2vcf [OPTIONS]

OPTIONS:
    -i, --input DIR     Input directory containing IDAT files (required)
    -o, --output DIR    Output directory (required)
    -c, --cores NUM     Number of cores to use (default: 1)
    --dry-run           Show what would be executed without running
    -h, --help          Show this help message

EXAMPLES:
    # Basic run
    idat2vcf -i /path/to/igm/snp_data -o /path/to/outdir

    # Use 4 cores
    idat2vcf -i /path/to/igm/snp_data -o /path/to/outdir -c 4
```
## Workflow 2: SNP imputation 

```bash
impute - Run genotype imputation workflow

USAGE:
    impute [OPTIONS]

OPTIONS:
    -i, --input FILE    VCF file to impute (can be used multiple times, required)
    -o, --output DIR    Output directory (required)
    -t, --token TOKEN   ImputationBot API token (required)
    -c, --cores NUM     Number of cores to use (default: 1)
    -n, --dry-run       Show what would be executed without running
    -h, --help          Show this help message

EXAMPLES:
    # Basic run with multiple files
    impute -i file1.vcf.gz -i file2.vcf.gz -o results -t your_api_token

    # Dry run to see what would be executed
    impute -i file1.vcf.gz -i file2.vcf.gz -o results -t your_api_token --dry-run

	# use 4 cores
	impute -i file1.vcf.gz -i file2.vcf.gz -o results -t your_api_token -c 4
```


## TODO

- [X] convert idat to gtc
- [X] convert gtc to vcf
- [X] liftover vcf to hg38 if necessary
- [X] submit to imputation server (see =https://github.com/marcoralab/imputePipeline)
- [ ] set BCFTOOLS tempdir
- [ ] add TopMed functionality
- [ ] add QC checks
   - [ ] check for sample duplicates
   - [ ] report callrates in idat2vcf
   - [ ] report imputation rates in impute

