# Snakemake pipeline to convert to IDAT files to vcf for imputation

TODO: 
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

