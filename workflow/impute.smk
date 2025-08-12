from pathlib import Path
from snakemake.utils import min_version
from array_config import read_samplesheet
import random, string

##### set minimum snakemake version #####
min_version("9.5.0")

outdir = config["outdir"]
vcfs = config["vcfs"]
if not isinstance(vcfs, list):
	vcfs = [vcfs]
print(f"vcfs: {vcfs}")
token = config["token"]
CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X"]
# CHROMOSOMES = ["22"] # for testing

rule combine_vcfs:
	input:
		vcfs = vcfs
	output:
		vcf = outdir + "/preimpute.vcf.gz"
	log:
		outdir + "/combine_vcfs.log"
	shell:
		"""
		bcftools merge -m all -O z -o {output.vcf} {input.vcfs} 2> {log}
		"""

rule chunk:
	input:
		vcf = rules.combine_vcfs.output.vcf if len(vcfs) > 1 else vcfs[0]
	output:
		vcf = outdir + "/chrs/preimpute.chr{chrom}.vcf.gz"
	log:
		outdir + "/chrs/preimpute.chunk.chr{chrom}.log"
	shell:
		"""
		bcftools view -r chr{wildcards.chrom} {input.vcf} -O z -o {output.vcf} 2> {log}
		"""


rule impute_on_server:
	input:
		vcf = expand(
			rules.chunk.output.vcf,
			chrom=CHROMOSOMES
		),
	output:
		vcfs = temp(expand(
			outdir + '/imputed/chr{chrom}.dose.vcf.gz',
			chrom=CHROMOSOMES
		)),
		infos = temp(expand(
			outdir + '/imputed/chr{chrom}.info.gz',
			chrom=CHROMOSOMES
		)),
		qc_report = outdir + "/imputed/qc_report.txt",
		quality_report = outdir + "/imputed/quality-control.html",
		job_json = outdir + "/imputed/job.json"
	params:
		refpanel = "1000g-phase3-deep",
		population = "all",
		build = "hg38",
		token = token,
		server = "michigan",
	log:
		outdir + "/impute.log"
	script:
		"impute_on_server.py"

rule concat_imputed:
	input:
		vcf = expand(
			outdir + '/imputed/chr{chrom}.dose.vcf.gz',
			chrom=CHROMOSOMES
		)
	output:
		unsorted_vcf = temp(outdir + "/imputed/unsorted.vcf.gz"),
		sorted_vcf = outdir + "/genotypes_imputed.vcf.gz"
	shell:
		"""
		bcftools concat -Oz -o {output.unsorted_vcf} {input.vcf}
		bcftools sort {output.unsorted_vcf} -Oz -o {output.sorted_vcf}
		tabix -p vcf {output.sorted_vcf}
		"""

rule all:
	input:
		rules.concat_imputed.output.sorted_vcf