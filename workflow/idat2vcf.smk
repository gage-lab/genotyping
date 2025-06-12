from pathlib import Path
from snakemake.utils import min_version
from array_config import read_samplesheet

##### set minimum snakemake version #####
min_version("9.5.0")

# read config from cli
outdir, igm_dir = config["outdir"], config["igm_dir"]
assert Path(igm_dir).exists(), f"igm_dir {igm_dir} does not exist"
Path(outdir).mkdir(parents=True, exist_ok=True)

# read samplesheet
samples, array_config = read_samplesheet(igm_dir)


rule idat2gtc:
    input:
        green_idats=array_config["Grn.idat"],
        red_idats=array_config["Red.idat"],
    output:
        green_idats=temp(outdir + "/green_idats.txt"),
        red_idats=temp(outdir + "/red_idats.txt"),
        outdir=temp(directory(outdir + "/gtc")),
        outfiles=temp(expand(outdir + "/gtc/{sample}.gtc", sample=samples["sentrix_id"])),
    params:
        **array_config,
    log:
        outdir + "/idat2gtc.log",
    shell:
        """
        export BCFTOOLS_PLUGINS="$CONDA_PREFIX/libexec/bcftools"
        mkdir -p {output.outdir}

        for i in {input.green_idats}; do echo $i >> {output.green_idats}; done
        for i in {input.red_idats}; do echo $i >> {output.red_idats}; done

        bcftools plugin $BCFTOOLS_PLUGINS/idat2gtc.so \
            --bpm {params.bpm_manifest_file} \
            --egt {params.egt_cluster_file} \
            --grn-idats {output.green_idats} \
            --red-idats {output.red_idats} \
            --output {output.outdir} 2> {log}
        """


rule gtc2vcf:
    input:
        gtc=rules.idat2gtc.output.outfiles,
    output:
        names=temp(outdir + "/reheader_names.txt"),
        gtc=temp(outdir + "/gtc_files.txt"),
        sentrix_vcf=outdir + "/genotypes_sentrix.vcf.gz",
        sentrix_vcf_idx=outdir + "/genotypes_sentrix.vcf.gz.csi",
        sample_vcf=outdir + "/genotypes_sample.vcf.gz",
        sample_vcf_idx=outdir + "/genotypes_sample.vcf.gz.csi",
    params:
        **array_config,
        reheader="\n".join(
            [
                f"{sent}\t{sample}"
                for sent, sample in zip(samples["sentrix_id"], samples["Sample_ID"])
            ]
        ),
    log:
        outdir + "/gtc2vcf.log",
    shell:
        """
        # write all stderr to log
        exec 2>> {log}

        # set bcftools plugins
        export BCFTOOLS_PLUGINS="$CONDA_PREFIX/libexec/bcftools"

        # write gtc files to file
        for i in {input.gtc}; do echo $i >> {output.gtc}; done

        bcftools plugin $BCFTOOLS_PLUGINS/gtc2vcf.so \
            --no-version -Ou \
            --bpm {params.bpm_manifest_file} \
            --csv {params.csv_manifest_file} \
            --egt {params.egt_cluster_file} \
            --fasta-ref {params.ref} \
            --gtcs {output.gtc} | \
        bcftools sort -Ou | \
        bcftools filter -i 'FORMAT/GT != "."' | \
        bcftools norm --no-version -o {output.sentrix_vcf} -Oz -c x -f {params.ref} --write-index

        # write reheader names
        echo -e "{params.reheader}" > {output.names}
        bcftools reheader -s {output.names} {output.sentrix_vcf} > {output.sample_vcf}
        bcftools index {output.sample_vcf}
        """


rule liftover:
    input:
        vcf=rules.gtc2vcf.output.sample_vcf,
        chain="resources/hg19ToHg38.over.chain",
        hg38="resources/Homo_sapiens_assembly38.fasta",
        chr_mapping="resources/chr_mapping.txt",
    output:
        chr_vcf=temp(outdir + "/genotypes_sample.hg19.chr.vcf.gz"),
        vcf=outdir + "/genotypes_sample.hg19_to_hg38.vcf.gz",
        rejected=outdir + "/liftover_rejected.vcf",
    params:
        tmpdir=outdir + "/liftover_tmp",
        max_records_in_ram=100000,
    log:
        outdir + "/liftover.log",
    shell:
        """
        exec 2>> {log}
        bcftools annotate --rename-chrs {input.chr_mapping} -o {output.chr_vcf} -Oz {input.vcf}

        gatk LiftoverVcf \
            -C {input.chain} \
            -I {output.chr_vcf} \
            -R {input.hg38} \
            --CREATE_INDEX \
            --TMP_DIR {params.tmpdir} \
            --REJECT {output.rejected} \
            --MAX_RECORDS_IN_RAM {params.max_records_in_ram} \
            -O {output.vcf} 
        """


if array_config["ref"] == "resources/Homo_sapiens_assembly38.fasta":
    final = rules.gtc2vcf.output.sample_vcf
elif array_config["ref"] == "resources/Homo_sapiens_assembly19.fasta":
    final = rules.liftover.output.vcf
else:
    raise ValueError(f"Reference fasta {array_config['ref']} not supported")


onsuccess:
    samples[["sentrix_id", "Sample_ID"]].to_csv(
        outdir + "/samplesheet.tsv", sep="\t", index=False
    )


rule all:
    input:
        final,
