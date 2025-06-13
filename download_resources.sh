#!/usr/bin/env bash
# Author: Mike Cuoco
# Created on: 6/12/25, 9:03 PM
#
# Download resources for the idat2vcf and impute workflows

# exit if any non-zero, exit if undefined var
set -euo pipefail

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source common functions
mkdir -p "${SCRIPT_DIR}"/resources && cd "${SCRIPT_DIR}"/resources

PREFIX="https://storage.googleapis.com/gcp-public-data--broad-references/"
files=(
	"${PREFIX}/hg38/v0/Homo_sapiens_assembly38.fasta"
	"${PREFIX}/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
	"${PREFIX}/hg38/v0/Homo_sapiens_assembly38.dict"
	"${PREFIX}/hg19/v0/Homo_sapiens_assembly19.fasta"
	"${PREFIX}/hg19/v0/Homo_sapiens_assembly19.fasta.fai"
	"${PREFIX}/hg19/v0/Homo_sapiens_assembly19.dict"
)

for file in "${files[@]}"; do
	echo "Downloading $file"
	curl -L -s -S -o "${file##*/}" "${file}" >/dev/null 2>&1
done

mkdir -p illumina_array_support && cd illumina_array_support

files=(
	"https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/global-screening-array-24/v3-0/GSA-24v3-0-A2-manifest-file-csv.zip"
	"https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/global-screening-array-24/v3-0/GSA-24v3-0-A2-manifest-file-bpm.zip"
	"https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/global-screening-array-24/v3-0/GSA-24v3-0-A1-cluster-file.zip"
	"https://webdata.illumina.com/downloads/productfiles/humancoreexome/24-v1-4/InfiniumCoreExome-24v1-4_A1-csv.zip"
	"https://webdata.illumina.com/downloads/productfiles/humancoreexome/24-v1-4/InfiniumCoreExome-24v1-4_A1-bpm.zip"
	"https://webdata.illumina.com/downloads/productfiles/humancoreexome/24-v1-4/InfiniumCoreExome-24v1-4_A1_ClusterFile.egt"
)

for file in "${files[@]}"; do
	echo "Downloading $file"
	curl -L -s -S -o "${file##*/}" "${file}"
	if [[ "${file##*/}" == *.zip ]]; then
	    echo "Unzipping ${file##*/}"
		unzip -o -q "${file##*/}"
		rm -f "${file##*/}"
	fi
done

