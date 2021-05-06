#!/bin/bash
now=$(date +"%Y-%m-%d-%H-%M")
echo "${now}" > bowtiebuild_${now}.log
bowtie2 --version >> bowtiebuild_${now}.log

genome="${SCRATCH}/genomes/human/hg38/fasta"
outdir="${SCRATCH}/genomes/human/hg38/bowtie_index"

mkdir -p ${outdir}

fastafiles=$(ls ${genome}/*.fa | tr '\n' ',')
echo "fasta files are $fastafiles"

bowtie2-build ${fastafiles} ${outdir}/hg38_ --threads 10 2 >> bowtiebuild_${now}.log
