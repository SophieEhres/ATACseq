#!/bin/bash
now=$(date +"%Y-%m-%d-%H-%M")
echo "${now}" > bowtiebuild_${now}.log

genome="${SCRATCH}/genomes/human/hg38/fasta"
outdir="${SCRATCH}/genomes/human/hg38/bowtie_index"

fastafiles=$(ls ${genome}/*.fa)

bowtie2-build ${fastafiles} ${outdir} 2 > bowtiebuild_${now}.log
