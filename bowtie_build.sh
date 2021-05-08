#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --job-name=bowtie_build
#SBATCH --output=bowtie_build.out
#SBATCH --error bowtie_build.err
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sophie.ehresmann@gmail.com
#SBATCH --account=def-sauvagm

now=$(date +"%Y-%m-%d-%H-%M")
echo "${now}" > bowtiebuild_${now}.log
bowtie2 --version >> bowtiebuild_${now}.log

genome="${SCRATCH}/genomes/human/hg38/fasta"
outdir="${SCRATCH}/genomes/human/hg38/bowtie_index"

mkdir -p ${outdir}

fastafiles=$(ls ${genome}/*.fa | tr '\n' ',')
echo "fasta files are $fastafiles"

bowtie2-build ${fastafiles} ${outdir}/hg38_ --threads 20 2 >> bowtiebuild_${now}.log
