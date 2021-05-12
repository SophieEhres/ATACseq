#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --job-name=getmitoreads
#SBATCH --output=getmitoreads.log
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sophie.ehresmann@gmail.com
#SBATCH --account=def-sauvagm

now=$(date +"%Y-%m-%d-%H-%M")
echo "${now}" >getmitoreads.log

align_dir="${SCRATCH}/ATAC/align"
sorted_dir="${SCRATCH}/ATAC/sorted"

mkdir -p ${sorted_dir}

files=$(ls ${align_dir}/*.sam)

for file in ${files}; do

    name=$(echo "${file}"| xargs basename | cut -d "_" -f1)
    mitoreads=$(samtools view -h ${file} | grep -e "chrM" | wc -l)
    totalreads=$(samtools view -h ${file} | wc -l)

    sum=$(expr ${mitoreads} * 100)
    percentage=$(expr ${sum} / ${totalreads})

    echo "For sample ${name}, there are ${mitoreads} mitochondrial reads for a total of ${totalreads}, accounting for ${percentage} % of reads" >> getmitoreads.log

    samtools view -@ 30 -h ${file} | grep -v "chrM" | samtools sort -@ 30 -o ${sorted_dir}/${name}_sorted-rmChrM.sam 2>> getmitoreads.log

done