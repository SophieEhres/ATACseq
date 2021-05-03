#!/bin/bash
now=$(date +"%Y-%m-%d-%H-%M")
echo "${now}" >log_mito_${now}

computational="/Users/ehresms/computational"
align_dir="${computational}/ATAC/align"
sorted_dir="${computational}/ATAC/sorted"

mkdir -p ${sorted_dir}

files=$(ls ${align_dir}/*.sam)

for file in ${files}; do

    name=$(echo "${file}"| xargs basename | cut -d "_" -f1)
    mitoreads=$(samtools view -h ${file} | grep -e "chrM" | wc -l)
    totalreads=$(samtools view -h ${file} | wc -l)

    sum=$(expr ${mitoreads} * 100)
    percentage=$(${sum} / ${totalreads})

    echo "there are ${mitoreads} mitochondrial reads for a total of ${totalreads}, accounting for ${percentage} % of reads" >> log_mito_${now}
    samtools view -@ 30 -h ${file} | grep -v "chrM" | samtools sort -@ 30 -o ${sorted_dir}/${name}_sorted-rmChrM.sam 2>>log_mito_${now}
done