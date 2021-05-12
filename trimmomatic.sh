#!/bin/bash
2>&1
date>log_trim.txt

fastq_dir="/Users/ehresms/computational/ATAC/FASTQ"
samples=$(ls ${fastq_dir} grep -e "fastq." | rev | cut -d "_" -f2 | rev | uniq)
trimmomatic="/Users/ehresms/computational/tools/Trimmomatic-0.39/trimmomatic-0.39.jar"
adapters="/Users/ehresms/computational/tools/Trimmomatic-0.39/adapters"
trim_dir="/Users/ehresms/computational/ATAC/trim"

mkdir -p ${trim_dir}


for sample in ${samples}; do
    
    files=$(ls ${fastq_dir}/*.gz | grep -e ${sample} | tr '\n' ' ')
    java -jar ${trimmomatic} PE -threads 50 -phred33 ${files} \
        ${trim_dir}/${sample}_forward_paired.fq.gz ${trim_dir}/${sample}_forward_unpaired.fq.gz \
        ${trim_dir}/${sample}_reverse_paired.fq.gz ${trim_dir}/${sample}_reverse_unpaired.fq.gz \
        ILLUMINACLIP:${adapters}/IlluminaNextera_PE.fa:2:30:10
done >>log_trim.txt
        

