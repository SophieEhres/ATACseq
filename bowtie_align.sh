#!/bin/bash
2>&1
date>log_trim.txt


computational="/Users/ehresms/computational"
genome="${computational}/genomes/human/bowtie2"
trim_dir="${computational}/ATAC/trim"
align_dir="${computational}/ATAC/align"

mkdir -p ${align_dir}

samples=$(ls ${trim_dir} | cut -d "_" -f1 | uniq | head -n 1)
echo "samples are ${samples}"

for sample in ${samples}; do
echo "aligning ${sample}"

file_1=$(ls ${trim_dir}/*.gz | grep -e ${sample} | grep -e "forward_paired")
file_2=$(ls ${trim_dir}/*.gz | grep -e ${sample} | grep -e "reverse_paired")

file_f=$(echo ${file_1} | cut -d "." -f1-2)
file_r=$(echo ${file_2} | cut -d "." -f1-2)

gunzip ${file_1}
gunzip ${file_2}

bowtie2 -x ${genome} -1 ${file_f} -2 ${file_r} \
--phred33 --very-sensitive -X 1000 --no-discordant \
-t --threads 50 -S ${align_dir}/${sample}_aligned.sam

rm ${file_f}
rm ${file_r}
echo "done aligning ${sample}"

done >>log_trim.txt

