#!/bin/bash
date>log_align_${date}.txt


computational="/Users/ehresms/computational"
genome="${computational}/genomes/human/bowtie2/GRCh38_noalt_as/GRCh38_noalt_as"
trim_dir="${computational}/ATAC/trim"
align_dir="${computational}/ATAC/align"

mkdir -p ${align_dir}

samples=$(ls ${trim_dir} | cut -d "_" -f1 | uniq )
echo "samples are ${samples}"

for sample in ${samples}; do

if ls ${align_dir} | grep -e ${sample}; then
    echo "${sample} is already aligned"
else 
echo "aligning ${sample}"

file_1=$(ls ${trim_dir}/*.gz | grep -e ${sample} | grep -e "forward_paired")
file_2=$(ls ${trim_dir}/*.gz | grep -e ${sample} | grep -e "reverse_paired")

file_f=$(echo ${file_1} | cut -d "." -f1-2)
file_r=$(echo ${file_2} | cut -d "." -f1-2)

gunzip -k ${file_1}
gunzip -k ${file_2}

bowtie2 -x ${genome} -1 ${file_f} -2 ${file_r} \
--phred33 --very-sensitive -X 1000 --no-discordant \
-t --threads 30 -S ${align_dir}/${sample}_aligned.sam 2>>log_align.txt

rm ${file_f}
rm ${file_r}
echo "done aligning ${sample}"

fi
done

