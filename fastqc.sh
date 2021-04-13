#!/bin/bash
2>&1
date>log_fastqc.txt

fastq_dir="/Users/ehresms/computational/ATAC/FASTQ"
files=$(ls ${fastq_dir} | grep -e ".gz")
report_dir="/Users/ehresms/computational/ATAC/fastqc_reports"
echo "files are $files"
mkdir $report_dir


for file in ${files}; do
    name=$(echo ${file} | rev | cut -d "_" -f1-2 | rev)
    echo "name is $name and file is $dir/$file"
    gunzip -c ${fastq_dir}/${file} | fastqc -o ${report_dir} -t 12 stdin:${name}

done >> log_fastqc.txt