#!/bin/bash
2>&1
1>log_checksums.txt
date

dir="/Users/ehresms/computational/ATAC/FASTQ/readSetLinks"

for file in $(ls ${dir} | grep -v ".md5" | grep -e ".gz"); do
	file_md5=$(md5 ${dir}/${file} | awk '{print $4}')
	md5=$(cat ${dir}/${file}.md5 | awk '{print $1}')
	
	if [[ "${file_md5}" == "${md5}" ]] ; then
		echo "${file} checked and ok"
	else 
		echo "${file} checked and not ok \\n file_md5 is ${file_md5} and md5 is ${md5}"
	fi
done

