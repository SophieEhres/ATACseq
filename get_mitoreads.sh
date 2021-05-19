#!/bin/bash

align_dir="${SCRATCH}/ATAC/align"
sorted_dir="${SCRATCH}/ATAC/sorted"
jobdir="${SCRATCH}/scripts/ATACseq/job_files/mitoreads"
outdir="${SCRATCH}/scripts/ATACseq/out_files/mitoreads"

mkdir -p ${sorted_dir}
mkdir -p ${jobdir}
mkdir -p ${outdir}

files=$(ls ${align_dir}/*.sam)

for file in ${files}; do

    name=$(echo "${file}"| xargs basename | cut -d "_" -f1) 

    if [ -f ${sorted_dir}/${name}_sorted-rmChrM.bam ]; then
        echo "${name} already mito read"
    else

echo "#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --job-name=${name}_getmitoreads
#SBATCH --output=${outdir}/${name}_getmitoreads.log
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=10000
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sophie.ehresmann@gmail.com
#SBATCH --account=def-sauvagm
    
echo "$(date +"%Y-%m-%d-%H-%M")" >> ${outdir}/${name}_getmitoreads.log
samtools --version >> ${outdir}/${name}_getmitoreads.log

samtools view -m 8G -h ${file} | grep -e "chrM" | wc -l > ${name}_mitoreads.txt
samtools view -m 8G -h ${file} | wc -l >> ${name}_mitoreads.txt

samtools view -m 8G -@ 20 -h ${file} | grep -v "chrM" | samtools sort -@ 20 -o ${sorted_dir}/${name}_sorted-rmChrM.bam -O bam 2>> ${outdir}/${name}_getmitoreads.log
samtools index -m 8G -@ 20 ${sorted_dir}/${name}_sorted-rmChrM.bam ${sorted_dir}/${name}_sorted-rmChrM.bai
" > ${jobdir}/${name}_getmitoreads.sh

    sbatch ${jobdir}/${name}_getmitoreads.sh

    fi

done