#!/bin/bash

align_dir="${SCRATCH}/ATAC/align"
sorted_dir="${SCRATCH}/ATAC/sorted"
jobdir="${SCRATCH}/scripts/ATACseq/job_files/mitoreads"
outdir="${SCRATCH}/scripts/ATACseq/out_files/mitoreads"

mkdir -p ${sorted_dir}
mkdir -p ${jobdir}

files=$(ls ${align_dir}/*.sam)

for file in ${files}; do

name=$(echo "${file}"| xargs basename | cut -d "_" -f1) 

echo "#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --job-name=${name}_getmitoreads
#SBATCH --output=${outdir}/${name}_getmitoreads.log
#SBATCH --ntasks=30
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sophie.ehresmann@gmail.com
#SBATCH --account=def-sauvagm
    
echo "$(date +"%Y-%m-%d-%H-%M")" >> ${outdir}/${name}_getmitoreads.log
samtools --version >> ${outdir}/${name}_getmitoreads.log
samtools view -h ${file} | grep -e "chrM" | wc -l > mitoreads
samtools view -h ${file} | wc -l > totalreads

echo "For name ${name}, there are ${mitoreads} mitochondrial reads for a total of ${totalreads}" >> ${outdir}/${name}_getmitoreads.log

samtools view -@ 30 -h ${file} | grep -v "chrM" | samtools sort -@ 30 -o ${sorted_dir}/${name}_sorted-rmChrM.bam -O bam 2>> ${outdir}/${name}_getmitoreads.log
" > ${jobdir}/${name}_getmitoreads.sh

sbatch ${jobdir}/${name}_getmitoreads.sh
done