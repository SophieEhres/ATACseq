#!/bin/bash

genome="${SCRATCH}/genomes/human/hg38/bowtie_index/hg38_"
trimdir="${SCRATCH}/ATAC/trim"
aligndir="${SCRATCH}/ATAC/align"
jobdir="${SCRATCH}/scripts/ATACseq/job_files/align_jobs"
outdir="${SCRATCH}/scripts/ATACseq/out_files"

mkdir -p ${aligndir}
mkdir -p ${jobdir}

samples=$(ls ${trimdir} | cut -d "_" -f1 | uniq)
echo "samples are ${samples}"

for sample in ${samples}; do

file_1=$(ls ${trimdir}/*.gz | grep -e "${sample}" | grep -e "forward_paired")
file_2=$(ls ${trimdir}/*.gz | grep -e "${sample}" | grep -e "reverse_paired")


file_f=$(echo ${file_1} | cut -d "." -f1-2)
file_r=$(echo ${file_2} | cut -d "." -f1-2)

echo "#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --job-name=bowtie_align_${sample}
#SBATCH --output=${outdir}/bowtie_align_${sample}.log
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sophie.ehresmann@gmail.com
#SBATCH --account=def-sauvagm

now=$(date +"%Y-%m-%d-%H-%M")
echo "${now}" > ${outdir}/bowtie_align_${sample}.log
bowtie2 --version >> ${outdir}/bowtie_align_${sample}.log

gunzip -k ${file_1}
gunzip -k ${file_2}

bowtie2 -x ${genome} -1 ${file_f} -2 ${file_r} \
--phred33 --very-sensitive -X 1000 --no-discordant \
-t --threads 30 -S ${aligndir}/${sample}_aligned.sam 2>> ${outdir}/bowtie_align_${sample}.log

rm ${file_f}
rm ${file_r}
echo "done aligning ${sample}"

" > ${jobdir}/${sample}_align.sh

done

for file in $(ls ${jobdir}); do
sbatch ${jobdir}/${file}
done