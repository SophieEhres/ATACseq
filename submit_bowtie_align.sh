#!/bin/bash

genome="${SCRATCH}/genomes/bowtie"
trim_dir="${SCRATCH}/ATAC/trim"
job_dir="${SCRATCH}/scripts/ATACseq/job_files/align"
align_dir="${SCRATCH}/ATAC/align"

mkdir -p ${job_dir}
mkdir -p ${align_dir}

samples=$(ls ${trim_dir} | cut -d "_" -f1 | uniq | head -n 1)
echo "samples are ${samples}"

for sample in ${samples}; do
file_1=$(ls ${trim_dir}/*.gz | grep -e ${sample} | grep -e "forward_paired")
file_2=$(ls ${trim_dir}/*.gz | grep -e ${sample} | grep -e "reverse_paired")

file_f=$(echo ${file_1} | cut -d "." -f1-2)
file_r=$(echo ${file_2} | cut -d "." -f1-2)

echo "#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --job-name=bowtie_align_${sample}
#SBATCH --output=${sample}_align.out
#SBATCH --error ${sample}_align.err
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sophie.ehresmann@gmail.com

module load python

gunzip ${file_1}
gunzip ${file_2}

python bowtie2 -x ${genome} -1 ${file_f} -2 ${file_r} \
--phred33 --very-sensitive -X 1000 --no-discordant \
-t --threads 20 -S ${align_dir}/${sample}_aligned.sam

rm ${file_f}
rm ${file_r}

" > ${job_dir}/${sample}_submit_align.sh

done

for job_file in $(ls ${job_dir}); do
    sbatch --account=def-sauvagm ${job_dir}/${job_file}
    wait 1
done