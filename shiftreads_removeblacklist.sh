#!/bin/bash

sortdir="${SCRATCH}/ATAC/sorted"
shiftdir="${SCRATCH}/ATAC/shift"
blacklist="${SCRATCH}/genomes/human/hg38/blacklist/lists/hg38-blacklist.v2.bed"
metricdir="${SCRATCH}/ATAC/shiftmetrics"
logdir="${SCRATCH}/scripts/ATACseq/out_files/shift"
jobdir="${SCRATCH}/scripts/ATACseq/job_files/shift"

mkdir -p ${logdir}
mkdir -p ${shiftdir}
mkdir -p ${metricdir}
mkdir -p ${jobdir}


for name in $(ls ${sortdir}| rev | cut -d "_" -f2- | rev ); do

    file=${sortdir}/${name}_sorted-rmChrM.bam

echo "#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --job-name=deeptools_shift_${name}
#SBATCH --output=${logdir}/${name}_shift.log
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sophie.ehresmann@gmail.com
#SBATCH --account=def-sauvagm

echo "$(date +"%Y-%m-%d-%H-%M")" >> ${logdir}/${name}_shift.log

module load python

deeptools --version >> ${logdir}/${name}_shift.log

alignmentSieve -b ${file} \
-p 20 --filterMetrics ${metricdir}/${name}.dTmetrics.txt \
--ATACshift -bl ${blacklist} -o ${shiftdir}/${name}_shift.bam 2 >>  ${logdir}/${name}_shift.log

" > ${jobdir}/${name}_shift.sh

sbatch --account=def-sauvagm ${jobdir}/${name}_shift.sh

done