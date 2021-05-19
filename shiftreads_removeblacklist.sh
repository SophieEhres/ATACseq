#!/bin/bash

sortdir="${SCRATCH}/ATAC/sorted"
shiftdir="${SCRATCH}/ATAC/shift"
blacklist="${SCRATCH}/genomes/human/hg38/blacklist/lists/hg38-blacklist.v2.bed"
logdir="${SCRATCH}/scripts/ATACseq/out_files/shift"
jobdir="${SCRATCH}/scripts/ATACseq/job_files/shift"
metricdir="${SCRATCH}/ATAC/shiftmetrics"


mkdir -p ${logdir}
mkdir -p ${shiftdir}
mkdir -p ${jobdir}
mkdir -p ${metricdir}


for name in $(ls ${sortdir}/*.bam | xargs -n 1 basename | rev | cut -d "_" -f2- | rev); do

    file=${sortdir}/${name}_sorted-rmChrM.bam

    if [ -f ${shiftdir}/${name}_shift.bam ]; then
        echo "${name} already shifted"
    elif [ -f ${shiftdir}/${name}_shift_sort.bam ]; then
        echo "${name} already shifted and sorted"
    else

echo "#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --job-name=deeptools_shift_${name}
#SBATCH --output=${logdir}/${name}_shift.log
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=6000
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sophie.ehresmann@gmail.com
#SBATCH --account=def-sauvagm

echo "$(date +"%Y-%m-%d-%H-%M")" >> ${logdir}/${name}_shift.log
source ~/.bash_profile

module load python

deeptools --version >> ${logdir}/${name}_shift.log

alignmentSieve -b ${file} \
-p 20 --filterMetrics ${metricdir}/${name}.dTmetrics.txt \
--ATACshift -bl ${blacklist} -o ${shiftdir}/${name}_shift.bam


" > ${jobdir}/${name}_shift.sh
    echo "submitting shift job for ${name}"
sbatch --account=def-sauvagm ${jobdir}/${name}_shift.sh

    fi

done