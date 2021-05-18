#!/bin/bash

cleandir="${SCRATCH}/ATAC/clean"
shiftdir="${SCRATCH}/ATAC/shift"
logdir="${SCRATCH}/scripts/ATACseq/out_files/dup"
jobdir="${SCRATCH}/scripts/ATACseq/job_files/dup"
metricdir="${SCRATCH}/ATAC/dupmetrics"


mkdir -p ${logdir}
mkdir -p ${cleandir}
mkdir -p ${jobdir}
mkdir -p ${metricdir}


for name in $(ls ${shiftdir}| rev | cut -d "_" -f2- | rev ); do

    file=${shiftdir}/${name}_shift.bam

    if [ -f ${cleandir}/${name}_dup.bam ]; then
        echo "${name} already removed duplicates"
    else

echo "#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --job-name=removeduplicates_${name}
#SBATCH --output=${logdir}/${name}_dup.log
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=6000
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sophie.ehresmann@gmail.com
#SBATCH --account=def-sauvagm

echo "$(date +"%Y-%m-%d-%H-%M")" >> ${logdir}/${name}_dup.log

module load picard java/11.0.2

samtools sort -@ 20 -o ${shiftdir}/${name}_shift_sort.bam -O bam ${file}

if [ -z $(samtools quickcheck ${shiftdir}/${name}_shift_sort.bam) ]; then

echo "file correctly sorted" >>  ${logdir}/${name}_dup.log
rm ${file}


java -jar /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/picard/2.23.3/picard.jar MarkDuplicates --version >> ${logdir}/${name}_dup.log

java -Xmx20G -jar /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/picard/2.23.3/picard.jar MarkDuplicates -I ${shiftdir}/${name}_shift_sort.bam -M ${metricdir}/${name}_dupmetrics.txt \
-O ${cleandir}/${name}_clean.bam \
--REMOVE_DUPLICATES true

else

echo "file not correctly sorted, exiting" >>  ${logdir}/${name}_dup.log

fi

" > ${jobdir}/${name}_dup.sh

sbatch --account=def-sauvagm ${jobdir}/${name}_dup.sh

    fi

done