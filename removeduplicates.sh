#!/bin/bash

markdupdir="${SCRATCH}/ATAC/markdup"
shiftdir="${SCRATCH}/ATAC/shift"
logdir="${SCRATCH}/scripts/ATACseq/out_files/dup"
jobdir="${SCRATCH}/scripts/ATACseq/job_files/dup"
metricdir="${SCRATCH}/ATAC/dupmetrics"
cleandir="${SCRATCH}/ATAC/cleanbam"

mkdir -p ${logdir}
mkdir -p ${markdupdir}
mkdir -p ${jobdir}
mkdir -p ${metricdir}
mkdir -p ${cleandir}


for name in $(ls ${shiftdir}| grep -v "sort" | rev | cut -d "_" -f2- | rev ); do

    file=${shiftdir}/${name}_shift.bam

    if [ -f ${cleandir}/${name}_clean.bam ]; then
        echo "${name} already removed duplicates"
    else

echo "#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --job-name=removeduplicates_${name}
#SBATCH --output=${logdir}/${name}_dup.log
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=12G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sophie.ehresmann@gmail.com
#SBATCH --account=def-sauvagm

echo "$(date +"%Y-%m-%d-%H-%M")" >> ${logdir}/${name}_dup.log

module load picard java/11.0.2

samtools sort -m 8G -@ 20 -o ${shiftdir}/${name}_shift_sort.bam -O bam ${file}

if [ -f ${shiftdir}/${name}_shift_sort.bam ]; then

echo "file sorted" >>  ${logdir}/${name}_dup.log
rm ${file}


java -jar /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/picard/2.23.3/picard.jar MarkDuplicates --version >> ${logdir}/${name}_dup.log

java -Xmx20G -jar /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/picard/2.23.3/picard.jar MarkDuplicates -I ${shiftdir}/${name}_shift_sort.bam -M ${metricdir}/${name}_dupmetrics.txt \
-O ${markdupdir}/${name}_markdup.bam \
--REMOVE_DUPLICATES false

samtools view -m 8G  -F 1804 ${markdupdir}/${name}_markdup.bam > ${cleandir}/${name}_clean.bam

else

echo "file not correctly sorted, exiting" >>  ${logdir}/${name}_dup.log

fi

" > ${jobdir}/${name}_dup.sh
    
echo "submitting dup job for ${name}"

sbatch --account=def-sauvagm ${jobdir}/${name}_dup.sh

    fi

done