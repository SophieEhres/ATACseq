#!/bin/bash

cleandir="${SCRATCH}/ATAC/cleanbam"
peakdir="${SCRATCH}/ATAC/peak"
logdir="${SCRATCH}/scripts/ATACseq/out_files/peak"
jobdir="${SCRATCH}/scripts/ATACseq/job_files/peak"


mkdir -p ${logdir}
mkdir -p ${jobdir}
mkdir -p ${peakdir}


for name in $(ls ${cleandir} | rev | cut -d "_" -f2- | rev); do

 file=${cleandir}/${name}_clean.bam

    if [ -f ${peakdir}/${name}/${name}_peaks.xls ]; then
        echo "${name} already peak called"
    else

echo "#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --job-name=peakcall_${name}
#SBATCH --output=${logdir}/${name}_peak.log
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=20G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sophie.ehresmann@gmail.com
#SBATCH --account=def-sauvagm

module load python scipy-stack

echo "$(date +"%Y-%m-%d-%H-%M")" >> ${logdir}/${name}_dup.log

macs2 --version >> ${logdir}/${name}_dup.log

macs2 callpeak -t ${file} -n ${name} --outdir ${peakdir}/${name} \
-f BAMPE --nomodel --keep-dup all --nolambda 2>>${logdir}/${name}_dup.log

" > ${jobdir}/${name}_peak.sh
    
echo "submitting peak job for ${name}"

sbatch --account=def-sauvagm ${jobdir}/${name}_peak.sh

    fi

done