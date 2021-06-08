#!/bin/bash

matrixdir="${SCRATCH}/ATAC/matrix"
wigdir="${SCRATCH}/ATAC/bigwig"
jobdir="${SCRATCH}/scripts/ATACseq/job_files/matrix"
outdir="${SCRATCH}/scripts/ATACseq/out_files/matrix"
consensus="${SCRATCH}/ATAC/diffbind/consensus_peakset.bed"

mkdir -p ${matrixdir}
mkdir -p ${jobdir}
mkdir -p ${outdir}

days=$(ls ${wigdir} | cut -d "-" -f2 | uniq)

for day in ${days}; do
    files=$(ls ${wigdir}/*.bw | grep -e "${day}" | tr '\n' ' ')


    if [ -f ${matrixdir}}/${day}_matrix.gz ] ; then
   
        echo "${day} already computed matrix"
    
    else
    echo "creating job to compute matrix of ${day}"


echo "#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --job-name=matrix_${day}
#SBATCH --output=${outdir}/matrix_${day}.log
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=6G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sophie.ehresmann@gmail.com
#SBATCH --account=def-sauvagm

computeMatrix --version >> ${outdir}/matrix_${day}.log
computeMatrix scale-regions -S ${files} -R ${consensus} \
-o ${matrixdir}/${day}_scaleRegions.gz -p 20 \
-b 1500 -a 1500

" > ${jobdir}/${day}_matrix.sh

    sbatch --account=def-sauvagm ${jobdir}/${day}_matrix.sh 

    fi

done
    

