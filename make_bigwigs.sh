#!/bin/bash

cleandir="${SCRATCH}/ATAC/cleanbam"
wigdir="${SCRATCH}/ATAC/bigwig"
jobdir="${SCRATCH}/scripts/ATACseq/job_files/bigwig"
outdir="${SCRATCH}/scripts/ATACseq/out_files/bigwig"

mkdir -p ${wigdir}
mkdir -p ${jobdir}
mkdir -p ${outdir}

for file in $(ls ${cleandir}/*.bam); do
    name=$(echo ${file} | xargs -n 1 basename | cut -d "." -f1)
    
    if [ -f ${wigdir}}/${name}_clean.bw ] ; then
   
        echo "${name} already made bigwig file"
    
    else
    echo "creating job to make bigwig of ${name}"


echo "#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --job-name=bigwig_${name}
#SBATCH --output=${outdir}/bigwig_${name}.log
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=6G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sophie.ehresmann@gmail.com
#SBATCH --account=def-sauvagm

bamCoverage --version >> ${outdir}/bigwig_${name}.log
bamCoverage -b ${file} \
-o ${wigdir}/${name}.bw -p 20 \
--normalizeUsing RPKM

" > ${jobdir}/${name}_bigwig.sh

    sbatch --account=def-sauvagm ${jobdir}/${name}_bigwig.sh 

    fi

done
    