 #!/bin/bash

fastq_dir="${SCRATCH}/ATAC/fastq"
trim_dir="${SCRATCH}/ATAC/trim"
tool_dir="${SCRATCH}/tools"
trimmomatic="${tool_dir}/Trimmomatic-0.39/trimmomatic-0.39.jar"
adapters="${tool_dir}/Trimmomatic-0.39/adapters"

mkdir -p ${trim_dir}

samples=$(ls ${fastq_dir} | rev | cut -d "_" -f2 | rev | uniq )


for sample in ${samples}; do 

    files="$(ls ${fastq_dir}/*.gz | grep -e ${sample} | tr '\n' ' ')"
    echo "files are ${files}" 
    
    echo "#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --job-name=Trimmomatic_${sample}
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-user=sophie.ehresmann@gmail.com
#SBATCH --mail-type=END

module load java

java -jar ${trimmomatic} PE -threads 20 -phred33 ${files} \
${trim_dir}/${sample}_forward_paired.fq.gz ${trim_dir}/${sample}_forward_unpaired.fq.gz \
${trim_dir}/${sample}_reverse_paired.fq.gz ${trim_dir}/${sample}_reverse_unpaired.fq.gz \
ILLUMINACLIP:${adapters}/IlluminaNextera_PE.fa:2:30:10
"> ${sample}_trimjob.sh 
sbatch --account=def-sauvagm ${sample}_trimjob.sh 

done