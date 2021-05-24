#!/bin/bash
statfile="$SCRATCH/ATAC/ATACstats.txt"

outdir="${SCRATCH}/scripts/ATACseq/out_files"
aligndir="${outdir}/align"
mitodir="${outdir}/mitoreads"
blacklistmetdir="${SCRATCH}/ATAC/shiftmetrics"
dupmetdir="${SCRATCH}/ATAC/dupmetrics"
peakdir="${SCRATCH}/ATAC/peak"

echo -e "Sample\tInitial_Fragment_Number\tAlignmentRate\tMitochondrial_Reads\tBlacklist_Reads\tDuplicated_Percentage\tRemaining_Fragment_Number\tPeak_Number" > ${statfile}

for name in $(ls ${peakdir} | head -n 1); do
    readnumber=$(cat ${aligndir}/bowtie_align_${name}.log | grep -e "reads" | awk '{print $1}')

    alignmentrate=$(cat ${aligndir}/bowtie_align_${name}.log | grep -e "rate" | awk '{print $1}')

    mitoreads=$(cat ${mitodir}/${name}_mitoreads.txt | head -n 1)

    blacklistreads=$(sed '3q;d' ${blacklistmetdir}/${name}.dTmetrics.txt | awk '{print $2}')
    totalreads=$(sed '3q;d' ${blackli stmetdir}/${name}.dTmetrics.txt | awk '{print $3}')
    removedblacklist=$(expr ${totalreads} - ${blacklistreads})
    removedfragments=$(expr ${removedblacklist} \\ 2)

    echo "${blacklistreads}, ${totalreads}, ${removedblacklist}, ${removedfragments}"
    dupreads=$(sed '8q;d' ${dupmetdir}/${name}_dupmetrics.txt | awk '{print $10}')

    peakfile=${peakdir}/${name}/*.xls
    fragments=$(sed '20q;d' ${peakfile} | awk '{print $6}')
    peaks=$(tail -n 1 ${peakfile} | awk '{print $10}' | rev | cut -d "_" -f1 | rev)


    echo -e "${name}\t${readnumber}\t${alignmentrate}\t${mitoreads}\t${removedblacklist}\t${dupreads}\t${fragments}\t${peaks}" >> ${statfile}
    echo "got stats for ${name}"
done

