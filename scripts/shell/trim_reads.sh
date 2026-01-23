#!/usr/bin/env bash

fastqs="/Users/abrahamquaye/myocd_rnaseq/raw_files/merged_fastqs"
trimdir="/Users/abrahamquaye/myocd_rnaseq/results/trimmedReads"

forReads=( $fastqs/siO_*_R1_merged.fastq.gz )
revReads=( $fastqs/siO_*_R2_merged.fastq.gz )

if [ ${#forReads[@]} -ne ${#revReads[@]} ]; then
    echo "Number of forward reads do not match number of reverse reads"
    exit 1
fi

mkdir -p $trimdir

for n in ${!forReads[@]}; do
    fname=$(echo ${forReads[$n]} | cut -d "/" -f 7 | cut -d "." -f 1)
    fread=${forReads[$n]}
    rname=$(echo ${revReads[$n]} | cut -d "/" -f 7 | cut -d "." -f 1)
    rread=${revReads[$n]}

    echo "Trimming $fname and $rname ..."

    trim_galore --phred33 -q 20 --cores 8 --gzip \
    --path_to_cutadapt cutadapt --no_report_file --paired $fread $rread -o $trimdir
done

secs=$SECONDS
hrs=$(( secs/3600 )); mins=$(( (secs-hrs*3600)/60 )); secs=$(( secs-hrs*3600-mins*60 ))
printf 'Script completed in: %02d:%02d:%02d\n' $hrs $mins $secs
