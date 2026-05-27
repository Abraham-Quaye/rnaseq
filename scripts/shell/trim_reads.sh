#!/opt/homebrew/bin/bash

fastqdir="/Users/abrahamquaye/bulk_slc38a9_rnaseq_aq/raw_files/raw_fastqs"
trimdir="/Users/abrahamquaye/bulk_slc38a9_rnaseq_aq/results/trimmedReads"

forReads=( $fastqdir/*_R1_001.fastq.gz )
revReads=( $fastqdir/*_R2_001.fastq.gz )

if [ ${#forReads[@]} -ne ${#revReads[@]} ]; then
    echo "Number of forward reads do not match number of reverse reads"
    exit 1
fi

mkdir -p $trimdir

for n in ${!forReads[@]}; do
    fname=$(echo ${forReads[$n]##*/} | cut -d "_" -f1,2)
    fread=${forReads[$n]}
    rname=$(echo ${revReads[$n]##*/} | cut -d "_" -f1,2)
    rread=${revReads[$n]}

    echo "Trimming $fname and $rname ..."

    trim_galore --phred33 -q 20 --cores 8 --gzip \
    --path_to_cutadapt cutadapt --no_report_file --paired $fread $rread -o $trimdir
done

secs=$SECONDS
hrs=$(( secs/3600 )); mins=$(( (secs-hrs*3600)/60 )); secs=$(( secs-hrs*3600-mins*60 ))
printf 'Script completed in: %02d:%02d:%02d\n' $hrs $mins $secs
