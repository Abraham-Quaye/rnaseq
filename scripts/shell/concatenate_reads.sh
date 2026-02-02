#!/usr/bin/env bash

raw_dir="/Users/abrahamquaye/myocd_rnaseq/raw_files"

fastq_dirs=( $raw_dir/siO_* )

mkdir -p $raw_dir/merged_fastqs

 echo "Concatenating same-strand reads for all samples ..."
for d in ${fastq_dirs[@]}; do
    dir_name=$(basename $d)
    echo $dir_name
    cat $d/${dir_name}_L00*_R1_001.fastq.gz > $raw_dir/merged_fastqs/${dir_name}_R1_merged.fastq.gz
    cat $d/${dir_name}_L00*_R2_001.fastq.gz > $raw_dir/merged_fastqs/${dir_name}_R2_merged.fastq.gz
done