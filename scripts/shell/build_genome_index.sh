#!/usr/bin/env zsh


anno_dir=results/annotations
genomedir=results/genome_files
idxdir=results/hisat_genome_index

################# COMMAND TO BUILD hg38 GENOMIC INDEX WITH HISAT2 #############
echo "hg38 genomic index ..."
hisat2-build -p 12 --ss $anno_dir/hg38.ss --exon $anno_dir/hg38.exons\
$genomedir/hg38_genome.fa.gz $idxdir/hisat_hg38_tran &&

echo "Index built successfully" || echo "Error in building index"
