#!/usr/bin/env zsh


filedir=results/annotations


################# EXTRACT EXONS FROM TURKEY GTF FILE 
echo "TURKEY extracting exons ..."
hisat2_extract_exons.py $filedir/turkey_genome.gtf > $filedir/turkey_genome.exons &&
echo "TURKEY exon extration complete"

cat $filedir/thev_predicted_genes.exons $filedir/turkey_genome.exons > $filedir/host_thev.exons

