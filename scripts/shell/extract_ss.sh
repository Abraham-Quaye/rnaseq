#!/usr/bin/env zsh


filedir=raw_files/annotations

################# EXTRACT SPLICESITES FROM THEV GTF FILE ################
echo "extracting THEV splice-sites ..."
hisat2_extract_splice_sites.py $filedir/thev_predicted_genes.gtf > $filedir/thev_predicted_genes.ss &&
echo "THEV splice-site extraction completed successfully"

################# EXTRACT SPLICESITES FROM THEV GTF FILE ################
echo "extracting TURKEY splice-sites ..."
hisat2_extract_splice_sites.py $filedir/turkey_genome.gtf > $filedir/turkey_genome.ss &&
echo "TURKEY splice-site extraction completed successfully"

cat $filedir/*.ss > $filedir/host_thev.ss
