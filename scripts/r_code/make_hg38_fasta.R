#!/usr/bin/env Rscript --vanilla

library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)

############################ make fasta file ################
hg38 <- Hsapiens

chr_names <- names(hg38)

chr_seqs <- base::lapply(chr_names, FUN = function(chr) hg38[[chr]])

# assign each sequence its corresponding chr name
base::names(chr_seqs) <- chr_names

# make it a DNAStringSet object for indexing
chr_seqs_for_fasta <- DNAStringSet(chr_seqs)

#sava fasta file for indexing
writeXStringSet(x = chr_seqs_for_fasta,
                filepath = "raw_files/genome_files/hg38_genome.fa",
                compress = F, format = "fasta")
