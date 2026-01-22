#!/usr/bin/env Rscript --vanilla

library(Rsubread)

buildindex(basename = "raw_files/rsubread_GRCh38_genome_index/grch38_ref",
           reference = "raw_files/genome_files/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
           memory = 45000, indexSplit = T, gappedIndex = F)

           
buildindex(basename = "raw_files/rsubread_hg38_genome_index/hg38_ref",
           reference = "raw_files/genome_files/hg38_genome.fa",
           memory = 45000, indexSplit = T, gappedIndex = F)
