#!/usr/bin/env Rscript --vanilla

# library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
# library(org.Hs.eg.db)
library(magrittr)
library(Rsubread)
library(rtracklayer)

hg38_trxpts <- TxDb.Hsapiens.UCSC.hg38.knownGene

# create GTF file
# ensb_genes <- genes(EnsDb.Hsapiens.v86)
# rtracklayer::export(ensb_genes, format = "gtf",
#                     con = "results/annotations/ens_hg38_known_genes.gtf")

# Extract all exon info and make SAF format file for aligner to know annotated splice sites
exon_sites_SAF <- exons(hg38_trxpts, columns = c("gene_id", "tx_id")) %>%
    # remove exons with unannotated gene_ids (with gene_id length of 0)
    subset(., lengths(gene_id) == 1) %>%
    as.data.frame() %>%
    dplyr::select(GeneID = gene_id, Chr = seqnames,
                  Start = start, End = end, Strand = strand)




# Begin mapping trimmed fastq reads to hg38 genomic index

map_reads <- function(for_reads, rev_reads){

    # bam_name <- unlist(strsplit(for_reads, split = "[_]"))[c(1, 2)] %>%
    #     paste0(., collapse = "_")

    mapped <- subjunc(index = "results/genome_index/hg38_ref",
            readfile1 = "results/trimmed_reads/LCS9697_BB1_Clean_Data1_val_1.fq.gz",#for_reads,
            readfile2 = "results/trimmed_reads/LCS9697_BB1_Clean_Data2_val_2.fq.gz",#rev_reads,
            nthreads = 12,
            useAnnotation = T,
            unique = T,
            annot.ext = exon_sites_SAF,
            isGTF = F, output_format = "BAM"#,
            # output_file = paste0(bam_name, ".bam")
            )
}
