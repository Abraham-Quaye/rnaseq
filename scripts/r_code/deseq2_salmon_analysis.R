#!/usr/bin/env Rscript

library(DESeq2)
library(apeglm)
library(ashr)
library(magrittr)
library(tximport)
library(AnnotationDbi)
library(Mus.musculus)
library(ggrepel)
library(org.Mm.eg.db)
library(pheatmap)
library(ggplotify)
library(ggforce)
library(patchwork)
library(tidyverse)

# ====================================================================
# TRANSCRIPT ABUNDANCE ANALYSIS WITH DEseq2
# ====================================================================
home_path <- "~/bulk_slc38a9_rnaseq_aq"
result_path <- paste0(home_path, "/results/")

# 1. Read in experiment metadata (colData) =============================
exp_metadata <- tibble(sample_name = dir(paste0(result_path,
                                                "salmon_quant")) %>%
                         str_remove(., "quant_"),
                       genotype = ifelse(str_detect(sample_name, "WT"),
                                          "WT", "KO"),
                       timepoint = str_replace(sample_name,
                                               "(\\w+)(WT|KO)\\d",
                                               "\\1")) %>%
  mutate(genotype = factor(genotype, levels = c("WT", "KO")),
         timepoint = factor(timepoint),
         condition = factor(base::paste0(timepoint, genotype))) %>%
  column_to_rownames("sample_name")

# 2. Read in salmon quantification files ==========================================
quant_files <- data.frame(
  files = list.files(list.dirs(paste0(result_path, "salmon_quant"),
                     recursive = F), pattern = "quant.sf",
                     full.names = T)) %>%
  mutate(sample_name = map_chr(files,
                               ~sub(".+quant_(\\w+)/quant\\.sf",
                                    "\\1", .x)))

quant_files <- pull(quant_files, files) %>%
  set_names(pull(quant_files, sample_name))

# 3. Make a transcript ID to gene ID lookup table ==============================
gtf <- rtracklayer::import(
  "raw_files/annotations/Mus_musculus.GRCm39.gtf.gz") %>%
  as.data.frame() 

tx_gene <- gtf %>%
  filter(type == "transcript") %>%
  select(transcript_id, gene_id)

# 3.1  Make a gene_id to gene_name lookup table ==============
gene_name_map <- gtf %>%
  filter(type == "gene") %>%
  select(gene_id, gene_name)

# 4. Load Salmon quantification data and summarize at gene level ====================================
count_matrix <- tximport(files = quant_files,
                         type = "salmon", tx2gene = tx_gene,
                         countsFromAbundance = "lengthScaledTPM",
                         geneIdCol = "gene_id",
                         txIdCol = "transcript_id",
                         ignoreTxVersion = TRUE)

# 5. Check matching samples in experimental metadata and count matrix ==========
stopifnot(
  all(rownames(exp_metadata) %in% colnames(count_matrix$counts)), # checks identity/presence
  all(rownames(exp_metadata) == colnames(count_matrix$counts)) # checks order
  )

# 6. Make Deseq dataset object =================
deseq_obj <- DESeqDataSetFromTximport(
  txi = count_matrix,
  colData = exp_metadata,
  design = ~ genotype + timepoint + genotype:timepoint)

# Find your smallest group size (e.g., 3 replicates)
# Keep genes with a count >= 10 in at least 3 samples
deseq_obj <- deseq_obj[rowSums(counts(deseq_obj, normalized = FALSE) >= 10) >= 3, ]

# 6.1 Include additional gene annotations ==============================
genes_in_data <- rownames(count_matrix$counts)

genes <- AnnotationDbi::select(
  Mus.musculus, keys = genes_in_data,
  columns = c("SYMBOL", "GENENAME",
              "ENTREZID", "DEFINITION"),
  keytype = "ENSEMBL") %>%
  distinct(ENSEMBL, .keep_all = T) %>% 
  as_tibble()

gene_annotations <- left_join(
  gene_name_map, genes, by = c("gene_id" = "ENSEMBL")) %>%
  as_tibble() %>%
  mutate(SYMBOL = ifelse(is.na(SYMBOL) & !is.na(gene_name),
                         gene_name, SYMBOL)) %>%
  select(gene_id, ENTREZID, SYMBOL, GENENAME, DEFINITION)

# Ensure row alignment perfectly matches deseq_obj before assignment
annotation_ordering <- match(rownames(deseq_obj), gene_annotations$gene_id)
sorted_annotations <- gene_annotations[annotation_ordering, ]

mcols(deseq_obj) <- DataFrame(mcols(deseq_obj), sorted_annotations)

# 6.2 Change reference sample for comparisons ==========================
# design_name <- as.character(design(deseq_obj))[[2]]
# 
# deseq_obj[[design_name]] <- relevel(deseq_obj[[design_name]],
#                                     ref = "WT_0hr")
# dds <- DESeq(deseq_obj)

deseq_obj$genotype <- relevel(deseq_obj$genotype, ref = "WT")
deseq_obj$timepoint <- relevel(deseq_obj$timepoint, ref = "NoAC")

dds <- DESeq(deseq_obj, test = "LRT", reduced = ~ genotype + timepoint)

normalized_counts <- DESeq2::counts(dds, normalized = TRUE) %>%
  as_tibble(., rownames = "gene_id")

source("scripts/r_code/DEG_plotting_functions.R")

get_sample_names <- function(contr_name, metadata){
  treatments <- str_split(contr_name, "_") %>%
    base::unlist(.) %>% .[c(1,3)]
  
  metadata %>% filter(condition %in% treatments) %>%
    rownames()
}

deseq_results <- tibble(
  comparisons = c(
    paste0(c("NoACKO", "12hrWT", "12hrKO", "24hrWT", "24hrKO"), "_vs_NoACWT"),
    "12hrKO_vs_12hrWT", "24hrKO_vs_24hrWT"),
  dds_contrasts = setNames(list("genotype_KO_vs_WT", "timepoint_12hr_vs_NoAC",
                       c("genotype_KO_vs_WT", "timepoint_12hr_vs_NoAC",
                         "genotypeKO.timepoint12hr"),
                       "timepoint_24hr_vs_NoAC",
                       c("genotype_KO_vs_WT", "timepoint_24hr_vs_NoAC",
                         "genotypeKO.timepoint24hr"),
                       c("genotype_KO_vs_WT", "genotypeKO.timepoint12hr"),
                       c("genotype_KO_vs_WT", "genotypeKO.timepoint24hr")),
                       nm = comparisons),
  dds_results = map(.x = dds_contrasts,
                      ~results(dds, contrast = list(as.character(.x)),
                               test = "Wald", alpha = 0.05)),
  lfc_results_tbl = map2(.x = dds_contrasts, .y = dds_results,
                     ~lfcShrink(dds, contrast = list(as.character(.x)),
                                res = .y, type = "ashr",
                                saveCols = c("ENTREZID", "SYMBOL",
                                             "GENENAME", "DEFINITION")) %>%
                       as_tibble(rownames = "gene_id") %>%
                       arrange(padj) %>% drop_na(padj)),
  counts_filter = map(names(dds_contrasts),
                          ~get_sample_names(.x, metadata = exp_metadata)),
  norm_counts = map(counts_filter, \(sample_cols){
    normalized_counts %>% select(gene_id, all_of(sample_cols))
    }),
  total_res = map2(norm_counts, lfc_results_tbl,
                   \(.x, .y) left_join(.x, .y, by = "gene_id") %>%
                     select(-c(baseMean, lfcSE)) %>%
                     select(gene_id, ENTREZID, SYMBOL,
                            GENENAME, DEFINITION, everything()) %>%
                     arrange(padj)),
  sig_res = map(total_res, ~filter(.x, padj <= 0.05) %>% arrange(padj)),
  volcano_plt = map2(lfc_results_tbl, names(dds_contrasts),
                     ~plot_volcano(lfc_res_tbl = .x, treatment = .y)),
  heatmaps = map2(sig_res, names(dds_contrasts),
                  ~plot_topGenes_heatmap(.x, .y)))

figs_path <- paste0(result_path, "r/figures/")
if(!dir.exists(figs_path)){dir.create(figs_path, recursive = T)}

tables_path <- paste0(result_path, "r/tables/")
if(!dir.exists(tables_path)){dir.create(tables_path, recursive = T)}

# save DEG tables =======================================
treatment_labels <- names(deseq_results$dds_contrasts)

# Save Significant DEGs
map2(.x = deseq_results$sig_res, .y = treatment_labels,
     ~write.csv(.x, file = paste0(tables_path, "significant_", .y, "_DEGs.csv"),
                row.names = F))

# Save Total Results with counts
map2(.x = deseq_results$total_res, .y = treatment_labels,
     ~write.csv(.x, file = paste0(tables_path, "total_", .y, "_DEGs.csv"),
                row.names = F))

# Save Data QC plots =========================
map2(.x = deseq_results$volcano_plt, .y = treatment_labels,
     ~ggsave(plot = .x,
             filename = paste0(figs_path, "volcano_", .y, ".pdf"),
             width = 8, height = 8))

map2(.x = deseq_results$heatmaps, .y = treatment_labels,
     ~ggsave(plot = .x,
             filename = paste0(figs_path, "heatmap_", .y, ".pdf"),
             width = 6.5, height = 8))

# Save Distance matrix for all samples
ggsave(plot = plot_sample_dists(
  dds = dds, color_grp_feature = "condition",
  row_labs_feature = "condition") + theme(plot.margin = margin_auto(10, 10)),
  filename = paste0(figs_path, "sample_distance_heatmap.pdf"),
  height = 8, width = 8)

## PCA plot for all samples
ggsave(plot = plot_PCA(dds, dds_design = "condition"),
       filename =  paste0(figs_path, "sample_PCA_complete.pdf"),
       height = 7, width = 8.7)
