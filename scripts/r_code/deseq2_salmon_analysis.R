#!/usr/bin/env Rscript

library(DESeq2)
library(apeglm)
library(ashr)
library(magrittr)
library(tximport)
library(AnnotationDbi)
library(Homo.sapiens)
library(ggrepel)
library(org.Hs.eg.db)
library(pheatmap)
library(ggplotify)
library(patchwork)
library(tidyverse)

# ====================================================================
# TRANSCRIPT ABUNDANCE ANALYSIS WITH DEseq2
# ====================================================================
home_path <- "~/myocd_rnaseq"
result_path <- paste0(home_path, "/results/")

# 1. Read in experiment metadata (colData) =============================
exp_metadata <- tibble(sample_name = dir(paste0(result_path, "salmon_quant")) %>%
                         str_remove(., "quant_"),
                       treatment = ifelse(str_detect(sample_name, "GFP"),
                                          "GFP", "MYOCD"),
                       sample_num = str_remove_all(
                         str_extract(sample_name, "_\\d_"), "_"
                         ),
                       condition = paste0(treatment, sample_num)) %>%
  select(-sample_num) %>%
  mutate(treatment = factor(treatment, levels = base::unique(treatment))) %>%
  base::as.data.frame() %>%
  set_rownames(.$condition)
  

# 2. Read in salmon quantification files ==========================================
quant_files <- data.frame(files = list.files(list.dirs(paste0(result_path, "salmon_quant"),
                                                       recursive = F),
                                             pattern = "quant.sf", full.names = T)) %>%
  mutate(sample_name = map_chr(files,
                               ~sub(".+quant_(siO_(GFP|MYOCD)_\\d_S\\d{2})/quant\\.sf",
                                    "\\1", .x))) %>%
  left_join(., exp_metadata, by = "sample_name") %>%
  select(condition, files)

quant_files <- pull(quant_files, files) %>%
  set_names(pull(quant_files, condition))

# 3. Make a transcript ID to gene ID lookup table ==============================
gtf <- rtracklayer::import("raw_files/annotations/Homo_sapiens.GRCh38.115.gtf") %>%
  as.data.frame() 

tx_gene <- gtf %>%
  filter(type == "transcript") %>%
  select(transcript_id, gene_id)

# 3.1  Make a gene_id to gene_name lookup table ==============
gene_name_map <- gtf %>%
  filter(type == "gene") %>%
  select(gene_id, gene_name)

# 4. Load Salmon quantification data and summarize at gene level ====================================
count_matrix <- tximport(files = quant_files, type = "salmon",
                       tx2gene = tx_gene, countsFromAbundance = "lengthScaledTPM",
                       geneIdCol = "gene_id", txIdCol = "transcript_id",
                       ignoreTxVersion = TRUE)

# 5. Check matching samples in experimental metadata and count matrix ==========
stopifnot(
  all(rownames(exp_metadata) %in% colnames(count_matrix$abundance)), # checks identity/presence
  all(rownames(exp_metadata) == colnames(count_matrix$abundance)) # checks order
  )

# 6. Make Deseq dataset object =================
deseq_obj <- DESeqDataSetFromTximport(txi = count_matrix,
                                colData = exp_metadata,
                                design = ~ treatment)

# 6.1 Include additional gene annotations ==============================
genes <- rownames(count_matrix$abundance)

genes <- AnnotationDbi::select(Homo.sapiens, keys = genes,
                               columns = c('SYMBOL','GENENAME', "ENTREZID"),
                               keytype = 'ENSEMBL') %>%
  distinct(ENSEMBL, .keep_all = T) %>% 
  as_tibble()

gene_annotations <- left_join(gene_name_map, genes,
                              by = c("gene_id" = "ENSEMBL")) %>%
  as_tibble() %>%
  mutate(SYMBOL = ifelse(is.na(SYMBOL) & !is.na(gene_name),
                         gene_name, SYMBOL)) %>%
  select(gene_id, SYMBOL, GENENAME, ENTREZID) %>%
  filter(gene_id %in% rownames(deseq_obj))

mcols(deseq_obj) <- DataFrame(mcols(deseq_obj), gene_annotations)

# 6.2 Change reference sample for comparisons ==========================
deseq_obj$treatment <- relevel(deseq_obj$treatment, ref = "GFP")

dds <- DESeq(deseq_obj)
resultsNames(dds)

# results from lfcShrink(dds) is better:
# Stabilizes LFC estimates: Raw LFCs from results() can be noisy, especially for low-count genes. Shrinking reduces exaggerated fold changes.
# Improves ranking: Helps prioritize genes by effect size rather than just statistical significance.
# Better visualization: Shrunken LFCs look cleaner in volcano plots or MA plots.
# More conservative estimates: Useful when reporting fold changes in publications or downstream analyses.

# DESeq2::counts:
# it gives you: Normalized raw counts using size factors to account for differences in sequencing depth across samples.
# Values: Still on the count scale (integers or decimals), not log-transformed.
# Use case: Good for downstream statistical modeling, but not ideal for visualization due to high variance and skew.

source("scripts/r_code/DEG_plotting_functions.R")

# Organize all analysis into one table
deseq_results <- tibble(local_dds = list(dds),
                        treatments = resultsNames(dds)[-1],
                        lfc_results = map2(local_dds, treatments,
                                          ~lfcShrink(.x, coef = .y,
                                                     type = "apeglm")),
                        lfc_results_tbl = map(lfc_results,
                                              ~as_tibble(.x, rownames = "gene_id") %>%
                                                dplyr::arrange(padj) %>%
                                                drop_na(padj) %>%
                                                left_join(., gene_annotations,
                                                          by = "gene_id") %>%
                                                select(gene_id, SYMBOL:ENTREZID,
                                                       everything())),
                        norm_counts = map(.x = local_dds,
                                           ~DESeq2::counts(.x, normalized = T) %>%
                                             as_tibble(., rownames = "gene_id")),
                        total_res = map2(norm_counts, lfc_results_tbl,
                                         \(.x, .y) inner_join(.x, .y,
                                                              by = "gene_id") %>%
                                           dplyr::select(-c(baseMean, lfcSE)) %>%
                                           dplyr::select(gene_id, ENTREZID,
                                                         SYMBOL, GENENAME,
                                                         everything()) %>%
                                           arrange(padj)),
                        sig_res = map(total_res, ~filter(.x, padj <= 0.05) %>%
                                        arrange(padj)),
                        volcano_plt = map2(lfc_results_tbl, treatments,
                            ~plot_volcano(lfc_res_tbl = .x, treatment = .y)),
                        heatmaps = map2(sig_res, treatments,
                                        ~plot_topGenes_heatmap(.x, .y)))

figs_path <- paste0(result_path, "r/figures/")
if(!dir.exists(figs_path)){dir.create(figs_path, recursive = T)}

tables_path <- paste0(result_path, "r/tables/")
if(!dir.exists(tables_path)){dir.create(tables_path, recursive = T)}

# save DEG tables =======================================
treatment_labels <- map(deseq_results$treatments,
                        ~str_remove(.x, as.character(dds@design)[[2]]))
# Save Significant DEGs
map2(.x = deseq_results$sig_res, .y = treatment_labels,
     ~write.csv(.x, file = paste0(tables_path, "significant", .y, "_DEGs.csv"),
                row.names = F))

# Save Total Results with counts
map2(.x = deseq_results$total_res, .y = treatment_labels,
     ~write.csv(.x, file = paste0(tables_path, "total", .y, "_DEGs.csv"),
                row.names = F))

# Save Data QC plots =========================
map2(.x = deseq_results$volcano_plt, .y = treatment_labels,
     ~ggsave(plot = .x,
             filename = paste0(figs_path, "volcano", .y, ".pdf"),
             width = 8, height = 8))

map2(.x = deseq_results$heatmaps, .y = treatment_labels,
     ~ggsave(plot = .x,
             filename = paste0(result_path, "r/figures/heatmap", .y, ".pdf"),
             width = 6.5, height = 8))

# Save Distance matrix for all samples
ggsave(plot = plot_sample_dists(
  dds = dds, dds_design = as.character(dds@design)[[2]],
  color_grp_feature = as.character(dds@design)[[2]],
  row_labs_feature = "condition"),
  filename = paste0(figs_path, "dists", treatment_labels, ".pdf"),
  height = 4, width = 5)

## PCA plot for all samples
ggsave(plot = plot_PCA(dds = dds, dds_design = as.character(dds@design)[[2]]),
       filename = paste0(figs_path, "pca", treatment_labels, ".pdf"),
       height = 5.5, width = 6.2)
