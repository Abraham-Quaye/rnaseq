#!/usr/bin/env Rscript

library(DESeq2)
library(apeglm)
library(ashr)
library(magrittr)
library(tximport)
library(AnnotationDbi)
library(Homo.sapiens)
library(readxl)
library(ggrepel)
library(org.Hs.eg.db)
library(pheatmap)
library(ggplotify)
library(patchwork)
library(tidyverse)

# ====================================================================
# TRANSCRIPT ABUNDANCE ANALYSIS WITH DEseq2
# ====================================================================
home_path <- "~/berges_rnaseq"
result_path <- paste0(home_path, "/results")

# 1. Read in experiment metadata (colData) =============================
exp_metadata <- read_xls(paste0(home_path, "/raw_files/metadata.xls"), sheet = "Sample Information",
                         range = "A30:I63") %>%
  dplyr::select(sample_num = "Sample Number", condition = "Sample Name",
         infection = "Group Name") %>%
  mutate(sample_name = paste0("BB", sample_num),
         condition = case_when(str_detect(condition, "^Control") ~
                                   sub("Control ", "mock_", condition),
                                 str_detect(condition, "^WildType") ~
                                   sub("WildType ", "wt_", condition),
                                 TRUE ~ condition),
         condition = str_replace_all(condition, "\\s", "_"),
         infection = case_match(infection,
                                "Control" ~ "mock",
                                "Treatment" ~ "wt",
                                "Treatment2" ~ "r77q"),
         timepoint = if_else(str_detect(condition, "^(wt_|R77Q_)"),
                             str_replace(condition,
                                         "^(wt_|R77Q_)(\\d{1,2}hr)_\\d",
                                         "\\2"),
                             "72hr")) %>%
  mutate(treatment = paste0(infection, timepoint)) %>% 
  base::as.data.frame() %>%
  set_rownames(.$condition) %>%
  dplyr::select(-sample_num) %>%
  mutate(treatment = factor(treatment, levels = base::unique(treatment)))

# 2. Read in salmon quantification files ==========================================
quant_files <- data.frame(files = list.files(list.dirs(paste0(result_path, "/salmon_quant"),
                                                       recursive = F),
                                             pattern = "quant.sf", full.names = T)) %>%
  mutate(sample_name = map_chr(files, ~sub(".+quant_(BB\\d{1,2})/quant\\.sf", "\\1", .x)),
         sample_num = map_dbl(sample_name, ~as.numeric(parse_number(.x)))) %>%
  dplyr::arrange(sample_num) %>% 
  as_tibble() %>% select(sample_name, files) %>%
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
deseq_obj$treatment <- relevel(deseq_obj$treatment, ref = "wt72hr")
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
deseq_results <- tibble(treatments = resultsNames(dds)[-1],
                        dds_results_df = map(treatments,
                                             ~results(dds, name = .x,alpha = 0.05) %>%
                                               as_tibble(., rownames = "gene_id") %>%
                                               dplyr::arrange(padj) %>%
                                               drop_na(padj)),
                        annot_dds_results = map(dds_results_df,
                                                ~left_join(.x, gene_annotations,
                                                           by = "gene_id")),
                        lfc_results = map(treatments, ~lfcShrink(dds, coef = .x,
                                                                 type = "apeglm")),
                        lfc_results_tbl = map(lfc_results,
                                              ~as_tibble(.x, rownames = "gene_id") %>%
                                                dplyr::arrange(padj) %>%
                                                drop_na(padj)),
                        counts_filter = map(treatments, ~get_sample_names(.x)),
                        norm_counts = map2(.x = replicate(length(treatments), dds,
                                                          simplify = F),
                                           .y = counts_filter,
                                           ~DESeq2::counts(.x, normalized = T) %>%
                                             as_tibble(., rownames = "gene_id") %>%
                                             select(gene_id, all_of(.y))),
                        total_res = map2(norm_counts, annot_dds_results,
                                         \(.x, .y) inner_join(.x, .y,
                                                              by = "gene_id") %>%
                                           dplyr::select(-c(baseMean, lfcSE)) %>%
                                           dplyr::select(gene_id, ENTREZID,
                                                         SYMBOL, GENENAME,
                                                         everything())),
                        sig_res = map(total_res, ~filter(.x, padj <= 0.05)),
                        volcano_plt = map2(lfc_results_tbl, treatments,
                            ~plot_volcano(lfc_res_tbl = .x, treatment = .y)),
                        heatmaps = map2(sig_res, treatments, ~plot_topGenes_heatmap(.x, .y)))


# Save Distance matrix for all samples
plot_sample_dists(
  dds = dds, dds_design = as.character(dds@design)[[2]],
  color_grp_feature = as.character(dds@design)[[2]],
  row_labs_feature = "condition", plot_height = 18, plot_width = 18,
  plot_name = paste0(result_path, "/r/figures/r77qProject_distMatrix_all.pdf")
  )


## PCA plot for all samples
ggsave(plot = plot_PCA(
  dds, dds_design = as.character(dds@design)[[2]],
  plot_title = "A PCA of HIV R77Q Vs Wild Type at 4, 8, 12, 24, and 72 Hour Time Points"),
  filename = paste0(result_path, "/r/figures/r77qProject_PCA_all.pdf"),
  height = 8, width = 10)

# save DEG tables =======================================
treatment_labels <- map(deseq_results$treatments,
                        ~str_remove(.x, as.character(dds@design)[[2]]))
# Save Significant DEGs
map2(.x = deseq_results$sig_res, .y = treatment_labels,
     ~write.csv(.x, file = paste0(result_path, "/r/tables/significant", .y, "_DEGs.csv"),
                row.names = F))

# Save Total Results with counts
map2(.x = deseq_results$total_res, .y = treatment_labels,
     ~write.csv(.x, file = paste0(result_path, "/r/tables/total", .y, "_DEGs.csv"),
                row.names = F))

# Save Data QC plots =========================
map2(.x = deseq_results$volcano_plt, .y = treatment_labels,
     ~ggsave(plot = .x,
             filename = paste0(result_path, "/r/figures/volcano", .y, ".pdf"),
             width = 8, height = 8))

map2(.x = deseq_results$heatmaps, .y = treatment_labels,
     ~ggsave(plot = .x,
             filename = paste0(result_path, "/r/figures/heatmap", .y, ".pdf"),
             width = 8, height = 10))
