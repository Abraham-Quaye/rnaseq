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
# library(ggforce)
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
                                   sub("Control", "mock", condition),
                                 str_detect(condition, "^WildType") ~
                                   sub("WildType ", "wt", condition),
                               str_detect(condition, "^R77Q") ~
                                 sub("R77Q ", "r77q", condition),
                                 TRUE ~ condition),
         condition = str_replace_all(condition, "\\s", "_"),
         infection = case_match(infection,
                                "Control" ~ "mock",
                                "Treatment" ~ "wt",
                                "Treatment2" ~ "r77q"),
         timepoint = if_else(str_detect(condition, "^(wt|r77q)"),
                             str_replace(condition,
                                         "^(wt|r77q)(\\d{1,2}hr)_\\d",
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

# 4. Load Salmon quantification data and summarize at gene level ====================
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


# 7. Run DE analysis
source("scripts/r_code/DEG_plotting_functions.R")
# results from lfcShrink(dds) is better:
# Stabilizes LFC estimates: Raw LFCs from results() can be noisy, especially for low-count genes. Shrinking reduces exaggerated fold changes.
# Improves ranking: Helps prioritize genes by effect size rather than just statistical significance.
# Better visualization: Shrunken LFCs look cleaner in volcano plots or MA plots.
# More conservative estimates: Useful when reporting fold changes in publications or downstream analyses.

# DESeq2::counts:
# it gives you: Normalized raw counts using size factors to account for differences in sequencing depth across samples.
# Values: Still on the count scale (integers or decimals), not log-transformed.
# Use case: Good for downstream statistical modeling, but not ideal for visualization due to high variance and skew.

process_dds_results <- function(dds_local, treatments){
  tibble(treatments = treatments,
         dds = replicate(n = NROW(treatments), expr = dds_local, simplify = F),
         dds_results_df = map(treatments,
                              ~results(dds_local, name = .x, alpha = 0.05) %>%
                                as_tibble(., rownames = "gene_id") %>%
                                dplyr::arrange(padj) %>%
                                drop_na(padj)),
         annot_dds_results = map(dds_results_df,
                                 ~left_join(.x, gene_annotations, by = "gene_id")),
         lfc_results = map(treatments, ~lfcShrink(dds = dds_local, coef = .x,
                                                  type = "apeglm")),
         lfc_results_tbl = map(lfc_results,
                               ~as_tibble(.x, rownames = "gene_id") %>%
                                 dplyr::arrange(padj) %>%
                                 drop_na(padj)),
         counts_filter = map(treatments, ~get_sample_names(.x)),
         norm_counts = map(.x = counts_filter,
                            ~DESeq2::counts(dds_local, normalized = T) %>%
                              as_tibble(., rownames = "gene_id") %>%
                              select(gene_id, all_of(.x))),
         total_res = map2(norm_counts, annot_dds_results,
                          \(.x, .y) inner_join(.x, .y, by = "gene_id") %>%
                            dplyr::select(-c(baseMean, lfcSE)) %>%
                            dplyr::select(gene_id, ENTREZID, SYMBOL,
                                          GENENAME, everything())),
         sig_res = map(total_res, ~filter(.x, padj <= 0.05)),
         volcano_plt = map2(lfc_results_tbl, treatments,
                            ~plot_volcano(lfc_res_tbl = .x, treatment = .y)),
         heatmaps = map2(sig_res, treatments, ~plot_topGenes_heatmap(.x, .y)),
         treatment_labels = map(.x = treatments,
                                ~str_remove(.x, as.character(dds_local@design)[[2]]))
         )
}

full_deseq <- tibble(
  comparisons = c("all_vs_mock72hr", "r77q4hr_vs_wt4hr", "r77q8hr_vs_wt8hr",
                  "r77q12hr_vs_wt12hr", "r77q24hr_vs_wt24hr",
                  "r77q72hr_vs_wt72hr"),
  ref_treatment = map_chr(comparisons, ~str_remove(.x, "\\w+vs_")),
  quant_file = map(comparisons,
                   ~get_quant_files(quant_files, .x)),
  exp_meta = map(comparisons,
                 ~get_metadata(exp_metadata, .x)),
  txi_matrix = map(quant_file,
                   ~tximport(files = .x, type = "salmon",
                             tx2gene = tx_gene, countsFromAbundance = "lengthScaledTPM",
                             geneIdCol = "gene_id", txIdCol = "transcript_id",
                             ignoreTxVersion = TRUE)),
  deseq_obj = map2(.x = txi_matrix, .y = exp_meta,
                   ~DESeqDataSetFromTximport(txi = .x,
                                             colData = .y,
                                             design = ~ treatment)),
  dds = map(deseq_obj, ~DESeq(.x)),
  dds_contrasts = map(dds, ~resultsNames(.x)[-1]),
  deseq_results = map2(.x = dds, .y = dds_contrasts,
                       ~process_dds_results(dds_local = .x,
                                            treatments = .y)),
  sample_corr_plt = map(.x = dds,
                        ~plot_sample_dists(
                          dds = .x,
                          dds_design = as.character(.x@design)[[2]],
                          color_grp_feature = as.character(.x@design)[[2]],
                          row_labs_feature = "condition")),
  pca_plt = map(.x = dds,
                ~plot_PCA(dds = .x, dds_design = as.character(.x@design)[[2]]))
)

# save DEG tables =======================================
# Save Significant DEGs
map(.x = full_deseq[["deseq_results"]],
    \(.x){
      map2(.x = .x$sig_res, .y = .x$treatment_labels,
           ~write.csv(
             .x,
             file = paste0(result_path,
                           "/r/tables/significant",
                           .y, "_DEGs.csv"),
                row.names = F))
    })


# Save Total Results with counts
map(.x = full_deseq[["deseq_results"]],
    \(.x){
      map2(.x = .x$total_res, .y = .x$treatment_labels,
           ~write.csv(
             .x,
             file = paste0(result_path,
                           "/r/tables/total",
                           .y, "_DEGs.csv"),
                row.names = F))
    })

# Save Data QC plots =========================
# save volcano plots
map(.x = full_deseq[["deseq_results"]],
    \(.x){
      map2(.x = .x$volcano_plt, .y = .x$treatment_labels,
           \(.x, .y){
             if(is.null(.x)){return(NULL)}
             ggsave(plot = .x,
             filename = paste0(result_path,
                               "/r/figures/volcano",
                               .y, ".pdf"),
             width = 8, height = 8)})
    })

# Save heatmaps
map(.x = full_deseq[["deseq_results"]],
    \(.x){
      map2(.x = .x$heatmaps, .y = .x$treatment_labels,
           \(.x, .y){
             if(is.null(.x)){return(NULL)}
             ggsave(plot = .x,
             filename = paste0(result_path,
                               "/r/figures/heatmap",
                               .y, ".pdf"),
             width = 8, height = 10)})
    })

# save PCA plots
map2(.x = full_deseq$pca_plt, .y = full_deseq$comparisons,
    ~ggsave(plot = .x,
            filename = paste0(result_path,
                               "/r/figures/PCA_",
                               .y, ".pdf"),
            width = 10, height = 8)
)

# save sample correlation plots
map2(.x = full_deseq$sample_corr_plt, .y = full_deseq$comparisons,
    ~ggsave(plot = .x,
            filename = paste0(result_path,
                               "/r/figures/sampleCorrelation_",
                               .y, ".pdf"))
)



