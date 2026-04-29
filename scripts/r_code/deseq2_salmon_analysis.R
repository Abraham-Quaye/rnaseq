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
library(patchwork)
library(tidyverse)

# ====================================================================
# TRANSCRIPT ABUNDANCE ANALYSIS WITH DEseq2
# ====================================================================
home_path <- "~/bm_fn_rnaseq"
result_path <- paste0(home_path, "/results/")

# 1. Read in experiment metadata (colData) =============================
exp_metadata <- tibble(sample_name = dir(paste0(result_path,
                                                "salmon_quant")) %>%
                         str_remove(., "quant_"),
                       treatment = ifelse(str_detect(sample_name, "BM"),
                                          "BM", "FN")) %>%
  mutate(treatment = factor(treatment,
                            levels = base::unique(treatment))) %>%
  base::as.data.frame() %>%
  set_rownames(.$sample_name)
  

# 2. Read in salmon quantification files ==========================================
quant_files <- data.frame(
  files = list.files(list.dirs(paste0(result_path, "salmon_quant"),
                     recursive = F), pattern = "quant.sf",
                     full.names = T)
) %>%
  mutate(sample_name = map_chr(files,
                               ~sub(".+quant_((BM|FN)\\d{1,4})/quant\\.sf",
                                    "\\1", .x))) %>%
  left_join(., exp_metadata, by = "sample_name")


quant_files <- pull(quant_files, files) %>%
  set_names(pull(quant_files, sample_name))

# 3. Make a transcript ID to gene ID lookup table ==============================
gtfpath <- "raw_files/annotations/Mus_musculus.GRCm39.gtf.gz"
gtf <- rtracklayer::import(gtfpath) %>%
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
  all(rownames(exp_metadata) %in% colnames(count_matrix$abundance)), # checks identity/presence
  all(rownames(exp_metadata) == colnames(count_matrix$abundance)) # checks order
  )

# 6. Make Deseq dataset object =================
deseq_obj <- DESeqDataSetFromTximport(txi = count_matrix,
                                colData = exp_metadata,
                                design = ~ treatment)

# 6.1 Include additional gene annotations ==============================
genes <- rownames(count_matrix$abundance)

genes <- AnnotationDbi::select(Mus.musculus, keys = genes,
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
deseq_obj$treatment <- relevel(deseq_obj$treatment, ref = "BM")

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
                                                arrange(padj) %>%
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
                                           select(-c(baseMean, lfcSE)) %>%
                                           select(gene_id, ENTREZID,
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

target_receptors <- c("Havcr2", "Calr", "Itgav", "Axl",
                      "Mfge8", "Cd300a", "Lrp1",
                      "Scarb1", "Gas6", "Mertk", "C1qa",
                      "Havcr1", "Cd36", "Itgb5", "Itgb3")

receptors_df <- deseq_results$total_res[[1]] %>% filter(SYMBOL %in% target_receptors) %>%
  arrange(desc(log2FoldChange)) 

receptors_scaled <- receptors_df %>%
  select(SYMBOL, BM1:FN4) %>%
  as.data.frame() %>%
  set_rownames(.$SYMBOL) %>%
  select(-SYMBOL) %>% 
  as.matrix() %>%
  t() %>%
  scale() %>%
  t() %>%
  as.data.frame() %>%
  as_tibble(rownames = "SYMBOL")

receptors_fc <- receptors_df %>% select(SYMBOL,
                                        log2FoldChange,
                                        padj)

receptors_plt_ready <- left_join(receptors_scaled,
                                 receptors_fc,
                                 by = "SYMBOL") %>%
  pivot_longer(BM1:FN4, names_to = "samples",
               values_to = "zscores") %>%
  mutate(sample_grp = ifelse(str_detect(samples, "BM"),
                             "Basement Membrane",
                             "Fibronectin"),
         SYMBOL = factor(SYMBOL,
                         levels = rev(base::unique(.$SYMBOL))),
         fc = 2^log2FoldChange)

plot_receptor_heat <- function(df, sample_num = 4){
  if(sample_num == 3){
    df <- df %>%
      filter(!(samples %in% c("BM4", "FN4")))
    
    text_x1 <- 2
    text_x2 <- 5
    text_y <- 15.2
    seg_x2 <- 3.6
    seg_y <- 14.8
    seg_xend1 <- 3.4
    seg_xend2 <- 6.4
  }else{
    text_x1 <- 2.5
    text_x2 <- 6.5
    text_y <- 15.2
    seg_x2 <- 4.6
    seg_y <- 14.8
    seg_xend1 <- 4.4
    seg_xend2 <- 8.4
  }
    
  plt <- df %>% 
  ggplot(aes(samples, SYMBOL, fill = zscores)) +
  geom_tile(color = "grey40") +
  scale_fill_gradient2(high = "#ff0000",
                       low = "#0000ff",
                       mid = "#ffffff") +
  annotate(geom = "text", x = c(text_x1, text_x2),
           y = text_y,
           label = c("Basement Membrane", " Fibronectin"),
           fontface = "bold", size = 10, size.unit = "pt") +
  annotate(geom = "segment", x = c(0.6, seg_x2),
           xend = c(seg_xend1, seg_xend2),
           y = seg_y, yend = seg_y,
           linewidth = 0.5, color = "#000000") +
  coord_cartesian(clip = "off", expand = FALSE) +
  labs(title = NULL,
       y = NULL, x = NULL, fill = "Z-score") +
  theme(plot.margin = margin_auto(5, 2),
        plot.title.position = "plot",
        plot.title = element_text(margin = margin(b = 20),
                                  hjust = 0.5, size = 14,
                                  face = "bold"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_text(size = 12, color = "#000000",
                                 margin = margin(r = 10, t = 5)),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(0.9, 1),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 10),
        legend.key.height = unit(20, "pt"))
  
  return(plt)
}

rec_plt1a <- plot_receptor_heat(receptors_plt_ready, 3)
rec_plt1b <- plot_receptor_heat(receptors_plt_ready)

ggsave(plot = rec_plt1a,
       filename = paste0(figs_path,
                         "receptors_heatmap1_threeSamples.pdf"),
       width = 5, height = 3.4)

ggsave(plot = rec_plt1b,
       filename = paste0(figs_path,
                         "receptors_heatmap1_allSamples.pdf"),
       width = 5, height = 3.4)

rec_plt2 <- receptors_df %>%
  pivot_longer(BM1:FN4, names_to = "samples",
               values_to = "tpm") %>%
  select(-c(gene_id, log2FoldChange, pvalue, padj)) %>% 
  mutate(sample_grp = ifelse(str_detect(samples, "BM"),
                             "Basement Membrane",
                             "Fibronectin"),
         SYMBOL = factor(SYMBOL, levels = rev(base::unique(.$SYMBOL)))
         ) %>%
  summarise(mean_tpm = round(mean(tpm), 1),
            .by = c(sample_grp, SYMBOL)) %>%
  mutate(col_lab = ifelse(mean_tpm > 50000,
                          "#000000", "#ffffff")) %>% 
  ggplot(aes(sample_grp, SYMBOL, fill = mean_tpm)) +
  geom_tile() +
  geom_text(aes(label = mean_tpm, color = col_lab),
            fontface = "bold",
            size = 12, size.unit = "pt") +
  annotate(geom = "text", x = c(1, 2),
           y = 15.2,
           label = c("Basement Membrane", "Fibronectin"),
           fontface = "bold", size = 10, size.unit = "pt") +
  scale_fill_viridis_c() +
  scale_color_identity() +
  coord_cartesian(clip = "off", expand = FALSE) +
  labs(title = NULL, y = NULL, x = NULL,
       fill = "Mean TPM\n(n = 4)") +
  theme(plot.margin = margin_auto(5, 5),
        plot.title.position = "plot",
        plot.title = element_text(margin = margin(b = 5),
                                  hjust = 0.5, size = 14,
                                  face = "bold"),
        axis.ticks = element_blank(),
        # axis.text.x.top = element_text(size = 10, face = "bold",
        #                                colour = "#000000",
        #                                margin = margin(b = 10)),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12, colour = "#000000"),
        legend.justification = c(0.8, 1),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.height = unit(20, "pt")
  )

ggsave(plot = rec_plt2,
       filename = paste0(figs_path,
                         "receptors_heatmap2.pdf"),
       width = 4, height = 6)

rec_plt2b <- rec_plt2 +
  theme(axis.text.y = element_blank())

rec_composite3 <- (rec_plt1a | rec_plt2b) +
  plot_layout(guides = "collect", axes = "collect_y",
              widths = c(1, 0.75))

ggsave(plot = rec_composite3,
       filename = paste0(figs_path,
                         "receptors_composite_heatmaps_threeSamples.pdf"),
       width = 8.5, height = 4)


rec_composite4 <- (rec_plt1b | rec_plt2b) +
  plot_layout(guides = "collect", axes = "collect_y",
              widths = c(1, 0.7))

ggsave(plot = rec_composite4,
       filename = paste0(figs_path,
                         "receptors_composite_heatmaps_allSamples.pdf"),
       width = 9, height = 4)


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
  row_labs_feature = as.character(dds@design)[[2]]),
  filename = paste0(figs_path, "dists", treatment_labels, ".pdf"),
  height = 4, width = 5)

## PCA plot for all samples
ggsave(plot = plot_PCA(dds = dds, dds_design = as.character(dds@design)[[2]]),
       filename = paste0(figs_path, "pca", treatment_labels, ".pdf"),
       height = 5.5, width = 6.2)
