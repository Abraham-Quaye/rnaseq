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
  design = ~ condition)

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

# Run DESeq2 Analysis

source("scripts/r_code/DEG_plotting_functions2.R")

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
         lfc_results_tbl = map(treatments,
                               ~lfcShrink(dds = dds_local, coef = .x,
                                          type = "apeglm") %>%
                             as_tibble(rownames = "gene_id") %>%
                             left_join(., gene_annotations, by = "gene_id") %>%
                               arrange(padj) %>% drop_na(padj)),
         counts_filter = map(treatments, ~get_sample_names(.x)),
         norm_counts = map(.x = counts_filter,
                           ~DESeq2::counts(dds_local, normalized = T) %>%
                             as_tibble(., rownames = "gene_id") %>%
                             select(gene_id, all_of(.x))),
         total_res = map2(norm_counts, annot_dds_results,
                          \(.x, .y) inner_join(.x, .y, by = "gene_id") %>%
                            dplyr::select(-c(baseMean, lfcSE)) %>%
                            dplyr::select(gene_id, ENTREZID, SYMBOL,
                                          GENENAME, DEFINITION, everything())),
         sig_res = map(total_res, ~filter(.x, padj <= 0.05 &
                                            abs(log2FoldChange) >= 1)),
         volcano_plt = map2(lfc_results_tbl, treatments,
                            ~plot_volcano(lfc_res_tbl = .x, treatment = .y)),
         heatmaps = map2(sig_res, treatments, ~plot_topGenes_heatmap(.x, .y)),
         treatment_labels = map(.x = treatments,
                                ~str_remove(.x, as.character(dds_local@design)[[2]]))
  )
}

full_deseq <- tibble(
  comparisons = c(
    paste0(c("NoACKO", "12hrWT", "12hrKO", "24hrWT", "24hrKO"), "_vs_NoACWT"),
    "12hrKO_vs_12hrWT", "24hrKO_vs_24hrWT"),
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
  dds = map2(
    .x = txi_matrix, .y = exp_meta,
    \(count_mat = .x, metadata = .y){
      
      dobj <- DESeqDataSetFromTximport(txi = count_mat,
                                       colData = metadata,
                                       design = ~ condition)
      
      dobj <- dobj[rowSums(counts(dobj,normalized = FALSE) >= 10) >= 3, ]
      
      if(sum(str_detect(metadata$condition, "NoACWT")) == 3){
        print(paste0("Samples are: ",
                     paste(rownames(metadata), collapse = ", "),
                     "\n Reference sample is set to \"NoACWT\""))
        dobj$condition <- relevel(dobj$condition, ref = "NoACWT")
      }else if(sum(str_detect(metadata$condition, "12hrWT")) == 3){
        
        print(paste0("Samples are: ",
                     paste(rownames(metadata), collapse = ", "),
                     "\n Reference sample is set to \"12hrWT\""))
        
        dobj$condition <- relevel(dobj$condition, ref = "12hrWT")
      }else if(sum(str_detect(metadata$condition, "24hrWT")) == 3){
        print(paste0("Samples are: ",
                     paste(rownames(metadata), collapse = ", "),
                     "\n Reference sample is set to \"24hrWT\""))
        
        dobj$condition <- relevel(dobj$condition, ref = "24hrWT")
      }else{
        stop("NO MATCHING REFERENCE SAMPLE FOUND")
      }
      
      dds <- DESeq(dobj)
      
      return(dds)
      }),
  dds_contrasts = map_chr(dds, ~resultsNames(.x)[-1]),
  deseq_results = map2(.x = dds, .y = dds_contrasts,
                       ~process_dds_results(dds_local = .x,
                                            treatments = .y)),
  sample_corr_plt = map(.x = dds,
                        ~plot_sample_dists(
                          dds = .x,
                          color_grp_feature = as.character(.x@design)[[2]],
                          row_labs_feature = "condition")),
  pca_plt = map(.x = dds,
                ~plot_PCA(dds = .x, dds_design = as.character(.x@design)[[2]]))
)

# Save data ----
figs_path <- paste0(result_path, "r/figures2/")
if(!dir.exists(figs_path)){dir.create(figs_path, recursive = T)}

tables_path <- paste0(result_path, "r/tables2/")
if(!dir.exists(tables_path)){dir.create(tables_path, recursive = T)}

map(.x = full_deseq[["deseq_results"]],
    \(res = .x){
      map2(.x = res$sig_res,
           .y = str_remove(res$treatment_labels, "condition"),
           ~write.csv(.x, file = paste0(
             tables_path, "significant", .y, "_DEGs.csv"),
             row.names = F))
    })

# Save Total Results with counts
map(.x = full_deseq[["deseq_results"]],
    \(res = .x){
      map2(.x = res$total_res,
           .y = str_remove(res$treatment_labels, "condition"),
           ~write.csv(.x, file = paste0(
             tables_path, "total", .y, "_DEGs.csv"),
             row.names = F))
    })

# Save Data QC plots
map(.x = full_deseq[["deseq_results"]],
    \(res = .x){
      map2(.x = res$volcano_plt,
           .y = str_remove(res$treatment_labels, "condition"),
           ~ggsave(filename = paste0(figs_path, "volcano", .y, ".png"),
                   plot = .x, width = 6.5, height = 6.5, dpi = 350))
    })

map(.x = full_deseq[["deseq_results"]],
    \(res = .x){
      map2(.x = res$heatmaps,
           .y = str_remove(res$treatment_labels, "condition"),
           ~ggsave(filename = paste0(figs_path, "heatmap", .y, ".png"),
                   plot = .x, width = 6, height = 7.5, dpi = 350))
    })

map2(.x = full_deseq$sample_corr_plt, .y = full_deseq$comparisons,
     ~ggsave(filename = paste0(figs_path, "sample_dists_", .y, ".png"),
             plot = .x, width = 6.5, height = 8, dpi = 350))

map2(.x = full_deseq$pca_plt, .y = full_deseq$comparisons,
     ~ggsave(filename = paste0(figs_path, "PCA_", .y, ".png"),
             plot = .x, width = 6.5, height = 6.5, dpi = 350))
