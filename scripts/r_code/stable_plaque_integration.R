
library(DESeq2)
library(magrittr)
library(ggrepel)
library(ggtext)
library(AnnotationDbi)
library(Homo.sapiens)
library(org.Hs.eg.db)
library(pheatmap)
library(ggplotify)
library(apeglm)
library(ashr)
library(ggforce)
library(patchwork)
library(tidyverse)


# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"

path <- paste(urld, "acc=GSE120521",
              "file=GSE120521_raw_counts_GRCh38.p13_NCBI.tsv.gz",
              sep = "&")

tbl <- as.matrix(data.table::fread(path, header = T,
                                   colClasses = "integer"),
                 rownames = 1) %>%
  as_tibble(rownames = "geneid")

exp_metadata <- tibble(plaque = paste0(rep(c("Stable", "Unstable"),
                                           each = 4)),
                       patient = rep(paste0("Pt_", 1:4), 2)) %>%
  arrange(patient) %>%
  mutate(accession = paste0("GSM34025", str_pad(4:11, width = 2,
                                                side = "left",
                                                pad = 0)),
         group = paste0(plaque, "_", patient))

namemap <- setNames(exp_metadata$accession, exp_metadata$group)

exp_metadata <- exp_metadata %>%
  mutate(plaque = factor(plaque)) %>%
  as.data.frame() %>% 
set_rownames(.$group)

count_matrix <- tbl %>% rename(all_of(namemap)) %>%
  column_to_rownames("geneid")

stopifnot(
  all(rownames(exp_metadata) %in% colnames(count_matrix)), # checks identity/presence
  all(rownames(exp_metadata) == colnames(count_matrix)) # checks order
)

deseq_obj <- DESeqDataSetFromMatrix(countData = count_matrix,
                                    colData = exp_metadata,
                                    design = ~plaque)

genes_in_data <- rownames(count_matrix)

gene_annotations <- AnnotationDbi::select(
  Homo.sapiens, keys = genes_in_data,
  columns = c("SYMBOL", "GENENAME",
              "ENSEMBL", "DEFINITION"),
  keytype = "ENTREZID") %>%
  distinct(ENTREZID, .keep_all = T) %>% 
  as_tibble()

# Ensure row alignment perfectly matches deseq_obj before assignment
annotation_ordering <- match(rownames(deseq_obj), gene_annotations$ENTREZID)
sorted_annotations <- gene_annotations[annotation_ordering, ]

mcols(deseq_obj) <- DataFrame(mcols(deseq_obj), sorted_annotations)

deseq_obj$plaque <- relevel(deseq_obj$plaque, ref = "Stable")
dds <- DESeq(deseq_obj)


source("scripts/r_code/DEG_plotting_functions2.R")

# Organize all analysis into one table
deseq_results <- tibble(local_dds = list(dds),
                        treatments = resultsNames(dds)[-1],
                        lfc_results = map2(local_dds, treatments,
                                           ~lfcShrink(.x, coef = .y,
                                                      type = "apeglm")),
                        lfc_results_tbl = map(lfc_results,
                                              ~as_tibble(.x, rownames = "ENTREZID") %>%
                                                arrange(padj) %>%
                                                drop_na(padj) %>%
                                                left_join(., gene_annotations,
                                                          by = "ENTREZID")),
                        norm_counts = map(.x = local_dds,
                                          ~DESeq2::counts(.x, normalized = T) %>%
                                            as_tibble(., rownames = "ENTREZID")),
                        total_res = map2(norm_counts, lfc_results_tbl,
                                         \(.x, .y) inner_join(.x, .y,
                                                              by = "ENTREZID") %>%
                                           select(-c(baseMean, lfcSE)) %>%
                                           arrange(padj)),
                        sig_res = map(total_res, ~filter(.x, padj <= 0.05) %>%
                                        arrange(padj)),
                        volcano_plt = map2(lfc_results_tbl, treatments,
                                           ~plot_volcano(lfc_res_tbl = .x, treatment = .y)),
                        heatmaps = map2(sig_res, treatments,
                                        ~plot_topGenes_heatmap(.x, .y)))

plot_genes_heatmap <- function(results, contr_name){
  
  # Construct the plot title
  feature_remove <- contr_name %>% str_split(., "_") %>%
    flatten() %>% pluck(1)
  
  title_name <- contr_name %>% 
    str_replace_all(., "_", " ") %>%
    str_remove(., feature_remove) %>% 
    toupper(.) %>%
    str_replace_all(., "HR", "hr") %>%
    paste0("Gene Expression for ", .,
           "\nby Adjusted P-values and Log2(Fold Change)")
  
  # Prepare Matrix to plot
  sample_cols <- contr_name %>%
    str_split(., "_") %>%
    base::unlist() %>%
    .[c(2, 4)]
  
  genes <- results %>%
    arrange(padj, log2FoldChange) %>%
    drop_na(SYMBOL) %>%
    mutate(isduplicate = duplicated(SYMBOL)) %>% 
    filter(!isduplicate) %>% 
    select(SYMBOL, starts_with(sample_cols)) %>%
    column_to_rownames("SYMBOL")
  
  genes_mat <-  genes[apply(genes, 1, var) != 0, ] %>% 
    as.matrix()
  
  hmap <- pheatmap(genes_mat,
                   # main = title_name,
                   cluster_rows = TRUE,
                   cluster_cols = F,
                   # treeheight_row = 0,
                   cellwidth = 30,
                   show_rownames = F,
                   show_colnames = TRUE,
                   color = colorRampPalette(
                     colors = c('blue','white','red'))(250),
                   scale = "row")
  
  return(as.ggplot(hmap))
}

plot_genes_heatmap(results = deseq_results$sig_res[[1]], contr_name = deseq_results$treatments)

plaque_sigres <- deseq_results$sig_res[[1]] %>%
  select(human_ENTREZID = ENTREZID, human_ENSEMBL = ENSEMBL, SYMBOL,
         GENENAME, Stable_Pt_1:padj) %>%
  arrange(desc(abs(log2FoldChange)), padj) %>%
  drop_na(SYMBOL)

bulk_data <- read_csv(paste0("~/bulk_slc38a9_rnaseq_aq/results/r/tables/",
                             "significant_12hrKO_vs_12hrWT_DEGs.csv")) %>%
  mutate(SYMBOL = toupper(SYMBOL)) %>%
  drop_na(SYMBOL) %>%
  select(-c(GENENAME, DEFINITION, pvalue)) %>%
  rename(mouse_ENSEMBL = gene_id, mouse_ENTREZID = ENTREZID,
         mouse_log2fc = log2FoldChange, mouse_padj = padj)

common_res <- inner_join(plaque_sigres, bulk_data, by = "SYMBOL")








