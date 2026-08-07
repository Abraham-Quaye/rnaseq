
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
library(ComplexHeatmap)
library(clusterProfiler)
library(enrichplot)
library(pathview)
library(fgsea)
library(ggarchery)
library(ggtangle)
library(GOplot)
library(circlize)
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
  
  genes_mat <- genes[apply(genes, 1, var) != 0, ] %>% 
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

plot_genes_heatmap(results = deseq_results$sig_res[[1]],
                   contr_name = deseq_results$treatments)

plaque_sigres <- deseq_results$sig_res[[1]] %>%
  select(human_ENTREZID = ENTREZID, human_ENSEMBL = ENSEMBL, SYMBOL,
         GENENAME, Stable_Pt_1:padj) %>%
  arrange(desc(abs(log2FoldChange)), padj) %>%
  drop_na(SYMBOL)

bulk_data <- read_csv(paste0("~/bulk_slc38a9_rnaseq_aq/results/r/tables/",
                             "significant_24hrKO_vs_24hrWT_DEGs.csv")) %>%
  mutate(SYMBOL = toupper(SYMBOL)) %>%
  drop_na(SYMBOL) %>%
  select(-c(GENENAME, DEFINITION, pvalue)) %>%
  rename(mouse_ENSEMBL = gene_id, mouse_ENTREZID = ENTREZID,
         mouse_log2fc = log2FoldChange, mouse_padj = padj)

common_res <- inner_join(plaque_sigres, bulk_data, by = "SYMBOL") 

# FUNCTIONAL ENRICHMENT ANALYSES =======================================

bggenes <- deseq_results$total_res[[1]] %>% drop_na(ENTREZID) %>% pull(ENTREZID)

common_dges <- common_res %>% select(ENTREZID = human_ENTREZID,
                                     SYMBOL, Stable_Pt_1:padj) %>%
  mutate(significant = padj <= 0.05,
         reg = case_when(significant & log2FoldChange > 0 ~ "up",
                         significant & log2FoldChange < 0 ~ "down",
                         TRUE ~ "not sig")) %>%
  select(-contains("table_"))


plot_dotplot <- function(res, labb){
  if(is.null(res) | nrow(as_tibble(res)) == 0){return(NULL)
  }else if(nrow(res) > 25){
    num_cat <- 25
  }else{
    num_cat <- nrow(res)
  }
  
  dotplot(object = res, showCategory = num_cat,
          title = paste0("Significant KEGG pathways for: ",
                         labb, "DEGs"),
          font.size = 10.5) +
    theme(plot.title = element_text(face = "bold",
                                    size = 15,
                                    hjust = 0.5))
}

make_pretty <- function(name_){
  
  words <- unlist(strsplit(as.character(name_), split = " "))
  
  if(length(words) > 5){
    labb <- str_replace(
      name_,
      "([-+,A-Za-z]+\\s[-+,A-Za-z]+)\\s([-+,A-Za-z]+\\s[-+,A-Za-z]+)\\s(.+)",
      "\\1\n\\2\n\\3")
    return(labb)
  }else if(length(words) > 3){
    
    labb <- str_replace(name_,
                        "([-+,A-Za-z]+\\s[-+,A-Za-z]+)\\s(.+)",
                        "\\1\n\\2")
    return(labb)
  }else if(length(words) >= 2 & any(nchar(words) >= 15)){
    
    labb <- str_replace(name_,
                        "([-+,A-Za-z]+)\\s(.+)",
                        "\\1\n\\2")
    return(labb)
    
  }else{return(name_)}
}

# Perform functional enrichment analyses
perform_fea <- function(sigdata_, bggenes_){
  
  enrich_result <- tibble(
    regulation = c("total_sig", "up", "down"),
    sig_res = list(sigdata_,
                   sigdata_ %>% filter(reg == "up"),
                   sigdata_ %>% filter(reg == "down")),
    kegg_res = map(sig_res,
                   ~enrichKEGG(
                     gene = .x$ENTREZID,
                     qvalueCutoff = 0.05,
                     pvalueCutoff = 0.05,
                     # universe = bggenes_,
                     organism = "hsa")),
    kegg_dotplot = map2(kegg_res, regulation,
                        ~plot_dotplot(res = .x, labb = .y)),
    go_res = map(sig_res, ~enrichGO(gene = .x$ENTREZID,
                                    keyType = "ENTREZID",
                                    OrgDb = org.Hs.eg.db,
                                    # universe = bggenes_,
                                    ont = "ALL",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = 0.05,
                                    qvalueCutoff = 0.05,
                                    readable = T)),
    go_dotplot = map2(go_res, regulation,
                      ~plot_dotplot(res = .x, labb = .y))
  )
  
  return(enrich_result)
}

common_enrich_result <- perform_fea(sigdata_ = common_dges,
                                    bggenes_ = bggenes)

total_kegg_res <- common_enrich_result$kegg_res[[1]] %>%
  setReadable(., OrgDb = org.Hs.eg.db, keyType = "ENTREZID") %>%
  as_tibble() %>%
  select(-c(category, subcategory, BgRatio, zScore, pvalue))

total_kegg_res <- total_kegg_res %>%
  mutate(Description = map_chr(total_kegg_res$Description, ~make_pretty(.x))) %>%
  slice_head(n = 12)


pathwaynames <- total_kegg_res$Description

pathway_mat <- matrix(NA_character_, nrow = nrow(common_dges),
                      ncol = nrow(total_kegg_res),
                      dimnames = list(common_dges$SYMBOL, pathwaynames))

for(i in seq_along(pathwaynames)){
  genelist <- rownames(pathway_mat)
  pathwaygenes <- total_kegg_res %>% filter(Description == pathwaynames[i]) %>%
    pull(geneID) %>% strsplit(., "/") %>% unlist()
  
  pathway_mat[genelist %in% pathwaygenes, i] <- pathwaynames[i]
}

genes_keep <- pathway_mat %>% rowAlls() %>% as_tibble(rownames = "gene")

pathway_mat <- as_tibble(pathway_mat, rownames = "gene") %>%
  left_join(., genes_keep, by = "gene") %>%
  drop_na(value) %>%
  column_to_rownames("gene") %>%
  select(-value) %>% 
  as.matrix()

## Steps to make heatmap for plaques

sub_common_res <- common_res %>%
  filter(SYMBOL %in% rownames(pathway_mat)) %>% 
  mutate(total_score = abs(log2FoldChange) + abs(mouse_log2fc)) %>%
  arrange(desc(total_score))

plaque_common_exp <- sub_common_res %>% select(SYMBOL, Stable_Pt_1:padj)

plaque_scaled_exp <- plaque_common_exp %>%
  select(SYMBOL, starts_with("Stable"), starts_with("Unstable")) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix() %>%
  t() %>%
  scale() %>%
  t()

z_col_fun <- colorRampPalette(colors = c('blue','white','red'))(250)

# z_col_fun <- colorRamp2(c(-2, 0, 2), c("#313695", "#FFFFBF", "#D73027"))

plaque_scaled_padj <- plaque_common_exp %>%
  select(SYMBOL, padj) %>%
  mutate(padj = -log10(padj)) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

padj_col_fun <- circlize::colorRamp2(c(0, 8, 16),   c("#4A148C", "#AAAAAA", "#D50000"))

# Column groups (Bottom color tracks)
anno_plaque <- HeatmapAnnotation(
  Group = c(rep("Stable", 4), rep("Unstable", 4)),
  col = list(Group = c("Stable" = "#5C6BC0", "Unstable" = "#EC407A")),
  show_annotation_name = F,
  show_legend = T
)

# Track 1: Heatmap for plaque dataset
ht_plaque <- Heatmap(plaque_scaled_exp,
                     col = z_col_fun,
                     bottom_annotation = anno_plaque,
                     show_row_names = F, 
                     show_column_names = F,
                     column_title = "Human Plaques",
                     show_heatmap_legend = F,
                     heatmap_legend_param = list(direction = "horizontal"))

# Track 2: Heatmap for plaque padj
ht_plaque_padj <- Heatmap(plaque_scaled_padj, 
                  name = "Human plaques\n-log10(FDR)", 
                  col = padj_col_fun, 
                  width = unit(4, "mm"), # Thin vertical bar
                  show_row_names = F, 
                  show_column_names = F,
                  row_names_gp = gpar(fontsize = 7.5),
                  heatmap_legend_param = list(direction = "horizontal"))

## Steps to make heatmap for plaques
mouse_common_exp <- sub_common_res %>% select(SYMBOL, `24hrKO1`:mouse_padj)

mouse_scaled_exp <- mouse_common_exp %>%
  select(SYMBOL, starts_with("24hr")) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix() %>%
  t() %>%
  scale() %>%
  t()

mouse_scaled_padj <- mouse_common_exp %>%
  select(SYMBOL, mouse_padj) %>%
  mutate(mouse_padj = -log10(mouse_padj)) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

mouse_padj_col_fun <- circlize::colorRamp2(
  c(0, 30, 60),c("#4A148C", "#AAAAAA", "#D50000"))

# Column groups (Bottom color tracks)
anno_mouse <- HeatmapAnnotation(
  Group = c(rep("WT", 3), rep("MKO", 3)),
  col = list(Group = c("WT" = "lightblue", "MKO" = "#FFA726")),
  show_annotation_name = F,
  show_legend = T,
  annotation_legend_param = list(direction = "horizontal"),
  which = "column"
)


# Track 1: Heatmap for plaque dataset
ht_mouse <- Heatmap(mouse_scaled_exp, 
                     name = "Z-Score", 
                     col = z_col_fun,
                     bottom_annotation = anno_mouse,
                     show_row_names = F, 
                     show_column_names = F,
                     column_title = "Mouse Macrophages",
                    heatmap_legend_param = list(direction = "horizontal"))

# Track 2: Heatmap for plaque padj
ht_mouse_padj <- Heatmap(mouse_scaled_padj, 
                          name = "Mouse Macs\n-log10(FDR)", 
                          col = mouse_padj_col_fun, 
                          width = unit(4, "mm"), # Thin vertical bar
                          show_row_names = T, 
                          show_column_names = F,
                         heatmap_legend_param = list(direction = "horizontal"))


pathway_mat <- pathway_mat[rownames(plaque_scaled_exp), , drop = FALSE]

# Discrete colors for the Pathway Grid columns
pathway_cols <- hcl.colors(n = length(pathwaynames)) %>%
  setNames(., pathwaynames)

# Heatmap 5: Functional Pathway Grid
ht_pathways <- Heatmap(pathway_mat, 
                       name = "KEGG Pathways", 
                       col = pathway_cols,
                       na_col = "#FFFFFF",
                       rect_gp = gpar(col = "#444444",
                                      lwd = 0.5),
                       show_column_names = FALSE, 
                       row_names_side = "right",
                       row_names_gp = gpar(fontsize = 9),
                       cluster_columns = FALSE,
                       # width = unit(30, "mm")
)
ht_list <- ht_plaque + ht_plaque_padj + ht_mouse + ht_mouse_padj + ht_pathways

pdf("~/bulk_slc38a9_rnaseq_aq/results/r/plaque_24hr_heatmap.pdf",
    width = 7, height = 7)

draw(ht_list, 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right",
     merge_legends = TRUE)

dev.off()
