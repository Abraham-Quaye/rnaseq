#!/usr/bin/env Rscript --vanilla

library(magrittr)
library(Homo.sapiens)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(pathview)
library(ggrepel)
library(ggtext)
library(tidyverse)

# the data is located in the "results/tables" folder
result_path <- "~/myocd_rnaseq/results/r/"

# load functions ==================
source("scripts/r_code/enrichment_analysis_functions.R")

# get all file paths

# read all data into one big dataframe
deg_data <- read_csv(file = paste0(
  result_path,"tables/significant_MYOCD_vs_GFP_DEGs.csv"
  )) %>% 
  mutate(regulation = case_when(log2FoldChange >= 0 ~ "up",
                                log2FoldChange < 0 ~ "down",
                                TRUE ~ NA_character_),,
         ENTREZID = as.character(ENTREZID)) %>%
  arrange(padj)

# FUNCTIONAL ENRICHMENT ANALYSES =======================================

# Perform functional enrichment analyses
enrich_result <- tibble(
  data = list(deg_data),
  contr_name = "MYOCD_vs_GFP",
  all_g_list = map(data, ~pull(.x, log2FoldChange) %>% # for plot_kegg_pathway
                     set_names(., .x$ENTREZID)),
  up_genes = map(data, ~get_subset_genes(genes_tbl = .x, reg = "up")),
  down_genes = map(data, ~get_subset_genes(genes_tbl = .x, reg = "down")),
  # KEGG results
  total_kegg = map(data, ~enrichKEGG(gene = .x$ENTREZID)),
  up_kegg = map(up_genes, ~enrichKEGG(gene = .x$ENTREZID)),
  down_kegg = map(down_genes, ~enrichKEGG(gene = .x$ENTREZID)),
  # KEGG dotplots
  total_kegg_dotplot = map2(total_kegg, contr_name,
                            ~plot_dotplot(res = .x, labb = .y)),
  up_kegg_dotplot = map2(up_kegg, contr_name,
                         ~plot_dotplot(res = .x, labb = .y)),
  down_kegg_dotplot = map2(down_kegg, contr_name,
                           ~plot_dotplot(res = .x, labb = .y)),
  # GO results ================
  # for goplot
  total_go_bp = map(data, ~get_go_enrich(.x, "BP")),
  total_go_bp_connects = map(total_go_bp, ~goplot(.x)),
  total_go_all = map(data, ~get_go_enrich(.x, "ALL")),
  up_go_all = map(up_genes, ~get_go_enrich(.x, "ALL")),
  down_go_all = map(down_genes, ~get_go_enrich(.x, "ALL")),
  # GO dotplots
  total_go_dotplot = map2(total_go_all, contr_name,
                          ~plot_dotplot(res = .x, labb = .y)),
  up_go_dotplot = map2(up_go_all, contr_name,
                       ~plot_dotplot(res = .x, labb = .y)),
  down_go_dotplot = map2(down_go_all, contr_name,
                         ~plot_dotplot(res = .x, labb = .y))
)

# Save enrichment results =====================
# Save all KEGG enrichment results
save_kegg_results(res = enrich_result$total_kegg,
                  contr_name = enrich_result$contr_name,
                  labb = "_total")

save_kegg_results(res = enrich_result$up_kegg,
                  contr_name = enrich_result$contr_name,
                  labb = "_up")

save_kegg_results(res = enrich_result$down_kegg,
                  contr_name = enrich_result$contr_name,
                  labb = "_down")

# Save all GO enrichment results
save_go_results(res = enrich_result$total_go_all,
                contr_name = enrich_result$contr_name,
                labb = "_total")

save_go_results(res = enrich_result$up_go_all,
                contr_name = enrich_result$contr_name,
                labb = "_up")

save_go_results(res = enrich_result$down_go_all,
                contr_name = enrich_result$contr_name,
                labb = "_down")

# Save Dotplots =========================
# KEGG dotplots
save_kegg_dotplots(dotplots = enrich_result$total_kegg_dotplot,
                   contr_name = enrich_result$contr_name,
                   labb = "_total")

save_kegg_dotplots(dotplots = enrich_result$up_kegg_dotplot,
                   contr_name = enrich_result$contr_name,
                   labb = "_up")

save_kegg_dotplots(dotplots = enrich_result$down_kegg_dotplot,
                   contr_name = enrich_result$contr_name,
                   labb = "_down")

# GO dotplots
save_go_dotplots(dotplots = enrich_result$total_go_dotplot,
                 contr_name = enrich_result$contr_name, 
                 labb = "_total")

save_go_dotplots(dotplots = enrich_result$up_go_dotplot,
                 contr_name = enrich_result$contr_name, 
                 labb = "_up")

save_go_dotplots(dotplots = enrich_result$down_go_dotplot,
                 contr_name = enrich_result$contr_name, 
                 labb = "_down")

# Plot KEGG Pathway Diagrams
kegg_to_plot <- enrich_result %>%
  pull(total_kegg) %>% pluck(1)

genes_to_plot <- enrich_result %>%
  pull(all_g_list) %>% unlist()

kegg_diagram_dir <- paste0(result_path, "figures/kegg_pathway_diagrams")

dir.create(kegg_diagram_dir, recursive = T)
orig_wd <- getwd()
setwd(kegg_diagram_dir)

map(kegg_to_plot$ID,
    ~safe_kegg_plotter(
      id = .x,
      gene_list = genes_to_plot,
      path_ = kegg_diagram_dir
      )
    )

setwd(orig_wd)
