#!/usr/bin/env Rscript --vanilla

library(magrittr)
library(Mus.musculus)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(pathview)
library(fgsea)
library(ggrepel)
library(ggtext)
library(ggarchery)
library(ggtangle)
library(GOplot)
library(circlize)
library(tidyverse)

# the data is located in the "results/tables" folder
result_path <- "~/bulk_slc38a9_rnaseq_aq/results/r/"

# load functions ==================
source("scripts/r_code/enrichment_analysis_functions.R")

# get all file paths
# write function to extract data needed for downstream analysis
deg_files <- list.files(paste0(result_path, "tables"),
                        pattern = "^significant_\\w+_DEGs\\.csv",
                        full.names = TRUE)


sig_data <- map(deg_files, ~read_csv(file = .x, id = "contr_name") %>%
                  select(-matches("(12hr|24hr|NoAC)"), -DEFINITION) %>%
                  mutate(contr_name = str_replace(contr_name,
                                                  "[/\\w]+significant_(\\w+)_DEGs\\.csv",
                                                  "\\1"),
                         timepoint = ifelse(
                           str_detect(contr_name, "^\\d{2}hr"),
                           parse_number(str_extract(contr_name, "^\\d{2}hr")) %>%
                             as.character(.),
                           "NoAC"),
                         genotype = str_remove(contr_name, "_vs_\\w+") %>%
                           str_extract(., "(KO|WT)")
                  )) %>%
  list_rbind() %>%
  mutate(ref_sample = ifelse(str_detect(contr_name, "NoAC"),
                             "Compared to WT NoAC", "Same Timepoint KO vs WT"),
         regulation = case_when(log2FoldChange >= 0 ~ "up",
                                log2FoldChange < 0 ~ "down",
                                TRUE ~ NA_character_))

bg_genes <- read_csv(paste0(result_path, "tables/total_12hrKO_vs_NoACWT_DEGs.csv")) %>%
  pull(ENTREZID) %>% as.character() %>% na.omit()


# FUNCTIONAL ENRICHMENT ANALYSES =======================================

# Perform functional enrichment analyses
enrich_result <- sig_data %>%
  nest(data = -contr_name) %>%
  mutate(all_g_list = map(data, ~pull(.x, log2FoldChange) %>% # for plot_kegg_pathway
                            set_names(., .x$ENTREZID)),
         up_genes = map(data, ~get_subset_genes(genes_tbl = .x, reg = "up")),
         down_genes = map(data, ~get_subset_genes(genes_tbl = .x, reg = "down")),
         # KEGG results
         total_kegg = map(data, ~enrichKEGG(gene = .x$ENTREZID,
                                            qvalueCutoff = 0.05,
                                            universe = bg_genes,
                                            organism = "mmu")),
         up_kegg = map(up_genes, ~enrichKEGG(gene = .x$ENTREZID,
                                             qvalueCutoff = 0.05,
                                             universe = bg_genes,
                                             organism = "mmu")),
         down_kegg = map(down_genes, ~enrichKEGG(gene = .x$ENTREZID,
                                                 qvalueCutoff = 0.05,
                                                 universe = bg_genes,
                                                 organism = "mmu")),
         # KEGG dotplots
         across(.cols = ends_with("_kegg"),
                .fns = ~map(.x, ~prettify_kegg_names(.x))),
         total_kegg_dotplot = map2(total_kegg, contr_name,
                                   ~plot_dotplot(res = .x, labb = .y)),
         up_kegg_dotplot = map2(up_kegg, contr_name,
                                ~plot_dotplot(res = .x, labb = .y)),
         down_kegg_dotplot = map2(down_kegg, contr_name,
                                  ~plot_dotplot(res = .x, labb = .y)),
         # GO results ================
         # for goplot
         total_go_bp = map(data, ~get_go_enrich(.x, "BP")),
         # total_go_bp_connects = map(total_go_bp, ~goplot(.x)),
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
enrichres_dir <- paste0(result_path, "tables/enrichment/")
enrichfig_dir <- paste0(result_path, "figures/enrichment/")

if(!dir.exists(enrichres_dir)){
  dir.create(enrichfig_dir)
}

if(!dir.exists(enrichfig_dir)){
  dir.create(enrichfig_dir)
}

save_kegg_results(res = enrich_result$total_kegg,
                  contr_name = enrich_result$contr_name,
                  labb = "_total_")

save_kegg_results(res = enrich_result$up_kegg,
                  contr_name = enrich_result$contr_name,
                  labb = "_up_")

save_kegg_results(res = enrich_result$down_kegg,
                  contr_name = enrich_result$contr_name,
                  labb = "_down_")

# Save all GO enrichment results
save_go_results(res = enrich_result$total_go_all,
                contr_name = enrich_result$contr_name,
                labb = "_total_")

save_go_results(res = enrich_result$up_go_all,
                contr_name = enrich_result$contr_name,
                labb = "_up_")

save_go_results(res = enrich_result$down_go_all,
                contr_name = enrich_result$contr_name,
                labb = "_down_")

# Save Dotplots =========================
# KEGG dotplots
save_kegg_dotplots(dotplots = enrich_result$total_kegg_dotplot,
                   contr_name = enrich_result$contr_name,
                   labb = "_total_")

save_kegg_dotplots(dotplots = enrich_result$up_kegg_dotplot,
                   contr_name = enrich_result$contr_name,
                   labb = "_up_")

save_kegg_dotplots(dotplots = enrich_result$down_kegg_dotplot,
                   contr_name = enrich_result$contr_name,
                   labb = "_down_")

# GO dotplots
save_go_dotplots(dotplots = enrich_result$total_go_dotplot,
                 contr_name = enrich_result$contr_name, 
                 labb = "_total_")

save_go_dotplots(dotplots = enrich_result$up_go_dotplot,
                 contr_name = enrich_result$contr_name, 
                 labb = "_up_")

save_go_dotplots(dotplots = enrich_result$down_go_dotplot,
                 contr_name = enrich_result$contr_name, 
                 labb = "_down_")

# # Plot KEGG Pathway Diagrams
# kegg_to_plot <- enrich_result %>%
#   pull(total_kegg) %>% pluck(1)
# 
# genes_to_plot <- enrich_result %>%
#   pull(all_g_list) %>% unlist()
# 
# kegg_diagram_dir <- paste0(result_path,
#                            "figures/kegg_pathway_diagrams")
# 
# dir.create(kegg_diagram_dir, recursive = T)
# orig_wd <- getwd()
# setwd(kegg_diagram_dir)
# 
# map(kegg_to_plot$ID,
#     ~safe_kegg_plotter(
#       id = .x,
#       gene_list = genes_to_plot,
#       path_ = kegg_diagram_dir
#       )
#     )
# 
# setwd(orig_wd)
