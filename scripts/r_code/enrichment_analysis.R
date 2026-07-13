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

## Read in Significant genes data
read_sig_data <- function(tabdir, file_pattern = "^significant_\\w+_DEGs\\.csv"){
  
  deg_files <- list.files(paste0(result_path, tabdir),
                          pattern = file_pattern,
                          full.names = TRUE)
  
  sig_data <- map(deg_files, ~read_csv(file = .x, id = "contr_name") %>%
                    select(-matches("(12hr|24hr|NoAC)"), -DEFINITION) %>%
                    mutate(contr_name = str_replace(contr_name,
                                                    "[/\\w]+(significant|total)_(\\w+)_DEGs\\.csv",
                                                    "\\2"),
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
  return(sig_data)
}

sigdata <- map(c("tables", "tables2"),
               ~read_sig_data(tabdir = .x)) %>%
  setNames(., c("sigdata_lrt", "sigdata_norm"))

strict_degs <- left_join(nest(sigdata$sigdata_lrt, lrt_data = -contr_name),
          nest(select(sigdata$sigdata_norm, -stat), norm_data = -contr_name),
          by = join_by("contr_name")) %>%
  mutate(intersection = map2(
    .x = lrt_data, .y = norm_data,
    \(lrt = .x, norm = .y){
      commgenes <- base::intersect(lrt$gene_id, norm$gene_id)
      norm %>% filter(gene_id %in% commgenes) %>%
        mutate(ENTREZID = as.character(ENTREZID)) %>% 
        arrange(desc(abs(log2FoldChange)), padj) %>%
        drop_na(ENTREZID) %>%
        filter(abs(log2FoldChange) >= 1)
      })) %>%
  select(-c(lrt_data, norm_data))


bulk_diffdge <- sigdata$sigdata_norm %>% filter(contr_name == "12hrKO_vs_12hrWT")

sn_diffdge <- read_csv("~/mouse_sn_rnaseq/results/r/tables/total_degs_macrophages.csv") %>%
  rename(sn_log2fc = log2fc, sn_significant = is_significant, sn_regulation = regulation)

inner_join(sn_diffdge, bulk_diffdge, by = join_by(gene == SYMBOL)) %>%
  filter(sn_significant) %>%
  write_csv(., paste0(result_path, "tables/enrichment/bulk.and.snrnaseqMacs_KOvsWT_dge.csv"))

sn_diffdge %>% filter(sn_significant) %>%
  write_csv(., paste0(result_path, "tables/enrichment/snrnaseqMacs_KOvsWT_sig_dge.csv"))
# get_uniq_genesets <- function(deg_data, contr1, contr2, kovswt_contr){
#   anti_join(
#     pull(filter(deg_data, contr_name == kovswt_contr),
#          intersection)[[1]],
#     pull(filter(deg_data, contr_name == contr1),
#          intersection)[[1]],
#     by = "gene_id"
#     ) %>%
#     anti_join(., pull(filter(deg_data, contr_name == contr2),
#                     intersection)[[1]], by = "gene_id") %>%
#     anti_join(., pull(filter(deg_data, contr_name == "NoACKO_vs_NoACWT"),
#                       intersection)[[1]], by = "gene_id")
# }
# 
# strict_final_degs <- tibble(res_name = c("uniq_AC12", "uniq_AC24"),
#                             contrs = list(
#                               c(contr1 = "12hrKO_vs_NoACWT",
#                                 contr2 = "12hrWT_vs_NoACWT",
#                                 kovswt_contr = "12hrKO_vs_12hrWT"),
#                               c(contr1 = "24hrKO_vs_NoACWT",
#                                 contr2 = "24hrWT_vs_NoACWT",
#                                 kovswt_contr = "24hrKO_vs_24hrWT")),
#                             uniq_dataset = map(contrs, ~get_uniq_genesets(
#                               deg_data = nest(sigdata$sigdata_norm,
#                                               intersection = -contr_name),
#                               contr1 = .x[["contr1"]],
#                               contr2 = .x[["contr2"]],
#                               kovswt_contr = .x[["kovswt_contr"]]))
#                             )

unq12WT <- anti_join(
  pull(filter(strict_degs, contr_name == "12hrKO_vs_NoACWT"),
       intersection)[[1]],
  pull(filter(strict_degs, contr_name == "12hrWT_vs_NoACWT"),
       intersection)[[1]],
  by = "gene_id"
)


## Read in Total background genes data ----
get_bg_genes <- function(tabdir_, file_patt = "^total_\\w+_DEGs\\.csv"){
  
  read_sig_data(tabdir = tabdir_, file_pattern = file_patt) %>%
  select(contr_name, gene_id, ENTREZID, SYMBOL, log2FoldChange) %>%
  nest(data = -contr_name) %>%
  mutate(bg_entrez = map(data, ~drop_na(.x, ENTREZID) %>%
                           pull(ENTREZID) %>%
                           as.character()),
         bg_logfc = map(data, \(.x){
           d <- drop_na(.x, ENTREZID)
           d %>% pull(log2FoldChange) %>% as.character() %>%
             setNames(., d$ENTREZID)
           }),
         bg_ensembl = map(data, ~drop_na(.x, gene_id) %>%
                            pull(gene_id) %>%
                            as.character()))
}

bggenes <- map(c("tables", "tables2"), ~get_bg_genes(tabdir_ = .x) %>%
                  rename(bg_data = data)) %>%
  setNames(., c("deseq_lrt", "deseq_norm"))

# FUNCTIONAL ENRICHMENT ANALYSES =======================================

# Perform functional enrichment analyses
perform_fea <- function(sigdata_, bggenes_){
  
  if("contr_name" %in% colnames(sigdata_)){
    enrich_result <- sigdata_ %>%
      nest(sig_data = -contr_name) %>%
      left_join(., bggenes_, by = join_by(contr_name))
  }else{
      enrich_result <- sigdata_
    }
  
  enrich_result <- enrich_result %>%
    mutate(sig_g_list = map(sig_data, \(.x){
      d <- drop_na(.x, ENTREZID)
      d %>% pull(log2FoldChange) %>% as.character() %>%
        setNames(., d$ENTREZID)}),
      up_genes = map(sig_data, ~get_subset_genes(.x, reg = "up")),
      down_genes = map(sig_data, ~get_subset_genes(.x, reg = "down")),
      # KEGG results
      total_kegg = map2(sig_data, bg_entrez,
                        ~enrichKEGG(
                          gene = .x$ENTREZID,
                          qvalueCutoff = 0.05,
                          universe = .y,
                          organism = "mmu")),
      up_kegg = map2(up_genes, bg_entrez,
                     ~enrichKEGG(
                       gene = .x$ENTREZID,
                       qvalueCutoff = 0.05,
                       universe = .y,
                       organism = "mmu")),
      down_kegg = map2(down_genes, bg_entrez,
                       ~enrichKEGG(
                         gene = .x$ENTREZID,
                         qvalueCutoff = 0.05,
                         universe = .y,
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
      # GO results
      # for goplot
      # total_go_bp = map(sig_data, ~get_go_enrich(.x, ontology = "BP")),
      # total_go_bp_connects = map(total_go_bp, ~goplot(.x)),
      total_go_all =  map(sig_data, ~get_go_enrich(.x, ontology = "BP")),
      up_go_all = map(up_genes, ~get_go_enrich(.x, ontology = "ALL")),
      down_go_all = map(down_genes, ~get_go_enrich(.x, ontology = "ALL")),
      # GO dotplots
      total_go_dotplot = map2(total_go_all, contr_name,
                              ~plot_dotplot(res = .x, labb = .y)),
      up_go_dotplot = map2(up_go_all, contr_name,
                           ~plot_dotplot(res = .x, labb = .y)),
      down_go_dotplot = map2(down_go_all, contr_name,
                             ~plot_dotplot(res = .x, labb = .y))
    )
  
  return(enrich_result)
}

enrich_result <- map2(sigdata, bggenes, ~perform_fea(sigdata_ = .x, bggenes_ = .y)) %>%
  list_rbind(names_to = "analysis_type")

strict_enrich <- perform_fea(sigdata_ = unnest(strict_degs, intersection),
                             bggenes_ = bggenes$deseq_lrt)

final_strict_enrich <- tibble(uniq12AC = list(unq12WT)) %>%
  mutate(sig_g_list = map(uniq12AC, ~pull(.x, log2FoldChange,
                                          name = ENTREZID)),
    up_genes = map(uniq12AC, ~get_subset_genes(.x, reg = "up")),
    down_genes = map(uniq12AC, ~get_subset_genes(.x, reg = "down")),
    # KEGG results
    total_kegg = map(uniq12AC,
                      ~enrichKEGG(
                        gene = .x$ENTREZID,
                        qvalueCutoff = 0.05,
                        universe = bggenes$deseq_norm$bg_entrez[[3]],
                        organism = "mmu")),
    up_kegg = map(up_genes,
                   ~enrichKEGG(
                     gene = .x$ENTREZID,
                     qvalueCutoff = 0.05,
                     universe = bggenes$deseq_norm$bg_entrez[[3]],
                     organism = "mmu")),
    down_kegg = map(down_genes,
                     ~enrichKEGG(
                       gene = .x$ENTREZID,
                       qvalueCutoff = 0.05,
                       universe = bggenes$deseq_norm$bg_entrez[[3]],
                       organism = "mmu"))
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

# Save Dot plots =========================
# KEGG dot plots
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
