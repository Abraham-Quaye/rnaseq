#!/usr/bin/env Rscrip --vanilla

library(tidyverse)
library(readxl)
library(magrittr)
library(ggrepel)

# Read in excel file and process it to include only the significant genes
exc_files <- list.files("raw_files",
                       pattern = "^(r77q|wt)_vs_(mock|wt)\\.xlsx",
                       full.names = T)
names(exc_files) <- sub("raw_files/((r77q|wt)_vs_(mock|wt))\\.xlsx",
                        "\\1", exc_files)

sheet_name <- paste0(c(4, 8, 12, 24, 72), "hr")
names(sheet_name) <- sheet_name

# the data is located in the "raw_files" folder
get_raw_expr_data <- function(file){
  map_dfr(sheet_name, ~read_xlsx(file, sheet = .x), .id = "timepoint") %>%
    dplyr::select(timepoint, gene_id = ensembl_gene_id,
                gene_name, qval = FDR, log2fc = logFC) %>%
    dplyr::mutate(significant = ifelse(qval <= 0.05, TRUE, FALSE),
           regulation = case_when(significant & log2fc >= 0 ~ "up",
                                  significant & log2fc < 0 ~ "down",
                                  TRUE ~ "not_sig"))
}

raw_data <- map_dfr(exc_files, ~get_raw_expr_data(.x), .id = "dataset")

############################# Volcano plot ##################
# We are doing this for only the 72hr timepoints

plot_volcano <- function(expr_data, most_sig){

  expr_data %>%
    ggplot(aes(log2fc, -log10(qval), colour = regulation)) +
    geom_point(alpha = 0.8, size = 2.5) +
    geom_label_repel(data = most_sig,
                    aes(log2fc, log_qval, label = gene_name),
                    min.segment.length = 0, fontface = "bold",
                    max.overlaps = 25, box.padding = 1.5,
                    point.padding = 0.5, size = 3.5,
                    show.legend = F) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_y_continuous(expand = c(0, 0.1)) +
    coord_cartesian(clip = "off") +
    scale_color_manual(values = c("blue", "grey", "red"),
                       breaks = c("down", "not_sig", "up"),
                       labels = c("Downregulated", "Not Significant", "Upregulated")) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5,
                                    margin = margin_auto(0, 0)),
          legend.margin = margin_auto(0, 0),
          plot.title.position = "plot",
          plot.margin = margin(10, 10, 10, 10),
          axis.title = element_text(size = 18, face = "bold"),
          axis.text = element_text(size = 16, colour = "#000000"),
          legend.title = element_blank(),
          legend.text = element_text(size = 12, face = "bold",
                                     margin = margin(r = 30, l = 0, b = 0, t = 0)),
          legend.position = "top",
          legend.direction = "horizontal",
          legend.key.size = unit(0.7, "cm"),
          legend.key = element_rect(fill = NA)) +
    guides(color = guide_legend(override.aes = list(size = 6)))
}

plt_data <- nest(raw_data, .by = dataset, .key = "expr_data") %>%
  dplyr::mutate(expr_data = map(expr_data, ~filter(.x, timepoint == "72hr")),
                most_sig = map(expr_data,
                               \(.x){dplyr::mutate(.x, log_qval = -log10(qval)) %>%
                                   dplyr::arrange(desc(log_qval)) %>%
                                   dplyr::slice_head(n = 15)}),
                volcano = map2(expr_data, most_sig, ~plot_volcano(.x, .y)),
                volcano_title = case_match(dataset,
                                           "r77q_vs_mock" ~
                                             "R77Q VS Mock at 72 Hours Post Infection",
                                           "r77q_vs_wt" ~
                                             "R77Q VS WT at 72 Hours Post Infection",
                                           "wt_vs_mock" ~
                                             "WT VS Mock at 72 Hours Post Infection",
                                           .default = NA_character_),
                volcano = map2(volcano, volcano_title,
                               \(.x, .y){.x +
                                   labs(title = .y,
                                   y = expression("-Log"[10]*"(P-adjusted Values)"),
                                   x = expression("Log"[2]*"(Fold Change)"))}
                ))

########################### FUNCTIONAL ENRICHMENT ANALYSES #################
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
library(Homo.sapiens)

sig_data <- raw_data %>%
  dplyr::filter(qval <= 0.05) %>%
  dplyr::mutate(regulation = ifelse(log2fc >= 0, "up", "down")) %>%
  # we do this next step because only 72hr timepoint has enough genes for analyses
  dplyr::filter(timepoint == "72hr") %>%
  dplyr::select(-c(timepoint, significant))

############## translate gene ids to entrez ids
ensembl_unq <- sig_data %>%
  distinct(gene_id, gene_name) %>%
  dplyr::pull(gene_id)

# find supported id types using keytypes(<OrgDb>)
keys <- keytypes(org.Hs.eg.db)

annot_sig_genes <- AnnotationDbi::select(Homo.sapiens,
                                         keys = ensembl_unq,
                                         columns = c("ENTREZID", "SYMBOL","GENENAME"),
                                         keytype = "ENSEMBL") %>%
  distinct(ENSEMBL, ENTREZID, .keep_all = T) %>%
  as_tibble() %>%
  rename_all(tolower)

# add Entrez ID data to significant gene table (these are all uniqe,
# unlike ensembl IDs which can map to multiple gene symbols)
sig_data <- sig_data %>% left_join(., annot_sig_genes,
                                 by = c("gene_id" = "ensembl"),
                                 relationship = "many-to-many") %>%
  dplyr::mutate(symbol = ifelse(str_detect(symbol, "^LOC") | is.na(symbol),
                                gene_name, symbol)) %>%
  dplyr::select(dataset, ensembl = gene_id, entrez = entrezid, symbol,
                genename, log2fc, qval, regulation, gene_name) %>%
  arrange(dataset, ensembl)

all_g_list <- sig_data %>% pull(log2fc) %>%
  set_names(sig_data$entrez)

# filter data for only upregulated data
up_genes <- sig_data %>% dplyr::filter(is.na(entrez))
  dplyr::filter(regulation == "up") %>%
  drop_na(entrez) %>%
  distinct(entrez, .keep_all = T) %>%
  dplyr::select(-regulation)

up_g_list <- up_genes %>% pull(log2fc) %>%
  set_names(up_genes$entrez)

# filter data for only downregulated data
down_genes <- sig_data %>%
  dplyr::filter(regulation == "down") %>%
  drop_na(entrez) %>%
  distinct(entrez, .keep_all = T) %>%
  dplyr::select(-regulation)

down_g_list <- down_genes %>% pull(log2fc) %>%
  set_names(down_genes$entrez)

up_kegg <- enrichKEGG(up_genes$entrez)

up_kegg_tab <- up_kegg %>%
  as_tibble(rownames = NULL)

down_kegg <- enrichKEGG(down_genes$entrez)

down_kegg_tab <- down_kegg %>%
  as_tibble(rownames = NULL)

all_kegg <- enrichKEGG(sig_data$entrez)
all_kegg_tab <- all_kegg %>% as_tibble(rownames = NULL)

# to run this function, make sure that you change the "kegg.dir" argument
# to a folder you want or just remove that entire argument completely
# to have the figures show up in your working directory (folder)

plot_kegg_pathway <- function(id){
  pathview(gene.data = all_g_list,
           species = "hsa",
           pathway.id = id,
           res = 500,
           kegg.dir = here("results/r/figures"),
           limit = list(gene = 1, cpd = 1),
           low = list(gene = "steelblue1", cpd = "steelblue1"),
           mid = list(gene = "grey", cpd = "grey"),
           high = list(gene = "red", cpd = "red"))
}


####################### GO Analysis #######################
get_go_enrich <- function(gene_list, ontology){
  enrichGO(gene = gene_list$ensembl,
           keyType = "ENSEMBL",
           OrgDb = org.Hs.eg.db,
           ont = ontology,
           universe = raw_data$gene_id,
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05,
           readable = T)
}

all_rich_go_all <- get_go_enrich(sig_data, "ALL")
all_go_all_table <- all_rich_go_all %>%
  as_tibble(rownames = NULL)

# View a figure of the interconnections of the GO terms
# This plot does not work for all GO terms together
all_rich_go_bp <- get_go_enrich(sig_data, "BP")
go_connections <- goplot(all_rich_go_bp)

# Results for only upregulated genes
up_rich_go_all <- get_go_enrich(up_genes, "ALL")

up_go_all_table <- up_rich_go_all %>%
  as_tibble(rownames = NULL)

# Results for only downregulated genes
down_rich_go_all <- get_go_enrich(down_genes, "ALL")

down_go_all_table <- down_rich_go_all %>%
  as_tibble(rownames = NULL)
