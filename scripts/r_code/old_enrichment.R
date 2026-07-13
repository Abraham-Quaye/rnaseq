#!/usr/bin/env Rscrip --vanilla

library(magrittr)
library(ggrepel)
library(ggtext)
library(Mus.musculus)
library(clusterProfiler)
library(org.Mm.eg.db)
library(pathview)
library(tidyverse)
# the data is located in the "results/tables" folder
deg_files <- list.files("results/tables",
                        pattern = "^sig\\w+_\\w+DEGs\\.csv$",
                        full.names = TRUE) %>%
  set_names(map_chr(., ~sub("^results/tables/significant_(\\w+)_DEGs\\.csv$",
                         "\\1", .x)))

deg_data <- map(deg_files, ~read_csv(file = .x, id = "contrasts")) %>%
  list_rbind() %>%
  dplyr::rename(log2fc = log2FoldChange) %>%
  dplyr::mutate(contrasts = map_chr(contrasts,
                                ~sub("^\\w+/\\w+/\\w+cant_(\\w+)_DEGs\\.csv$",
                                 "\\1", .x)),
                regulation = case_when(log2fc >= 0 ~ "up",
                                log2fc < 0 ~ "down",
                                TRUE ~ "not_sig"))

############################# DEG bar plot ##################
deg_bar_plt <- deg_data %>%
  dplyr::summarise(num_genes = n(),
                   .by = c(regulation, contrasts)) %>%
  dplyr::mutate(contrasts = factor(contrasts, levels = c(names(deg_files)[1],
                                                         names(deg_files)[3],
                                                         names(deg_files)[2],
                                                         names(deg_files)[4]))) %>%
  ggplot(aes(contrasts, num_genes, fill = regulation)) +
  geom_col(position = position_dodge(0.9)) +
  geom_text(aes(label = num_genes), position = position_dodge(width = 0.9),
            vjust = -0.5, fontface = "bold", size = 6.5) +
  scale_fill_manual(values = c("#0000FF", "#ff0000"),
                    breaks = c("down", "up"),
                    labels = c("Downregulated", "Upregulated")) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0.16, 0.16),
                   breaks = names(deg_files),
                   labels = c("KO+AC vs WT+AC", "KO+AC vs WT",
                              "KO vs WT", "WT+AC vs WT")
                   ) +
  coord_cartesian(clip = "off", ylim = c(0, 549)) +
  labs(title = "DEGs of SLC38A9 KO Macrophages Cultured with or Without ACs <br>Compared to WT",
       y = "Number of Genes",
       x = NULL,
       fill = NULL) +
  theme_classic() +
  theme(plot.margin = margin_auto(10, 10),
        plot.title.position = "plot",
        plot.title = element_markdown(face = 'bold', size = 18,
                                      colour = 'black', hjust = 0,
                                      lineheight = 1.3, vjust = 0.5),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major.y = element_line(color = 'grey50',
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.minor = element_line(linewidth = 0.1, colour = "grey50",
                                        linetype = "dashed"),
        #adjust axis
        axis.text.x = element_markdown(size = 12, colour = 'black', face = 'bold',
                                       margin = margin(t = 10)),
        axis.text.y = element_text(size = 15, colour = 'black', face = 'bold',
                                   margin = margin(l = 10, r = 4)),
        axis.title.y = element_text(size = 18,
                                    face = 'bold',
                                    color = 'black'),
        axis.ticks.y = element_blank(),
        axis.ticks.length.x = unit(6, "pt"),
        legend.justification = c(0,1),
        legend.position = "inside",
        legend.position.inside = c(0.05, 0.8685),
        legend.box.background = element_rect(color = "grey60",
                                             linewidth = 0.1),
        legend.text = element_text(face = 'bold', size = 15),
        legend.background = element_blank(),
        legend.margin = margin_auto(4, 4),
        legend.key.size = unit(0.9, "cm"),
        legend.key.spacing.y = unit(0.3, "cm")
  )

ggsave(plot = deg_bar_plt, "results/figures/deg_bar.pdf",
       width = 8, height = 8)


########################### FUNCTIONAL ENRICHMENT ANALYSES #################

############## translate gene ids to entrez ids

# Find supported id types using keytypes(<OrgDb>)
keys <- keytypes(org.Mm.eg.db)

# Function for fetching and processing KEGG enrichment results
fetch_kegg_enrichment <- function(data){

  entrez_list <- data$entrezid

  kegg_res <- enrichKEGG(entrez_list, organism = "mmu") %>%
    as_tibble() %>%
    dplyr::mutate(overlap_genes = map(geneID ,~str_split(.x, "/") %>%
                                        unlist()),
                  gene_names = map_chr(overlap_genes,
                                       ~AnnotationDbi::select(
                                         Mus.musculus, keys = .x,
                                         columns = c("SYMBOL"),
                                         keytype = "ENTREZID") %>%
                                         distinct(ENTREZID, .keep_all = T) %>%
                                         pull(SYMBOL) %>%
                                         base::unlist() %>%
                                         paste0(., collapse = ", ")))
  return(kegg_res)
}

# Function for fetching and processing GO enrichment results
fetch_go_enrichment <- function(data, ontology = "ALL"){

  entrez_list <- data$entrezid

  enrichGO(gene = entrez_list,
           keyType = "ENTREZID",
           OrgDb = org.Mm.eg.db,
           ont = ontology,
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05,
           readable = T)
}

# 1. Separate data by condition
# 2. Get other gene annotations for significant DEGs
# 3-5. Join annotation info to DEG expression info
# 6-8. Extract log2fc values for each corresponding gene id as its name
# 9-11. Fetch significant KEGG pathways
# 12-14. Fetch "ALL" significant GO Terms

sig_genes_by_cond <- deg_data %>%
  nest(.by = contrasts, .key = "all_sig") %>%
  dplyr::mutate(all_sig = map(all_sig, ~select_if(.x, ~!any(is.na(.x)))),
                gene_annot = map(all_sig,
                                 ~AnnotationDbi::select(Mus.musculus,
                                                        keys = .x$gene_id,
                                                        columns = c("SYMBOL",
                                                                    "GENENAME",
                                                                    "ENTREZID"),
                                                        keytype = "ENSEMBL") %>%
                                   distinct(ENSEMBL, .keep_all = T) %>%
                  rename_all(tolower)),
                annot_all_sig = map2(all_sig, gene_annot,
                                     ~left_join(.x, .y,
                                                by = c("gene_id" = "ensembl")) %>%
                                     dplyr::select(-c(stat, pvalue)) %>%
                                     dplyr::select(gene_id, entrezid, symbol,
                                                   genename,
                                                   starts_with(c("KO", "WT")),
                                                   log2fc, padj, regulation)),
                annot_up_sig = map(annot_all_sig,
                                   ~dplyr::filter(.x, regulation == "up") %>%
                                     dplyr::select(-regulation)),
                annot_down_sig = map(annot_all_sig,
                                     ~dplyr::filter(.x, regulation == "down") %>%
                                       dplyr::select(-regulation)),
                all_sig_list = map(annot_all_sig,
                                     ~pull(.x, log2fc) %>%
                                     set_names(.x$entrezid)),
                up_sig_list = map(annot_up_sig,
                                    ~pull(.x, log2fc) %>%
                                    set_names(.x$entrezid)),
                down_sig_list = map(annot_down_sig,
                                    ~pull(.x, log2fc) %>%
                                      set_names(.x$entrezid)),
                all_kegg_res = map(annot_all_sig, fetch_kegg_enrichment),
                up_kegg_res = map(annot_up_sig, fetch_kegg_enrichment),
                down_kegg_res = map(annot_down_sig, fetch_kegg_enrichment),
                all_go_res = map(annot_all_sig, ~fetch_go_enrichment(.x) %>%
                                   as_tibble()),
                up_go_res = map(annot_up_sig, ~fetch_go_enrichment(.x) %>%
                                  as_tibble()),
                down_go_res = map(annot_down_sig, ~fetch_go_enrichment(.x) %>%
                                    as_tibble())
                )

find_target_DEGs <- function(ko_sig_df, wt_sig_df, ko_ref_wt_df){

  only_ko <- anti_join(ko_sig_df, wt_sig_df, by = "gene_id")

  targets <- inner_join(ko_ref_wt_df, only_ko, by = "gene_id") %>%
    dplyr::select(gene_id, entrezid = entrezid.x, symbol = symbol.x,
                  genename = genename.x, log2fc_vs_AC = log2fc.x,
                  log2fc_vs_noAC = log2fc.y)

  return(targets)
}

# downregulated targets
down_targets <- find_target_DEGs(ko_sig_df = sig_genes_by_cond$annot_down_sig[[2]],
                                 wt_sig_df = sig_genes_by_cond$annot_down_sig[[4]],
                                 ko_ref_wt_df = sig_genes_by_cond$annot_down_sig[[1]])
# interesting genes: Pip5k1c


# upregulated targets
up_targets <- find_target_DEGs(ko_sig_df = sig_genes_by_cond$annot_up_sig[[2]],
                               wt_sig_df = sig_genes_by_cond$annot_up_sig[[4]],
                               ko_ref_wt_df = sig_genes_by_cond$annot_up_sig[[1]])
# for follow-up: sestrin1 (sesn1),

# function to fetch KEGG pathway diagrams with DEG expression overlay

view_kegg_pathway <- function(g_name_lfc, pathway){

  id <- pathway$ID
  pathview(gene.data = g_name_lfc,
           species = "mmu",
           pathway.id = id,
           res = 500,
           kegg.dir = "results/figures",
           limit = list(gene = 1, cpd = 1),
           low = list(gene = "blue", cpd = "blue"),
           mid = list(gene = "grey", cpd = "grey"),
           high = list(gene = "red", cpd = "red"))
}


map2(sig_genes_by_cond$all_sig_list,
     sig_genes_by_cond$up_kegg_res,
     view_kegg_pathway)



####################### GO Analysis #######################

# View a figure of the interconnections of the GO terms
# This plot does not work for "ALL" GO terms together
all_rich_go_bp <- fetch_go_enrichment(sig_genes_by_cond$annot_all_sig[[1]], "BP")
go_connections <- goplot(all_rich_go_bp)

# Results for only upregulated genes
up_rich_go_all <- fetch_go_enrichment(up_genes, "ALL") # no result

up_go_all_table <- up_rich_go_all %>%
  as_tibble(rownames = NULL)

# Results for only downregulated genes
down_rich_go_all <- fetch_go_enrichment(down_genes, "ALL")

down_go_all_table <- down_rich_go_all %>%
  as_tibble(rownames = NULL)
