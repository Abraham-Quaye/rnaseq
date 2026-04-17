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
result_path <- "~/bm_fn_rnaseq/results/r/"

contr <- "FN_vs_BM"
# load functions ==================
source("scripts/r_code/enrichment_analysis_functions.R")

# get all file paths

# read all data into one big dataframe
deg_data <- read_csv(file = paste0(
  result_path,"tables/significant_" , contr, "_DEGs.csv"
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
  contr_name = contr,
  all_g_list = map(data, ~pull(.x, log2FoldChange) %>% # for plot_kegg_pathway
                     set_names(., .x$ENTREZID)),
  up_genes = map(data, ~get_subset_genes(genes_tbl = .x, reg = "up")),
  down_genes = map(data, ~get_subset_genes(genes_tbl = .x, reg = "down")),
  # KEGG results
  total_kegg = map(data, ~enrichKEGG(gene = .x$ENTREZID,
                                     organism = "mmu")),
  up_kegg = map(up_genes, ~enrichKEGG(gene = .x$ENTREZID,
                                      organism = "mmu")),
  down_kegg = map(down_genes, ~enrichKEGG(gene = .x$ENTREZID,
                                          organism = "mmu")),
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

##### Make Plots of Target Pathways ----
enrichkegg <- setReadable(enrich_result$total_kegg[[1]],
                          OrgDb = org.Mm.eg.db,
                          keyType = "ENTREZID")

sub_enrichkegg <- enrichkegg
target_pathways <- c("Efferocytosis", "Focal adhesion", "Lysosome biogenesis",
                     "Lipid and atherosclerosis",
                     "Phagosome", "Integrin signaling",
                     "Fluid shear stress and atherosclerosis",
                     "Rheumatoid arthritis")

sub_enrichkegg@result <- sub_enrichkegg@result %>%
  mutate(Description = str_remove_all(Description,
                                      " - Mus musculus \\(house mouse\\)"),
         p.adjust = as.numeric(p.adjust)) %>% 
  filter(Description %in% target_pathways)

p_lims <- pull(sub_enrichkegg@result, p.adjust)

interest_kegg_dot <- dotplot(object = sub_enrichkegg,
                             showCategory = 10,
        title = "Significantly Enriched KEGG Pathways",
        font.size = 10.5) +
  scale_fill_gradientn(name = "P adjusted",
                        colors = c("blue", "white", "red"),
                        values = scales::rescale(
                          c(max(p_lims), 0, min(p_lims))
                        )) +
  labs(size =  "Gene\nCount") +
  theme(plot.title = element_text(face = "bold",
                                  size = 15, hjust = 0.5), 
        legend.title = element_text(face = "bold")) 

ggsave(plot = interest_kegg_dot,
       paste0(result_path,
              "figures/kegg_FN_vs_BM_interestPathways_dotplot.pdf"),
       width = 8, height = 5.5)

gsea_list <- deg_data$log2FoldChange %>%
  setNames(., deg_data$ENTREZID) %>%
  sort(decreasing = T)

gsea_list <- gsea_list[!is.na(names(gsea_list))]

gsea_list <- gsea_list[!duplicated(names(gsea_list))]

cnet1 <- cnetplot(sub_enrichkegg, foldChange = gsea_list,
                  showCategory = 8, color_category = "black",
                  node_label = "item", curvature = 0.2) +
  geom_cnet_label(node_label = "category",
                  size = 4.5, fontface = "bold",
                  box.padding = 0, max.overlaps = 40,
                  color  = "#000000") +
  scale_color_gradientn(colors = c("blue", "grey", "red"),
                        values = scales::rescale(
                          c(min(gsea_list), 0, max(gsea_list))
                        ),
                        name = "log2FC",
                        breaks = c(-1, 0, 2, 4),
                        labels = c(-1, 0, 2, 4)) +
  theme(panel.background = element_rect(fill = NA),
        plot.margin = margin_auto(10)) +
  guides(size = guide_legend(title = "Gene\nCount"))

ggsave(plot = cnet1,
       paste0(result_path,
              "figures/cnet_enrich_pathways.pdf"),
       width = 12, height = 10)

# ChordPlot
target_plabs <- c("Efferocytosis",
                  "Fluid shear stress\nand atherosclerosis",
                  "Focal adhesion",
                  "Integrin signaling",
                  "Lipid and\natherosclerosis",
                  "Lysosome biogenesis",
                  "Phagosome", 
                  "Rheumatoid arthritis")

kegg_sub_df <- sub_enrichkegg@result %>% as_tibble() %>%
  select(term = Description, ID, genes = geneID, adj_pval = p.adjust) %>%
  mutate(genes = str_replace_all(genes, "/", ", ") %>%
           toupper() %>% trimws(),
         category = "KEGG") %>%
  select(category, ID, term, genes, adj_pval) %>%
  as.data.frame()

genes_df <- deg_data %>%
  select(ID = SYMBOL, logFC = log2FoldChange) %>%
  mutate(ID = trimws(toupper(ID)),
         logFC = as.numeric(logFC)) %>%
  distinct(ID, .keep_all = TRUE) %>% 
  as.data.frame()

circ_data <- circle_dat(kegg_sub_df, genes_df) %>%
  group_by(term) %>%
  arrange(desc(abs(logFC)), .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  as.data.frame()

# chord_mat <- chord_dat(circ_data)
chord_mat <- chord_dat(data = circ_data,
                       genes = unique(circ_data$genes),
                       process = unique(circ_data$term))

# 3. THE FIX: Manually inject the logFC if it's missing
if (!"logFC" %in% colnames(chord_mat)) {
  # Create a lookup vector from your genes_df
  lfc_lookup <- genes_df$logFC
  names(lfc_lookup) <- genes_df$ID
  
  # Match the logFC to the row names of your chord matrix
  chord_mat <- cbind(chord_mat, logFC = lfc_lookup[rownames(chord_mat)])
}

chordplot <- GOChord(chord_mat, space = 0.01,
        gene.order = 'logFC',
        gene.size = 3,
        gene.space = 0.2,
        lfc.col = c("red", "white", "blue"),
        lfc.max = 10, lfc.min = -10,
        border.size = 0,
        ribbon.col = viridis::viridis(n = 8)) +
  annotate(geom = "text",
           x = c(0.38, 0.98, 1.17, 1.35, 1.32, 1.28, 0.9, 0.45),
           y = c(1.12, 0.95, 0.65, 0.25, -0.2, -0.6, -0.9, -1.15),
           label = target_plabs,
           fontface = "bold", size = 6) +
  coord_equal(clip = "off") +
  scale_fill_gradientn(name = "log2FC\n(FN vs BM)",
                       colors = c("red", "grey", "blue"),
                       values = scales::rescale(
                         c(max(circ_data$logFC), 0,
                           min(circ_data$logFC))
                       )) +
  theme(plot.margin = margin(r = 100, l = 0,
                             t = 0, b = 0),
        legend.box.margin = margin(l = 15),
        legend.title.position = "top",
        legend.text = element_text(hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.12, 0.1),
        text = element_text(face = "bold")) +
  guides(shape = "none",
         size = "none")

chordplot$layers[[2]]$aes_params$colour <- NULL
chordplot$layers[[5]]$aes_params$colour <- NULL
chordplot$layers[[6]]$aes_params$colour <- NULL
chordplot$layers[[1]]$aes_params$colour <- NULL

ggsave(plot = chordplot,
       paste0(result_path,
              "figures/kegg_FN_vs_BM_chordplot.pdf"),
       width = 12, height = 12)

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

kegg_diagram_dir <- paste0(result_path,
                           "figures/kegg_pathway_diagrams")

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
