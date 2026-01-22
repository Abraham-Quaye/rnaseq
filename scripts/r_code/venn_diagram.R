#!/usr/bin/env Rscript

library(ggVennDiagram)
library(ggvenn)
library(patchwork)

source("scripts/r_code/extract_deg_geneIDs.R")

up <- all_deg_tables[c("up_degIDs_4hrs", "up_degIDs_12hrs",
                       "up_degIDs_24hrs")]

v_up <- list(`4-hpi` = up$up_degIDs_4hrs %>% pull(gene_id),
             `12-hpi` = up$up_degIDs_12hrs %>% pull(gene_id),
             `24-hpi` = up$up_degIDs_24hrs %>% pull(gene_id))

down <- all_deg_tables[c("down_degIDs_4hrs", "down_degIDs_12hrs",
                         "down_degIDs_24hrs", "down_degIDs_72hrs")]

v_down <- list(`4-hpi` = down$down_degIDs_4hrs %>% pull(gene_id),
               `12-hpi` = down$down_degIDs_12hrs %>% pull(gene_id),
               `24-hpi` = down$down_degIDs_24hrs %>% pull(gene_id),
               `72-hpi` = down$down_degIDs_72hrs %>% pull(gene_id))


draw_venn_plott <- function(input_sets, tittl){
  
  pl <- ggvenn(input_sets, stroke_size = 0.1,
               set_name_size = 10, show_percentage = F,
               text_size = 10
               ) +
    coord_cartesian(clip = "off") +
    labs(title = paste0(tittl, "regulated DEGs")) +
    theme(plot.title = element_text(size = 22, face = "bold",
                                      colour = "black", hjust = 0.5))

  return(pl)
}

plts <- map2(list(v_up, v_down), list("Up", "Down"),
             ~draw_venn_plott(input_sets = .x,tittl = .y)) %>%
  set_names(c("up", "down"))

fplot <- (plts[["up"]] | plot_spacer() | plts[["down"]]) +
  plot_layout(widths = c(1, 0.05, 1), tag_level = "new") +
  plot_annotation(title = "Overlaps of Differentially Expressed Genes Over Time",
                  tag_levels = 1, tag_prefix = "C",
                  theme = theme(plot.title = element_text(size = 25,
                                                          hjust = 0.5,
                                                          face = "bold",
                                                          margin = margin(t = 10,
                                                                          b = 10)))
                  ) &
  theme(plot.tag = element_text(size = 22, face = "bold",
                                margin = margin(b = -20)))
