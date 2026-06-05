#!/usr/bin/env Rscript --vanilla

library(magrittr)
library(ggtext)
library(tidyverse)

# the data is located in the "results/tables" folder
result_path <- "~/bulk_slc38a9_rnaseq_aq/results/r/" 

# write function to extract data needed for downstream analysis
read_sig_data <- function(tabdir, file_pattern = "^significant_\\w+_DEGs\\.csv"){
  
  deg_files <- list.files(paste0(result_path, tabdir),
                          pattern = file_pattern,
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
  return(sig_data)
}

sig_data1 <- read_sig_data("tables")
sig_data2 <- read_sig_data("tables2")

############################# DEG bar plot ##################
plot_deg_bar <- function(data){
  deg_bar_plt <- data %>%
    summarise(num_genes = n(),
              .by = c(contr_name, ref_sample, regulation,
                      genotype, timepoint)) %>%
    mutate(timepoint = factor(timepoint, levels = sort(unique(timepoint)))) %>%
    ggplot(aes(contr_name, num_genes, fill = regulation)) +
    geom_col(position = position_dodge(0.9), width = 0.85) +
    geom_text(aes(label = num_genes), position = position_dodge(width = 0.9),
              vjust = -0.5, fontface = "bold", size = 4) +
    scale_fill_manual(values = c(down = "#0000FF", up = "red3"),
                      labels = c("Downregulated", "Upregulated")) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0.1, 0.1),
                     breaks = unique(data$contr_name),
                     labels = str_replace_all(unique(data$contr_name),
                                              "_", "\n")) +
    coord_cartesian(clip = "off") +
    labs(title = paste0("Differentially Expressed Genes of *Slc38a9* ",
                        "KO vs WT<br>at Different Timepoints Post Efferocytosis"),
         y = "Number of Genes",
         x = NULL,
         fill = NULL) +
    theme_classic() +
    theme(plot.margin = margin_auto(4, 4),
          plot.title.position = "plot",
          plot.title = element_markdown(face = 'bold', size = 14,
                                        colour = 'black', hjust = 0.5,
                                        lineheight = 1.3, vjust = 1, 
                                        margin = margin(2, 0, 20, 0)),
          panel.background = element_rect(fill = 'white'),
          panel.grid.major.y = element_line(color = 'grey50',
                                            linewidth = 0.035,
                                            linetype = 1),
          panel.grid.minor = element_line(linewidth = 0.1, colour = "grey50",
                                          linetype = "dashed"),
          panel.spacing.x = unit(5, "pt"),
          #adjust axis
          axis.text.x = element_text(size = 12, colour = 'black', face = 'bold',
                                     margin = margin(t = 3)),
          axis.text.y = element_text(size = 12, colour = 'black', face = 'bold',
                                     margin = margin(l = 10, r = 3)),
          axis.title.y = element_text(size = 14,
                                      face = 'bold',
                                      color = 'black'),
          axis.ticks = element_line(linewidth = 0.3, color = "black"),
          axis.ticks.length = unit(4, "pt"),
          axis.line = element_line(linewidth = 0.3),
          legend.justification = c(0,1),
          legend.position = "inside",
          legend.position.inside = c(0.01, 1.05),
          legend.box.background = element_rect(color = "black", linewidth = 0.3),
          legend.text = element_text(face = 'bold', size = 9),
          legend.background = element_blank(),
          legend.margin = margin(t = 3, r = 3, b = 3, l = 3),
          legend.key.size = unit(0.5, "cm"),
          legend.key.spacing.y = unit(0.2, "cm")
    )
  return(deg_bar_plt)
}

barplots <- list(barplot1 = plot_deg_bar(sig_data1),
                 barplot2 = plot_deg_bar(sig_data2))

map2(.x = barplots, .y = names(barplots),
     ~ggsave(plot = .x, filename = paste0(
       result_path, "figures/DEG_levels_", .y, ".png"),
       width = 7.5, height = 6.8, dpi = 350))

plot_data <- map(list(sig_data1, sig_data2),
                 \(.x){.x %>%
                     summarise(num_genes = n(),
                               .by = c(contr_name, ref_sample, regulation,
                                       genotype, timepoint)) %>%
                     mutate(total = sum(num_genes), .by = contr_name) %>%
                     mutate(
                       percent = round((num_genes/total) * 100, 1),
                       contr_name = factor(
                         contr_name, levels = unique(contr_name),
                         labels = str_replace_all(unique(contr_name), "_", " ")))
                 })

plot_deg_pie1 <- function(data){
  data %>%
    ggplot(aes(regulation, percent, fill = regulation, group = contr_name)) +
    geom_col(width = 1, linewidth = 0.3) +
    geom_text(aes(label = paste0(num_genes, "\n(", percent, "%)"), y = percent/2),
              vjust = 0.5, hjust = 0.5, fontface = "bold", color = "#ffffff",
              size = 11, size.unit = "pt") +
    scale_fill_manual(name = NULL,
                      values = c(down = "#0000FF", up = "red3"),
                      labels = c("Downregulated", "Upregulated")) +
    coord_radial(expand = F) +
    facet_wrap(~contr_name, ncol = 3, scales = "free") +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          strip.text = element_text(size = 10, face = "bold",
                                    margin = margin(b  = -10)),
          strip.clip = "off",
          panel.spacing = unit(0, "pt"),
          strip.background = element_blank(),
          legend.position = "inside",
          legend.position.inside = c(0.7, 0.2),
          legend.background = element_rect(colour = "black", linewidth = 0.3,
                                           fill = NA),
          legend.key = element_rect(color = NA),
          legend.key.spacing.y = unit(10, "pt"),
          legend.key.height = unit(5, "pt"),
          legend.key.width = unit(10, "pt"),
          legend.text = element_text(size = 12, face = "bold",
                                     margin = margin(l = 3)))
}

pie1 <- list(piechart1a = plot_deg_pie1(plot_data[[1]]),
             piechart2a = plot_deg_pie1(plot_data[[2]]))

map2(.x = pie1, .y = names(pie1),
     ~ggsave(plot = .x, filename = paste0(
       result_path, "figures/DEG_levels_", .y, ".png"),
       width = 7.5, height = 6.8, dpi = 350))

plot_deg_pie2 <- function(data){
  data %>%
    ggplot(aes(x = 1, percent, fill = regulation)) +
    geom_col(width = 1) +
    geom_text(aes(label = paste0(num_genes, "\n(", percent, "%)"), y = percent/2),
              x = rep(c(0.05, 0.95), 7),
              vjust = 0.5, hjust = 0.5, fontface = "bold", color = "#ffffff",
              size = 11, size.unit = "pt") +
    coord_radial(expand = F, theta = "y") +
    facet_wrap(~contr_name) +
    scale_fill_manual(name = NULL,
                      values = c(down = "#0000FF", up = "red3"),
                      labels = c("Downregulated", "Upregulated")) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          strip.text = element_text(size = 10, face = "bold",
                                    margin = margin(b  = -10)),
          strip.clip = "off",
          strip.background = element_blank(),
          legend.position = "inside",
          legend.position.inside = c(0.7, 0.2),
          legend.background = element_rect(colour = "black", linewidth = 0.3,
                                           fill = NA),
          legend.key = element_rect(color = NA),
          legend.key.spacing.y = unit(10, "pt"),
          legend.key.height = unit(5, "pt"),
          legend.key.width = unit(10, "pt"),
          legend.text = element_text(size = 12, face = "bold",
                                     margin = margin(l = 3)))
}

pie2 <- list(piechart1b = plot_deg_pie1(plot_data[[1]]),
             piechart2b = plot_deg_pie1(plot_data[[2]]))

map2(.x = pie2, .y = names(pie2),
     ~ggsave(plot = .x, filename = paste0(
       result_path, "figures/DEG_levels_", .y, ".png"),
       width = 7.5, height = 6.8, dpi = 350))
