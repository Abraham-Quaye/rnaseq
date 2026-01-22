#!/usr/bin/env Rscript --vanilla

library(magrittr)
library(ggtext)
library(tidyverse)

# the data is located in the "results/tables" folder
result_path <- "~/berges_rnaseq/results/r/"

# get all file paths
deg_files <- list.files(paste0(result_path, "tables"),
                        pattern = "^significant_\\w+\\d{1,2}hr_DEGs\\.csv",
                        full.names = TRUE)

# write function to extract data needed for downstream analysis
get_deg_tables <- function(csv){
  df <- read_csv(file = csv, id = "contr_name") %>%
    select(-matches("(r77q|wt|mock)"), -stat) %>%
    mutate(contr_name = str_replace(contr_name,
                                    "[/\\w]+significant_(\\w+)_DEGs\\.csv",
                                    "\\1"))
  if(str_detect(df$contr_name[[1]], "mock")){
    df <- df %>%
      mutate(timepoint = str_remove(contr_name, "_vs_mock72hr") %>%
               str_extract(., "\\d{1,2}hr") %>%
             parse_number(.))
  }else{
    df <- df %>%
      mutate(timepoint = str_extract(contr_name, "\\d{1,2}hr") %>%
               parse_number(.))
  }
  return(df)
}

# read all data into one big dataframe
deg_data <- map(deg_files, ~get_deg_tables(.x)) %>%
  list_rbind() %>%
  mutate(regulation = case_when(log2FoldChange >= 0 ~ "up",
                                log2FoldChange < 0 ~ "down",
                                TRUE ~ NA_character_),
         ref_sample = ifelse(str_detect(contr_name, "mock"),
                             "All Samples VS Mock", "R77Q VS WT"))

############################# DEG bar plot ##################
deg_bar_plt <- deg_data %>%
  summarise(num_genes = n(),
                   .by = c(ref_sample, regulation, timepoint)) %>%
  # this fills in combinations with zero counts
  complete(data = ., ref_sample, regulation, timepoint,
           fill = list(num_genes = 0)) %>%
  mutate(timepoint = factor(timepoint, levels = sort(unique(timepoint)))) %>%
  ggplot(aes(timepoint, num_genes, fill = regulation)) +
  geom_col(position = position_dodge(0.9)) +
  facet_wrap(~ref_sample, scales = "free") +
  geom_text(aes(label = num_genes), position = position_dodge(width = 0.9),
            vjust = -0.5, fontface = "bold", size = 4) +
  scale_fill_manual(values = c("#0000FF", "#ff0000"),
                    breaks = c("down", "up"),
                    labels = c("Downregulated", "Upregulated")) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0.16, 0.16),
                   breaks = c(4, 8, 12, 24, 72),
                   labels = c("4 hours", "8 hours", "12 hours",
                              "24 hours", "72 hours")) +
  coord_cartesian(clip = "off") +
  labs(title = "Differentially Expressed Genes of R77Q HIV-infected Cells<br>Compared to Mock- or Wild Type HIV-infected Cells",
       y = "Number of Genes",
       x = NULL,
       fill = NULL) +
  theme_classic() +
  theme(plot.margin = margin_auto(4, 4),
        plot.title.position = "plot",
        plot.title = element_markdown(face = 'bold', size = 18,
                                      colour = 'black', hjust = 0.5,
                                      lineheight = 1.3, vjust = 0.5),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major.y = element_line(color = 'grey50',
                                          linewidth = 0.035,
                                          linetype = 1),
        panel.grid.minor = element_line(linewidth = 0.1, colour = "grey50",
                                        linetype = "dashed"),
        panel.spacing.x = unit(5, "pt"),
        # adjust the facet strip
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x.top = element_text(face = "bold", size = 15,
                                        colour = "black",
                                        margin = margin(b = 15, t = 15)),
        strip.clip = "off",
        #adjust axis
        axis.text.x = element_markdown(size = 12, colour = 'black', face = 'bold',
                                       margin = margin(t = 10)),
        axis.text.y = element_text(size = 15, colour = 'black', face = 'bold',
                                   margin = margin(l = 10)),
        axis.title.y = element_text(size = 18,
                                    face = 'bold',
                                    color = 'black'),
        axis.ticks = element_blank(),
        legend.justification = c(0,1),
        legend.position = "inside",
        legend.position.inside = c(0.05, 0.925),
        legend.box.background = element_rect(color = "grey60"),
        legend.text = element_text(face = 'bold', size = 12),
        legend.background = element_blank(),
        legend.margin = margin(t = 3, r = 3, b = 3, l = 3),
        legend.key.size = unit(0.8, "cm"),
        legend.key.spacing.y = unit(0.2, "cm")
  )

ggsave(plot = deg_bar_plt,
       filename = paste0(result_path, "figures/DEG_levels_barplot.pdf"),
       width = 11, height = 8)
