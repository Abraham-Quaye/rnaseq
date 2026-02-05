#!/usr/bin/env Rscript --vanilla

library(magrittr)
library(ggtext)
library(tidyverse)

# the data is located in the "results/tables" folder
result_path <- "~/myocd_rnaseq/results/r/"

# write function to extract data needed for downstream analysis

sig_data <- read_csv(file = paste0(
  result_path, "tables/significant_MYOCD_vs_GFP_DEGs.csv"
  )) %>%
    select(-matches("(MYOCD|GFP)")) %>%
  mutate(regulation = case_when(log2FoldChange >= 0 ~ "up",
                                log2FoldChange < 0 ~ "down",
                                TRUE ~ NA_character_))

############################# DEG bar plot ##################
deg_bar_plt <- sig_data %>%
  summarise(num_genes = n(),
                   .by = c(regulation)) %>% 
  ggplot(aes(regulation, num_genes, fill = regulation)) +
  geom_col(show.legend = F) +
  geom_text(aes(label = num_genes),
            vjust = -0.5, fontface = "bold",
            size = 15, size.unit = "pt") +
  scale_fill_manual(values = c(down = "#0000FF", up = "#ff0000"),
                    labels = c("Downregulated", "Upregulated")) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(breaks = c("down", "up"),
                   labels = c("Downregulated", "Upregulated")) +
  coord_cartesian(clip = "off") +
  labs(title = "Differentially Expressed Genes <br>of MYOCD VS GFP Controls",
       y = "Number of Genes", x = NULL) +
  theme_classic() +
  theme(plot.margin = margin_auto(4, 4),
        plot.title.position = "plot",
        plot.title = element_markdown(face = 'bold', size = 16,
                                      colour = 'black', hjust = 0.5,
                                      lineheight = 1.3, vjust = 0.5,
                                      margin = margin(b = 25, t = 3)),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major.y = element_line(color = 'grey50', linewidth = 0.1,
                                          linetype = 1),
        panel.grid.minor = element_line(linewidth = 0.1, colour = "grey50",
                                        linetype = "dashed"),
        panel.spacing.x = unit(5, "pt"),
        #adjust axis
        axis.text.x = element_markdown(size = 12, colour = 'black', face = 'bold',
                                       margin = margin(t = 10)),
        axis.text.y = element_text(size = 15, colour = 'black', face = 'bold',
                                   margin = margin(l = 10, r = 5)),
        axis.title.y = element_text(size = 18,
                                    face = 'bold',
                                    color = 'black'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = 'grey50', linewidth = 0.1,
                                    linetype = 1),
        axis.ticks.length.y = unit(5, "pt")
  )

ggsave(plot = deg_bar_plt,
       filename = paste0(result_path, "figures/DEG_levels_barplot.pdf"),
       width = 6, height = 7.5)
