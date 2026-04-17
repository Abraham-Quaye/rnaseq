#!/usr/bin/env Rscript --vanilla

library(magrittr)
library(ggtext)
library(tidyverse)

# the data is located in the "results/tables" folder
result_path <- "~/bm_fn_rnaseq/results/r/" 

# write function to extract data needed for downstream analysis

sig_data <- read_csv(file = paste0(
  result_path, "tables/significant_FN_vs_BM_DEGs.csv"
  )) %>%
    select(-matches("(FN|BM)")) %>%
  mutate(regulation = case_when(log2FoldChange >= 0 ~ "up",
                                log2FoldChange < 0 ~ "down",
                                TRUE ~ NA_character_))

############################# DEG bar plot ##################
deg_bar_plt <- sig_data %>%
  summarise(num_genes = n(), .by = c(regulation)) %>% 
  ggplot(aes(regulation, num_genes, fill = regulation)) +
  geom_col(show.legend = F, color = "#000000", linewidth = 0.3) +
  geom_text(aes(label = num_genes),
            vjust = -0.5, fontface = "bold",
            size = 15, size.unit = "pt") +
  scale_fill_manual(values = c(down = "#0000FF", up = "red3"),
                    labels = c("Downregulated", "Upregulated")) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(breaks = c("down", "up"),
                   labels = c("Downregulated", "Upregulated")) +
  coord_cartesian(clip = "off") +
  labs(title = "Differentially Expressed Genes <br>of FN VS BM",
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




p2 <- sig_data %>%
  summarise(num_genes = n(), .by = c(regulation)) %>%
  mutate(total = sum(num_genes)) %>%
  mutate(percent = round((num_genes/total) * 100, 2)) %>% 
  ggplot(aes(regulation, percent, fill = regulation)) +
  geom_col(width = 1, linewidth = 0.3) +
  geom_text(aes(label = paste0(num_genes, "\n(", percent, "%)"), y = percent/2),
            vjust = 0.5, hjust = 0.5, fontface = "bold", color = "#ffffff",
            size = 15, size.unit = "pt") +
  scale_fill_manual(name = NULL,
                    values = c(down = "#0000FF", up = "red3"),
                    labels = c("Downregulated", "Upregulated")) +
  coord_radial(expand = F) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.5, 0.95),
        legend.background = element_blank(),
        legend.key = element_rect(color = NA),
        legend.key.spacing.x = unit(30, "pt"),
        legend.key.height = unit(5, "pt"),
        legend.key.width = unit(10, "pt"),
        legend.text = element_text(size = 12, face = "bold",
                                   margin = margin(l = 3)),
        legend.direction = "horizontal",
        legend.text.position = "right")

ggsave(plot = p2,
       filename = paste0(result_path, "figures/DEG_levels_piechart1.pdf"),
       width = 4, height = 4)

plot_data <- sig_data %>%
  summarise(num_genes = n(), .by = c(regulation)) %>%
  mutate(total = sum(num_genes)) %>%
  mutate(percent = round((num_genes/total) * 100, 2)) 

p3 <- plot_data %>% 
  ggplot(aes(x = 1, percent, fill = regulation)) +
  geom_col(width = 1) +
  annotate(geom = "text", x = 1, y = c(25, 75),
           label = c(paste0(plot_data$num_genes[[1]],
                          "\n(", plot_data$percent[[1]], "%)"),
                     paste0(plot_data$num_genes[[2]],
                            "\n(", plot_data$percent[[2]], "%)")),
           fontface = "bold", color = "#ffffff",
           size = 15, size.unit = "pt"
           ) +
  scale_fill_manual(name = NULL,
                    values = c(down = "#0000FF", up = "red3"),
                    labels = c("Downregulated", "Upregulated")) +
  coord_radial(expand = F, theta = "y") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.5, 0.95),
        legend.background = element_blank(),
        legend.key = element_rect(color = NA),
        legend.key.spacing.x = unit(30, "pt"),
        legend.key.height = unit(5, "pt"),
        legend.key.width = unit(10, "pt"),
        legend.text = element_text(size = 12, face = "bold",
                                   margin = margin(l = 3)),
        legend.direction = "horizontal",
        legend.text.position = "right")

ggsave(plot = p3,
       filename = paste0(result_path, "figures/DEG_levels_piechart2.pdf"),
       width = 4.5, height = 4.5)


