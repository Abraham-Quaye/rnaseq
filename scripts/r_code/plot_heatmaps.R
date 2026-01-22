#!/usr/bin/env Rscript

library(pheatmap)
library(patchwork)
library(ggplotify)

## select columns to use for making heatmap
heat_data <- sig_72_data %>%
  arrange(qval) %>%
  arrange(desc(abs(log2fc))) %>%
  slice_head(n = 80) %>%
  arrange(regulation) %>%
  dplyr::select(entrez:symbol, wt_1:r77q_3) %>%
  distinct(symbol, .keep_all = T)

r_names <- heat_data %>%
  pull(symbol)


c_names <- c("WT 72hr1", "WT 72hr2", "WT 72hr3", "R77Q 72hr1",
             "R77Q 72hr2", "R77Q 72hr3")


heat_matrix <- heat_data %>%
  dplyr::select(wt_1:r77q_3) %>%
  data.matrix(.) %>%
  set_rownames(r_names) %>%
  set_colnames(c_names) %>%
  t(.) %>%
  scale(.) %>%
  t(.)

heat72 <- pheatmap(heat_matrix,
         main = "Top 80 Significant DEGs at 72 hours Scaled by Z-Scores",
         cellwidth = 50,
         cellheight = 10,
         color = colorRampPalette(
           colors = c('blue','white','red'))(250),
         angle_col = 45,
         fontsize = 15,
         fontsize_row = 10,
         fontsize_col = 15,
         treeheight_row = 0
         )

# save heatmaps in a tibble with appropriate identifier (timepoint)
ggsave(plot = heat72, filename = here("results/figures/heatmap72.pdf"),
       width = 9.5, height = 10.5)

