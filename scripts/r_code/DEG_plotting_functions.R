
plot_volcano <- function(lfc_res_tbl, treatment){
  
    # extract sample names for plot title
    samplenames <- str_split(treatment, "_") %>%
        base::unlist(.) %>% .[c(2,4)]
    
    plt_title <- "MYOCD VS GFP Control"
    
    # use lfcShrink results for volcano plot
    dff <- lfc_res_tbl %>%
        drop_na(padj) %>%
        mutate(sig = padj <= 0.05,
               reg = case_when(sig & log2FoldChange > 0 ~ "up",
                               sig & log2FoldChange < 0 ~ "down",
                               TRUE ~ "normal"))
    
    most_sig <- dff %>%
      drop_na(SYMBOL) %>%
      filter(sig & (padj >= quantile(padj, probs = 0.99) |
                      abs(log2FoldChange) >= quantile(abs(log2FoldChange),
                                                      probs = 0.99))) %>%
      dplyr::arrange(padj, log2FoldChange)
    
    if(nrow(most_sig) > 20){
        most_sig <- most_sig %>%
            slice_head(n = 20)
    }
    
    plt_bounds <- dff %>% summarise(xmin = min(log2FoldChange, na.rm = T),
                                    xmax = max(log2FoldChange, na.rm = T))
    
    plt <- dff %>%
        ggplot(aes(log2FoldChange, -log10(padj), colour = reg)) +
        geom_point(alpha = 0.7, size = 2.5) +
        geom_text_repel(data = most_sig,
                        aes(log2FoldChange, -log10(padj),
                            label = SYMBOL), max.overlaps = 40,
                        min.segment.length = 0, fontface = "bold",
                        box.padding = 1.5, point.padding = 0.5,
                        show.legend = F, size = 3.5) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
        scale_y_continuous(expand = c(0.025, 0.025)) +
        scale_x_continuous(expand = c(0, 0),
                           breaks = seq(from = floor(plt_bounds$xmin),
                                        to = ceiling(plt_bounds$xmax), by = 1),
                           labels = seq(from = floor(plt_bounds$xmin),
                                        to = ceiling(plt_bounds$xmax), by = 1),
                           limits = c((plt_bounds$xmin - 0.2), (plt_bounds$xmax + 0.2))) +
        coord_cartesian(clip = "off") +
        labs(
            title = plt_title,
            y = expression("-Log"[10]*"(P-adjusted Values)"),
            x = expression("Log"[2]*"(Fold Change)")) +
        scale_color_manual(values = c(down = "blue", normal = "grey", up = "red"),
                           labels = c("Downregulated", "Not Significant", "Upregulated")) +
        theme_bw() +
        theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
              plot.title.position = "plot",
              axis.title = element_text(size = 18, face = "bold"),
              axis.text = element_text(size = 16, colour = "#000000"),
              legend.title = element_blank(),
              legend.text = element_text(size = 15, face = "bold",
                                         margin = margin(r = 30, l = 0)),
              legend.position = "top",
              legend.direction = "horizontal",
              legend.key.size = unit(0.7, "cm"),
              legend.key = element_rect(fill = NA)) +
        guides(color = guide_legend(override.aes = list(size = 7)))
    
    return(plt)
}

## Distance Matrix
plot_sample_dists <- function(dds, dds_design, color_grp_feature, row_labs_feature){
    
    # rlog or vst: to stabilize variance across samples:
    stable_var <- rlog(dds, blind = T)
    dists <- dist(t(assay(stable_var)))
    dist_mat <- as.matrix(dists)
    
    col_heat <- colorRampPalette(hcl.colors(n = 12, palette = "Blues"), alpha = 1)(255)
    
    col_grp_labs <- base::as.data.frame(colData(dds)[, color_grp_feature, drop = F])
    
    col_grp_features <- pull(col_grp_labs, color_grp_feature) %>% base::unique(.)
    
    ann_colors <- setNames(
        object = list(c("blue", "red") %>%
                          set_names(col_grp_features)),
        nm = color_grp_feature
    )
    
    readable_labels <- base::as.data.frame(colData(dds)[, row_labs_feature,
                                                        drop = F]) %>%
        pull(condition) %>% toupper(.)
    
    hmap <- pheatmap(dist_mat,
                     annotation_col =  col_grp_labs,
                     labels_row = readable_labels,
                     labels_col = readable_labels,
                     show_colnames = T,
                     show_rownames = T,
                     annotation_colors = ann_colors,
                     clustering_distance_cols = dists,
                     clustering_distance_rows = dists,
                     color = col_heat,
                     treeheight_row = 30,
                     treeheight_col = 30,
                     fontsize = 12,
                     cellwidth = 20,
                     cellheight = 20)
    
    return(as.ggplot(hmap))
}

## PCA plot
plot_PCA <- function(dds, dds_design){

  # rlog or to stabilize variance across samples:
  stable_var <- rlog(dds, blind = F)
  pcaData <- plotPCA(stable_var,
                     intgroup = dds_design,
                     returnData = T,
                     ntop = base::nrow(dds))
  pvar <- round(attr(pcaData, "percentVar") * 100)

  plt <- pcaData %>%
    ggplot(aes(PC1, PC2, fill = group)) +
    geom_point(size = 4, shape = 21, color = "#000000", stroke = 0.2) +
    coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("blue", "red"),
                       name = str_to_title(dds_design),
                       breaks = levels(pcaData[[dds_design]]),
                      labels = toupper(levels(pcaData[[dds_design]]))
                         ) +
    labs(x = paste0("PC1: ", pvar[[1]], "% Variance"),
         y = paste0("PC2: ", pvar[[2]], "% Variance")) +
    theme_bw() +
    theme(plot.title = element_text(size = 18.5, face = "bold", hjust = 0),
          plot.title.position = "plot",
          panel.grid.major.y = element_line(color = "#000000", linewidth = 0.2,
                                            linetype = "dashed"),
          axis.title = element_text(size = 18, face = "bold"),
          axis.text = element_text(size = 16, colour = "#000000"),
          legend.title = element_blank(),
          legend.text = element_text(size = 12, face = "bold",
                                     margin = margin(r = 5, l = -5)),
          legend.position = "top",
          legend.direction = "horizontal",
          legend.key.size = unit(1, "cm"),
          legend.margin = margin(b = -10),
          legend.key.spacing.y = unit(-10, "pt"))
  
  return(plt)
}


## heatmap for top genes in all samples ===============================
# The assay function is used to extract the matrix of normalized values.
# it gives you: Regularized log-transformed counts using the rlog function, which stabilizes variance across the range of expression values.
# Values: On a log2 scale, with low-count genes shrunk toward the mean to reduce noise.
# Use case: Ideal for visualization (e.g., PCA, heatmaps), clustering, and exploratory analysis.
# If you're building heatmaps or PCA plots, go with rlog or vst. If you're doing statistical testing or exporting count tables, use the normalized counts.
plot_topGenes_heatmap <- function(sig_results, contr_name){
  
  if(NROW(sig_results) > 50){
    n_topgenes <- 50
  }else if(NROW(sig_results) >= 5){
    n_topgenes <- NROW(sig_results)
  }else{return(NULL)}
  # Construct the plot title
  feature_remove <- contr_name %>% str_split(., "_") %>%
    flatten() %>% pluck(1)
  
  title_name <- contr_name %>% 
    str_remove(., paste0(feature_remove, "_")) %>%
    toupper(.) %>%
    str_replace_all(., "HR", "hr") %>%
    paste0("Top", n_topgenes, " DEGs for ", .,
           "\nby Adjusted P-values and Log2(Fold Change)")
  
  # Prepare Matrix to plot
  sample_cols <- contr_name %>%
    str_remove(., paste0(feature_remove, "_")) %>%
    str_split(., "_vs_") %>%
    base::unlist()
  
  top_genes <- sig_results %>%
    arrange(padj, log2FoldChange) %>%
    drop_na(SYMBOL) %>%
    mutate(isduplicate = duplicated(SYMBOL)) %>% 
    filter(!isduplicate) %>% 
    select(SYMBOL, starts_with(sample_cols)) %>%
    as.data.frame() %>%
    set_rownames(.$SYMBOL) %>%
    select(-SYMBOL) %>%
    rename_all(toupper)
  
  top_genes_mat <-  top_genes[apply(top_genes, 1, var) != 0, ] %>% 
    head(., n_topgenes) %>%
    as.matrix()
  
  hmap <- pheatmap(top_genes_mat,
                   main = title_name,
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   cellwidth = 30,
                   show_rownames = TRUE,
                   show_colnames = TRUE,
                   color = colorRampPalette(
                     colors = c('blue','white','red'))(250),
                   scale = "row")
  
  return(as.ggplot(hmap))
}
