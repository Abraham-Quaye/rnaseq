#!/usr/bin/env Rscript

library(ggsci)
library(patchwork)

break_sentence <- function(str, n_words){
  pat_match <- paste0("(\\b[-/\\w]+(?:\\s+[-/\\w]+){",
                      n_words - 1, "})\\s+")

  out <- gsub(pattern = pat_match,
              "\\1\\\n", str, perl = T)

  return(out)
}

separate_feature <- function(df, feature){

  if(!is.null(df[["ONTOLOGY"]])){
    df <- df %>% dplyr::filter(ONTOLOGY == feature)
    go_term <- df[["ONTOLOGY"]] %>% base::unique(.)
    go <- T
  }else{go <- F}

  sub_df <- df %>%
    dplyr::select(Description, Count, GeneRatio, padj = p.adjust) %>%
    dplyr::arrange(padj) %>%
    mutate(padj = formatC(padj, digits = 2, format = "e") %>% as.numeric(),
           GeneRatio = map_dbl(GeneRatio, ~base::eval(base::parse(text = .x))),
           Description = map_chr(Description, ~break_sentence(.x, n_words = 3))) %>%
    slice_head(n = 25)

  if(go){
    sub_df <- sub_df %>%
      tibble::add_column(ontology = go_term, .before = "Description")
  }
  return(sub_df)
}

plot_enrichment <- function(df){

  if(!is.null(df[["ontology"]])){
  term <- df %>% pull(ontology) %>% unique(.)
  term_lab <- paste0("GO Term - " , term)
  }else{
    term_lab <- "KEGG Pathway"
  }

  df %>%
    ggplot(aes(GeneRatio, reorder(Description, GeneRatio),
               group = padj, color = padj, size = Count)) +
    geom_point() +
    theme_classic() +
    labs(y = term_lab, x = "Rich Factor") +
    # scale_y_continuous(breaks = c(seq(0, 1, 0.05)),
    #                    labels = c(seq(0, 1, 0.05))) +
    scale_color_gradient(low = "red", high = "blue", name = "Adjusted P-value",
                         # breaks = c(0, NA),
                         # labels = c(0, NA)
                         ) +
    scale_size(name = "Number of Genes", range = c(2, 10)) +
    theme(plot.margin = margin(4, 4, 4, 4),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          panel.grid.major = element_line(colour = "gray", linewidth = 0.2),
          axis.title = element_text(size = 18, face = "bold", colour = "#000000"),
          axis.text = element_text(size = 14, face = "bold", colour = "#000000"),
          axis.ticks.y = element_line(linewidth = 0.2, color = "gray"),
          axis.ticks.length.y = unit(8, "pt"),
          legend.background = element_rect(colour = "grey", linewidth = 0.15),
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 8, face = "bold"),
          legend.justification = c(0, 0),
          legend.position = "inside",
          legend.position.inside = c(0.8, 0.1))
}


tab_res <- tibble(category = c("all_bp", "all_cc", "all_mf",
                               "down_bp", "down_cc", "down_mf",
                               "all_kegg", "up_kegg", "down_kegg"
                               ),
                  term = str_replace(category, "[a-z]+_([a-z]+)", "\\1")) %>%
  dplyr::mutate(term = ifelse(term == "kegg", NA_character_, term) %>%
                  toupper(.),
                raw_data = list(all_go_all_table, all_go_all_table, all_go_all_table,
                                down_go_all_table, down_go_all_table, down_go_all_table,
                                all_kegg_tab, up_kegg_tab, down_kegg_tab),
                plt_rdy_results = map2(.x = raw_data, .y = term,
                                      ~separate_feature(.x, .y)),
                dotplots = map(plt_rdy_results, plot_enrichment),
                plt_title_spec = ifelse(is.na(term), "KEGG Pathways",
                                        paste0("GO-", term, " Terms")),
                dotplot_title = paste0("Most Significantly Enriched ",
                                       plt_title_spec, " at 72 hours"),
                titled_dotplots = map2(dotplots, dotplot_title,
                                      \(plt, ttl) plt + labs(title = ttl)),
                dotplot_name = paste0("results/figures/dotplot_",
                                     category, ".pdf"))

# save plots
map2(tab_res$titled_dotplots, tab_res$dotplot_name,
     ~ggsave(plot = .x, filename = here(.y),
             width = 12, height = 18))
