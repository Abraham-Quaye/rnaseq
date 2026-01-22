######### Extract data needed for downstream analysis
load_deg_results <- function(csv){
  df <- read_csv(file = csv, id = "contr_name") %>%
    select(-matches("(r77q|wt|mock)"), -stat, -pvalue) %>%
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

######## Subset significant results for just up or down regulated genes
get_subset_genes <- function(genes_tbl, reg){
  genes_tbl %>%
    drop_na(ENTREZID) %>%
    filter(regulation == reg) %>%
    distinct(ENTREZID, .keep_all = T) %>%
    select(-regulation)
}

####################### For GO Analysis #######################
get_go_enrich <- function(genes_tbl, ontology){
  enrichGO(gene = genes_tbl$gene_id,
           keyType = "ENSEMBL",
           OrgDb = org.Hs.eg.db,
           ont = ontology,
           # universe = raw_data$gene_id,
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05,
           readable = T)
}

####### Dotplots ==============================
plot_dotplot <- function(res, labb){
  if(is.null(res) | nrow(as_tibble(res)) == 0){return(NULL)
    }else if(nrow(res) > 25){
    num_cat <- 25
  }else{
    num_cat <- nrow(res)
  }
  
  dotplot(
    object = res,
    showCategory = num_cat,
    title = paste0("Top ", num_cat, " Enriched for ",
                   str_replace_all(str_replace_all(toupper(labb), "HR", "hr"),
                            "_", " ")),
    font.size = 10.5
  ) +
    theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5))
}

save_kegg_results <- function(res, contr_name, labb){
  map2(.x = res, .y = contr_name,
     \(.x, .y){
       as_tibble(.x) %>%
         write.csv(., file = paste0(
           result_path, "tables/kegg_",
           .y, labb, "DEG_sigPathways.csv"),
           row.names = F)
       }
     )
}

save_go_results <- function(res, contr_name, labb){
  map2(.x = res, .y = contr_name,
     \(.x, .y){
       as_tibble(.x) %>%
         write.csv(., file = paste0(
           result_path, "tables/go_",
           .y, labb, "DEG_sig.csv"),
           row.names = F)
       }
     )
}

save_kegg_dotplots <- function(dotplots, contr_name, labb){
  map2(.x = dotplots, .y = contr_name,
     \(.x, .y){
       if(is.null(.x)){
         return(NULL)
       }
       ggsave(plot = .x,
             filename = paste0(
               result_path, "figures/kegg_", .y,
               labb, "DEG_sigPathways_dotplot.pdf"),
             width = 10, height = 10)
       }
     )
}

save_go_dotplots <- function(dotplots, contr_name, labb){
  map2(.x = dotplots, .y = contr_name,
     \(.x, .y){
       if(is.null(.x)){
         return(NULL)
       }
       ggsave(plot = .x,
             filename = paste0(
               result_path, "figures/go_", .y,
               labb, "DEG_sig_dotplot.pdf"),
             width = 10, height = 10)
       }
     )
}

############## Function to plot KEGG Pathway Diagrams #############
plot_kegg_pathway <- function(id, gene_list, path_){
  pathview(gene.data = gene_list,
           species = "hsa",
           pathway.id = id,
           res = 500,
           kegg.dir = path_,
           limit = list(gene = 2, cpd = 2),
           low = list(gene = "blue", cpd = "blue"),
           mid = list(gene = "grey", cpd = "grey"),
           high = list(gene = "red", cpd = "red"),
           map.null = T)
  return("Complete")
}
safe_kegg_plotter <- safely(plot_kegg_pathway)
