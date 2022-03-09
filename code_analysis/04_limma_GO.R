my_ont = "BP" 
my_showCategory = 10 
my_logFC_threshold = 0.15


i = 2

for (i in c(1:2)){
  
regions[i] -> my_region

# make big df with gene, brain region and qval
limma_list<- readRDS(glue("results_RNAseqRDS/limma_{my_region}_noIEGs.RDS")) %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) 


gettopGO <- function(limma_df,my_showCategory = 5){
  go_df_up <- limma_df %>% 
    filter(logFC>0) %>% 
    left_join(grcm38 %>% dplyr::select(symbol, entrez)) %>% 
    mutate(entrez = as.character(entrez)) %>% 
    dplyr::select(entrez, logFC) %>% 
    filter(!is.na(entrez)) %>% 
    arrange(desc(logFC))
  
  
  go_df_down <- limma_df %>% 
    filter(logFC<0) %>% 
    left_join(grcm38 %>% dplyr::select(symbol, entrez)) %>% 
    mutate(entrez = as.character(entrez)) %>% 
    dplyr::select(entrez, logFC) %>% 
    filter(!is.na(entrez)) %>% 
    arrange(desc(logFC))
  
  ggo_up <- enrichGO(gene = go_df_up$entrez %>% unique(),
                     OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                     keyType = "ENTREZID",
                     ont = my_ont,
                     readable = T,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.1,
                     qvalueCutoff  = 0.50)
  
  ggo_down <- enrichGO(gene = go_df_down$entrez %>% unique() ,
                       OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                       keyType = "ENTREZID",
                       ont = my_ont,
                       readable = T,
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.1, 
                       qvalueCutoff  = 0.50) 
  
  fortify(
    ggo_up,
    showCategory = my_showCategory,
    by = "Count",
    split = NULL,
    includeAll = TRUE
  ) %>% 
    arrange(desc(GeneRatio)) %>% 
    mutate(direction = "Up") -> temp1
  
  
  fortify(
    ggo_down,
    showCategory = my_showCategory,
    by = "Count",
    split = NULL,
    includeAll = TRUE
  ) %>% 
    arrange(desc(GeneRatio)) %>% 
    mutate(direction = "Down") -> temp2
  return(rbind(temp1,temp2))
  
}


limma_list %>% 
  map(gettopGO) -> GO_result_list

GO_result_list %>% 
  map(~select(., geneID, Description, direction)) %>% 
  map2_df(.,names(.), ~mutate(.x, contrast = .y)) -> GO_results

write.csv(GO_results, glue("results_tables/GO_results_{my_region}_noIEGs.csv"), row.names = F)

}



# get top genes
for (i in 1:2){
  my_region <- regions[i]
  
  limma_list<- readRDS(glue("results_RNAseqRDS/limma_{my_region}_noIEGs.RDS")) %>% 
    map(~distinct(.)) %>% 
    map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
    map(~filter(.,P.Value <0.05)) 
  
  #colnames(limma_list$status)
  limma_list %>% 
    map(~summarise(.,Up = sum(logFC>0),
                   Down = sum(logFC<0))) %>% 
    map(~mutate(.,Total = Up + Down)) %>% print
  
  limma_list %>% 
    map(~arrange(.,desc(logFC))) %>% 
    map(~head(.,8)) %>% 
    map(~.$symbol) %>%  print
  
  limma_list %>% 
    map(~arrange(.,logFC)) %>% 
    map(~head(.,9)) %>% 
    map(~.$symbol) %>%  print
  
}



