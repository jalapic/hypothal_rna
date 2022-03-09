my_region = "mPOA"
my_region = "vHYP"


my_ont = "BP"
my_showCategory = 10
my_power = 4 

#module<- allcolors[1]  #for test

gettop10GO_WGCNA <- function(module,my_showCategory = 10){
  
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  colnames(datExpr)[moduleColors==module] -> module_gene
  
  
  grcm38 %>% 
    dplyr::filter(ensgene %in% module_gene) %>% 
    filter(!is.na(entrez)) %>% 
    dplyr::select(entrez) -> go_df_wgcna
  
  
  ggo <- enrichGO(gene = go_df_wgcna$entrez %>% unique(),
                  OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                  keyType = "ENTREZID",
                  ont = my_ont,
                  readable = T,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.2,
                  qvalueCutoff  = 0.50)
  
  
  fortify(
    ggo,
    showCategory = my_showCategory,
    by = "Count",
    split = NULL,
    includeAll = TRUE
  ) %>% 
    arrange(desc(GeneRatio)) %>% 
    mutate(module = module) -> temp1
  
  return(rbind(temp1))
  
}



for(r in 1:c(length(regions))){   
  
  regions[r] -> my_region
  
  
  lnames = load(glue("results_WGCNA/{my_region}-networkConstruction-auto_power{my_power}_noIEGs.RData"))
  #The variable lnames contains the names of loaded variables.
  lnames

  modNames = substring(names(MEs), 3)
  
  moduleColors %>% unique() -> allcolors
  
  WGCNA_GOs <- vector('list', length(allcolors))

  for(i in 1:length(allcolors)){
    gettop10GO_WGCNA(allcolors[i],my_showCategory) -> WGCNA_GOs[[i]]
  }
  WGCNA_GOs %>% 
    rbindlist() -> wgcna_all_gos
  
    write.csv(wgcna_all_gos,glue("results_tables/wgcna_all_gos_noIEGs_{my_region}.csv"),row.names = F)  
}

