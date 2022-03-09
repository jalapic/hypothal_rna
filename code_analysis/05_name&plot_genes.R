# how many genes are significant? ####

limma_list <-list()
for(i in 1: length(regions)){
  my_region <- regions[i]
  limma_list[[my_region]] <- readRDS(glue("results_RNAseqRDS/limma_{my_region}_noIEGs.RDS")) %>% 
    map(distinct)
}
names(limma_list)

my_region = "mPOA" 
my_region = "vHYP"
my_logFC_threshold = 0.15

limma_list[[my_region]] %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(distinct) %>% 
  map(~filter(.,P.Value <0.05))  -> limma_list_sig


limma_list_sig %>% 
  map(~mutate(.,updown = ifelse(logFC>0, "Up","Down"))) %>% 
  map(~group_by(.,updown)) %>% 
  map(count)

limma_list_sig %>% 
  map(~arrange(.,desc(logFC))) %>% 
  map(~head(.,10))

limma_list_sig %>% 
  map(~arrange(.,logFC)) %>% 
  map(~head(.,10))



# plot gene expression ####
# functions are in code_functions/functions_rnaseq.r

my_region = "mPOA"
my_gene_list <- c("Npy", "Avp", "Oxt", "Nos1")

my_region = "vHYP"
my_gene_list <- c("Npy", "Agrp", "Pomc", "Gal")

for(i in 1:length(my_gene_list)){
  my_gene_list[i] -> gene
  png(filename = glue("results_figures/exp_{my_region}_{gene}_new.png"),
      width = 7.5, height = 6, units = "cm", res = 600)
  plot_any_gene(gene,my_region)
  invisible(dev.off())
  
}


i = 1
i = 2
i = 3
i = 4
name_any_gene(my_gene_list[i],my_region)

  



