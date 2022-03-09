

my_networkType = "signed hybrid"
my_region = "mPOA"
my_power= "4"
my_module = c("blue","brown","red", "pink", "magenta", "greenyellow", "lightcyan", "royalblue")
 

my_networkType = "signed hybrid"
my_region = "vHYP"
my_power= "4"
my_module = c("blue","green", "black")


lnames = load(glue("results_WGCNA/{my_region}-networkConstruction-auto_power{my_power}_noIEGs.RData"))
#The variable lnames contains the names of loaded variables.
lnames
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


for (i in 1:length(my_module)){
  
  mod <- my_module[i]
  glue("ME{mod}") -> modx
  
  
  
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  orderMEs(MEs0) %>% 
    rownames_to_column("sampleID") %>%
    left_join(sample_data %>% left_join(sample_data) %>% 
                select(sampleID,status, ds, cort_post, bw_change_rate) %>% 
                mutate_if(is.numeric,scale)) %>% 
    relocate(status, ds, cort_post, bw_change_rate) %>% 
    mutate(statusx = factor(status, levels = c("Subdominant","Alpha","Subordinate"))) %>% 
    relocate(status, statusx) %>% 
    select(-sampleID)-> ME_df
  
  
  ME_df[,c(modx,"status")] -> plotdf
  colnames(plotdf)[1] <- "ME"
  
  max(plotdf$ME)*1.4 -> my_max
  min(plotdf$ME) -> my_min
  plotdf %>% 
    ggplot(aes(status,ME, color = status, fill = status))+
    geom_boxjitter(outlier.color = NA, jitter.shape = 21, 
                   jitter.color = NA,
                   alpha = 0.3,
                   jitter.height = 0.02, jitter.width = 0.07, 
                   errorbar.draw = TRUE,
                   position = position_dodge(0.85)) +
    scale_color_manual(values = c("purple4", "#21908CFF","orange")) +
    scale_fill_manual(values = c("purple4", "#21908CFF","orange")) +
    theme_bw(base_size = 6)+
    theme(legend.position = "none",
          plot.title = element_text(size = 7)) +
    ylim(my_min, my_max)+
    labs(x = "",
         y = "Module eigengene",
         title = glue("{my_region}: {mod} module")) -> temp_p2
  
  png(filename = glue("results_figures/ME_{my_region}_{mod}_noIEGs.png"),
      width = 5, height = 4.5, units = "cm", res = 600)
  print(temp_p2)
  invisible(dev.off())
  
}

