# Preset
my_logFC_threshold = 0.15
i = 1
j = 1
my_comparison <- c("Alpha vs. Subdom", "Alpha vs. Sub", "Subdom vs. Sub")

neuroendo_mouse_geneset <-readRDS("results_RNAseqRDS/neuroendo_mouse_geneset.RDS")

# for now, mPOA and VMH
limma_list <-list()
for(i in c(1,2)){
  my_region <- regions[i]
  limma_list[[my_region]] <- readRDS(glue("results_RNAseqRDS/limma_{my_region}_noIEGs.RDS"))
}


# 1. Neuroendorcine gene set = Donkelaar 2020
neuroendo_mouse_geneset %>% 
  as_tibble() %>% 
  .$symbol -> my_gs_sets
my_gs_names <- c("Neuroendocrine gene set")


# now plot it 

compdat = data.frame(x = c("Alpha","Alpha","Subdom"),
                     y = c("Subdom","Sub","Sub"))


ARG_list$rapid_PRGs %>% 
  rename(symbol = `Gene ID`) %>% 
  filter(!is.na(symbol)) %>% 
  .$symbol -> IEG_symbol

for (i in 1:length(limma_list)){
  
  for(j in 1:length(limma_list[[i]])){
    

    
    
    limma_list[[i]][[j]]  %>% 
      filter(symbol %in% my_gs_sets) %>% 
      filter(symbol %notin% IEG_symbol) %>% 
      dplyr::select(symbol, logFC, P.Value) %>% 
      dplyr::mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                                 ifelse(logFC>my_logFC_threshold,"Dominant genes","Subordinate genes"))) %>%
      dplyr::mutate(Sig = factor(Sig, levels = c("Dominant genes","Subordinate genes", "N.S."))) %>% 
      dplyr::mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>% 
      unique() -> df
    
    
    keyvals <- ifelse(
      df$logFC < -my_logFC_threshold & df$P.Value<0.05, 'orange',
      ifelse(df$logFC > my_logFC_threshold & df$P.Value<0.05, 'purple4',
             'grey'))
    
    table(keyvals) %>% 
      as.data.frame() %>% 
      mutate(prop = glue("{round(Freq/sum(Freq),3)*100}")) %>% 
      mutate(text = glue("{Freq} genes ({prop}%)")) -> my_texts
    
    if(is.na(my_texts[2,4])){my_texts[2,4] <- "0 genes"}
    if(is.na(my_texts[3,4])){my_texts[3,4] <- "0 genes"}
    
    keyvals[is.na(keyvals)] <- 'grey'
    names(keyvals)[keyvals == 'purple4'] <- 'high'
    names(keyvals)[keyvals == 'grey'] <- 'mid'
    names(keyvals)[keyvals == 'orange'] <- 'low'
    
    
    df %>% 
      filter(P.Value <0.05) %>%
      filter(abs(logFC) > my_logFC_threshold) %>%
      arrange(logFC) %>%
      mutate(effects = abs(sin(logFC)*-log10(P.Value))) %>%
      top_n(40) %>% 
      # top_frac(.,1) %>%
      # filter(P.Value <0.01 | abs(logFC) >0.5) %>%
      .$symbol -> for_label
    
    png(filename = glue("results_figures/EnVol/{my_gs_names}_{names(limma_list)[i]}_{names(limma_list[[i]][j])}_noIEGs.png"),
        width = 8, height = 7.3, units = "cm", res = 600)
    
    
    EnhancedVolcano(df,
                    lab = df$symbol,
                    selectLab= for_label,
                    x = 'logFC',
                    y = 'P.Value',
                    title = glue("{my_gs_names}: {names(limma_list)[i]}: {my_comparison[j]}"),
                    pCutoff = 0.05,
                    FCcutoff = 0.0,
                    cutoffLineType = 'blank',
                    # drawConnectors = TRUE, # just for a few plots
                    # widthConnectors = 0.1, # just for a few plots
                    vline = c(-0.15, 0.15),
                    vlineCol = c('grey90'),
                    vlineType = c( 'dashed'),
                    vlineWidth = c(0.3),
                    hline = c(0.05),
                    hlineCol = c('grey90'),
                    hlineType = c( 'dashed'),
                    hlineWidth = c(0.3),
                    colCustom = keyvals,
                    shape = 21,
                    pointSize = 2.5,
                    labCol = "black",
                    labSize = 2.1)+
      annotate("text", x = 0.8, y = 4.4,color = "purple4" , size = 3,
               label = glue("Upregulated in {compdat[j,1]} \n{my_texts[3,4]}"))+
      annotate("text", x = -0.8, y = 4.4,color = "orange" , size = 3,
               label = glue("Upregulated in {compdat[j,2]} \n{my_texts[2,4]}"))+
      scale_x_continuous(limits = c(-1.5,1.5),breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5))+
      scale_y_continuous(limits = c(-0.1,4.8),breaks = c(0,1,2,3,4),expand=expansion(mult=c(0.0005,0.0)))+
      theme_bw(base_size = 7)+
      labs(color = "",
           caption = paste0('total = ', nrow(df), ' genes'),
           y = bquote(~-Log[10]~italic(eFDR)))+
      theme(legend.position = "none",
            plot.title = element_text(size = 7),
            plot.subtitle = element_blank()) -> temp
    print(temp)
    invisible(dev.off())
    
    
    
  }
  
}
