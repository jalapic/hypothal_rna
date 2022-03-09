get_DEG_results <- function(dds, LFC_threshold = 0.2, pvalue_threshold = 0.05) # I won't set padj yet 
{
  
  temp <- results(dds, contrast=c("domgroup","Alpha","Subordinate")) %>% 
    as.data.frame() %>% 
    rownames_to_column('ensgene') %>% 
    as_tibble() %>% 
    select(ensgene, log2FoldChange, pvalue, padj) %>% 
    mutate(contrast = "Alpha - Subordinate") %>% 
    filter(abs(log2FoldChange)>=LFC_threshold) %>% 
    filter(pvalue <=pvalue_threshold) 
  
  
  temp %>% 
    arrange(desc(abs(log2FoldChange))) -> result
  
  return(result)
}


# name the gene 


limma_deg_table <- function(gene = "Avp", region = "vHYP"){
  limma_list[[region]] %>% 
    map(~filter(.,symbol == gene)) %>% 
    map2_df(., names(.), ~mutate(.x, comparison = .y)) %>% 
    select(-symbol, - adj.P.Val) %>% 
    relocate(comparison) %>% 
    mutate(Significant = ifelse(P.Value < 0.05, "YES", "")) %>% 
    knitr::kable(., "simple")
}

name_any_gene <-  function (gene = "Avp", region = "vHYP"){
  grcm38 %>% 
    filter(symbol == gene) %>% 
    .$description -> desc
  strrep("-",str_length(glue("{gene}: {desc}"))) -> mycat # I am making a lot of effort for trivial aesthetic detail...
  # gene info 
  cat(glue("{gene}: {desc} \n"+glue("{mycat}\n"+glue("in {region}"))))
  
  limma_deg_table (gene, region)
}



name_any_gene_multiregion <- function(gene = "Avp", 
                                      regions_of_interest = regions){
  grcm38 %>% 
    filter(symbol == gene) %>% 
    .$description -> desc
  strrep("-",str_length(glue("{gene}: {desc}"))) -> mycat # I am making a lot of effort for trivial aesthetic detail...
  # gene info 
  cat(glue("{gene}: {desc} \n"+glue("-\n"+glue(" "))))
  
  temp_p_list <- list()
  for(i in 1:length(regions_of_interest)){
    cat(glue("------------{mycat}  \n in {regions_of_interest[i]} \n "))
    print(limma_deg_table (gene, regions_of_interest[i]))
  }
}


# plot the expression 

limma_exp_plot <- function(gene = "Avp", region = "vHYP"){
  grcm38 %>% 
    filter(symbol == gene) %>% 
    .$ensgene %>% 
    unique() -> goi 
  
  if (goi %in% rownames(y_list[[region]]$E) == F)
  {temp_p <- ggplot()+ggtitle(glue("No detectable gene expression of {gene} in {region}"))}
  else {
    
    suppressMessages(
      y_list[[region]]$E[goi,] %>%
        as.data.frame(col.names = c("Expression")) %>% 
        rename(Expression = ".") %>%  
        rownames_to_column("sampleID") %>% 
        left_join(y_list[[region]]$target %>% 
                    rownames_to_column("sampleID") %>% 
                    select(sampleID, group)) -> temp)
        max(temp$Expression)*1.05 -> my_max
      min(temp$Expression) -> my_min
      temp %>% 
        ggplot(aes(group,Expression, color = group, fill = group)) + 
        geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                       alpha = 0.3,
                       jitter.size = 2,
                       jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                       position = position_dodge(0.85)) +
        ylim(my_min, my_max)+
        scale_color_manual(values = c("purple4","#21908CFF","orange"))+
        scale_fill_manual(values = c("purple4","#21908CFF","orange"))+
        labs(title = glue("{region}: {gene}"),
             x = "",
             y = "Normalized gene expression")+
        theme_bw(base_size = 10)+        ## this changes font size
        theme(legend.position = "none")  
    
  }
}


limma_exp_exclude <- function(gene = "Avp", region = "vHYP", exclude = "Subdominant"){
  grcm38 %>% 
    filter(symbol == gene) %>% 
    .$ensgene %>% 
    unique() -> goi 
  
  if (goi %in% rownames(y_list[[region]]$E) == F)
  {temp_p <- ggplot()+ggtitle(glue("No detectable gene expression of {gene} in {region}"))}
  else {
    
    suppressMessages(
      temp <- y_list[[region]]$E[goi,] %>%
        as.data.frame(col.names = c("Expression")) %>% 
        rename(Expression = ".") %>%  
        rownames_to_column("sampleID") %>% 
        left_join(y_list[[region]]$target %>% 
                    rownames_to_column("sampleID") %>% 
                    select(sampleID, group)) %>% filter(group != "Subdominant"))
      
      max(temp$Expression)*1.4 -> my_max
      min(temp$Expression) -> my_min
      temp %>% 
        ggplot(aes(group,Expression, color = group, fill = group)) + 
        geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                       alpha = 0.3,
                       jitter.size = 2,
                       jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                       position = position_dodge(0.85)) +
        scale_color_manual(values = c("purple4","orange"))+
        scale_fill_manual(values = c("purple4","orange"))+
        ylim(my_min,my_max)+
        labs(title = glue("{region}: {gene}"),
             x = "",
             y = "Normalized gene expression")+
        theme_bw(base_size = 6)+
        theme(legend.position = "none")  
    
  }
}

plot_any_gene <- function (gene = "Avp", region = "vHYP", exclude = F){
  grcm38 %>% 
    filter(symbol == gene) %>% 
    .$description -> desc
  strrep("-",str_length(glue("{gene}: {desc}"))) -> mycat # I am making a lot of effort for trivial aesthetic detail...
  # gene info 
  cat(glue("{gene}: {desc} \n"+glue("{mycat}\n"+glue("in {region}"))))
  if(exclude == F){
    limma_exp_plot(gene, region) -> temp_p
  }else{limma_exp_exclude(gene, region) -> temp_p}
  
  print(temp_p)
}

plot_any_gene_multiregion <- function(gene = "Avp", 
                                      regions_of_interest = regions,
                                      ncol = 2){
  grcm38 %>% 
    filter(symbol == gene) %>% 
    .$description -> desc
  strrep("-",str_length(glue("{gene}: {desc}"))) -> mycat # I am making a lot of effort for trivial aesthetic detail...
  # gene info 
  cat(glue("{gene}: {desc} \n"+glue("{mycat}\n"+glue(" in "))))
  
  temp_p_list <- list()
  for(i in 1:length(regions_of_interest)){
    cat(glue("\n {regions_of_interest[i]} \n "))
    temp_p_list[[i]] <- (limma_exp_plot(gene,regions_of_interest[i]))
  }
  
  grid.arrange(grobs = temp_p_list[1:length(regions_of_interest)],
               ncol = ncol)
}



# by CORT ===============================================================

limma_cort_plot <- function(gene = "Avp", region = "vHYP"){
  grcm38 %>% 
    filter(symbol == gene) %>% 
    .$ensgene %>% 
    unique() -> goi 
  
  if (goi %in% rownames(y_list[[region]]$E) == F)
  {temp_p <- ggplot()+ggtitle(glue("No detectable gene expression of {gene} in {region}"))}
  else {
    
    suppressMessages(
      y_list[[region]]$E[goi,] %>%
        as.data.frame(col.names = c("Expression")) %>% 
        rename(Expression = ".") %>%  
        rownames_to_column("sampleID") %>% 
        left_join(y_list[[region]]$target %>% 
                    rownames_to_column("sampleID") %>% 
                    select(sampleID, group)) %>% 
        left_join(sample_id) %>% 
        left_join(sample_data) %>%
        ggplot(aes(cort_post,Expression, color = group, fill = group)) + 
        geom_smooth(method = 'lm', alpha = 0.1)+
        geom_point(shape = 21, alpha = 0.3, size = 2.5)+
        scale_color_manual(values = c("purple4","#21908CFF","orange"))+
        scale_fill_manual(values = c("purple4","#21908CFF","orange"))+
        labs(title = glue("{region}: {gene}"),
             x = "Corticosterone (GD14) (ng/ml)",
             y = "Normalized gene expression",
             color = "Status",
             fill = "Status")+
        theme_bw()+
        theme(legend.position = "bottom")  
    )
  }
}



limma_cort_exclude <- function(gene = "Avp", region = "vHYP", exclude = "Subdominant"){
  grcm38 %>% 
    filter(symbol == gene) %>% 
    .$ensgene %>% 
    unique() -> goi 
  
  if (goi %in% rownames(y_list[[region]]$E) == F)
  {temp_p <- ggplot()+ggtitle(glue("No detectable gene expression of {gene} in {region}"))}
  else {
    
    suppressMessages(
      y_list[[region]]$E[goi,] %>%
        as.data.frame(col.names = c("Expression")) %>% 
        rename(Expression = ".") %>%  
        rownames_to_column("sampleID") %>% 
        left_join(y_list[[region]]$target %>% 
                    rownames_to_column("sampleID") %>% 
                    select(sampleID, group)) %>% 
        left_join(sample_id) %>% 
        left_join(sample_data) %>%
        filter(group != exclude) %>% 
        ggplot(aes(cort_post,Expression, color = group, fill = group)) + 
        geom_smooth(method = 'lm', alpha = 0.1)+
        geom_point(shape = 21, alpha = 0.3, size = 2.5)+
        scale_color_manual(values = c("purple4","orange"))+
        scale_fill_manual(values = c("purple4","orange"))+
        labs(title = glue("{region}: {gene}"),
             x = "Corticosterone (GD14) (ng/ml)",
             y = "Normalized gene expression",
             color = "Status",
             fill = "Status")+
        theme_bw()+
        theme(legend.position = "bottom")  
    )
  }
}



plot_any_gene_cort <- function (gene = "Avp", region = "vHYP", exclude = F){
  grcm38 %>% 
    filter(symbol == gene) %>% 
    .$description -> desc
  strrep("-",str_length(glue("{gene}: {desc}"))) -> mycat # I am making a lot of effort for trivial aesthetic detail...
  # gene info 
  cat(glue("{gene}: {desc} \n"+glue("{mycat}\n"+glue("in {region}"))))
  
  if(exclude == F){
    limma_cort_plot(gene, region) -> temp_p  
  }
  else{limma_cort_exclude(gene = gene, region = region, exclude = exclude) -> temp_p  }
  print(temp_p)
}


plot_any_gene_cort_multiregion <- function(gene = "Avp",
                                      regions_of_interest = regions,
                                      exclude = F,
                                      ncol = 2){
  grcm38 %>%
    filter(symbol == gene) %>%
    .$description -> desc
  strrep("-",str_length(glue("{gene}: {desc}"))) -> mycat # I am making a lot of effort for trivial aesthetic detail...

  # gene info
  cat(glue("{gene}: {desc} \n"+glue("{mycat}\n"+glue(" in "))))
  temp_p_list <- list()


  if(exclude == F ){
      for(i in 1:length(regions_of_interest)){
        cat(glue("\n {regions_of_interest[i]} \n "))
        temp_p_list[[i]] <- (limma_cort_plot(gene,regions_of_interest[i])+
          theme(legend.position = 'none'))}
    }else{
      for(i in 1:length(regions_of_interest)){
        cat(glue("\n {regions_of_interest[i]} \n "))
        temp_p_list[[i]] <- limma_cort_exclude(gene,regions_of_interest[i])+
          theme(legend.position = "none")}
    }

  for_legend <- limma_cort_plot(gene,regions_of_interest[i])
  cowplot::get_legend(for_legend) -> my_legend

  temp_p_list[1:length(regions_of_interest)]/my_legend

}


# Crisscross ===========================================================
y_list <-list()
for(i in 1: length(regions)){
  my_region <- regions[i]
  y_list[[my_region]] <- readRDS(glue("results_RNAseqRDS/limma_vdl_{my_region}"))
}



crisscross <- function(two_genes = c("Pomc","Pomc"),
                       two_regions = c("mPOA","vHYP"),
                       by_status = FALSE){
  
  grcm38 %>% 
    filter(symbol == two_genes[1]) %>% 
    .$ensgene %>% 
    unique() -> goi1 
  
  grcm38 %>% 
    filter(symbol == two_genes[2]) %>% 
    .$ensgene %>% 
    unique() -> goi2
  
  
  suppressMessages(y_list[[two_regions[1]]]$E[goi1,] %>%
                     as.data.frame %>% 
                     rename(Expression1 = ".") %>%  #error here because column names "." doesnt exist.
                     rownames_to_column("sampleID") %>% 
                     left_join(sample_data %>% select(status, sampleID,subjectID)) %>% 
                     select(-sampleID) -> plotdf1)
  
  suppressMessages(y_list[[two_regions[2]]]$E[goi2,] %>%
                     as.data.frame %>% 
                     rename(Expression2 = ".") %>%  
                     rownames_to_column("sampleID") %>% 
                     left_join(sample_data %>% select(status, sampleID,subjectID)) %>% 
                     select(-sampleID) -> plotdf2)
  
  suppressMessages(full_join(plotdf1,plotdf2) -> plotdf)
  
  if(by_status == FALSE){
    suppressMessages(plotdf %>% 
                       ggplot(aes(Expression1,Expression2))+
                       geom_smooth(method = 'lm', alpha = 0.1)+
                       geom_point(shape = 21, alpha = 0.3, size = 2.5)+
                       labs(x = glue("{two_genes[1]} in {two_regions[1]}"),
                            y = glue("{two_genes[2]} in {two_regions[2]}"))+
                       theme_bw())
    
  }else{
    
    suppressMessages(plotdf %>% 
                       ggplot(aes(Expression1,Expression2,color = status, fill = status))+
                       geom_smooth(method = 'lm', alpha = 0.1)+
                       geom_point(shape = 21, alpha = 0.3, size = 2.5)+
                       scale_color_manual(values = c("purple4","#21908CFF","orange"))+
                       scale_fill_manual(values = c("purple4","#21908CFF","orange"))+
                       labs(x = glue("{two_genes[1]} in {two_regions[1]}"),
                            y = glue("{two_genes[2]} in {two_regions[2]}"),
                            color = "Status",
                            fill = "Status")+
                       theme_bw()+
                       theme(legend.position = "none")  )
  } 
}







crisscross_quick_lm <- function(two_genes = c("Avp","Oxt"),
                       two_regions = c("vHYP","vHYP")){
  
  grcm38 %>% 
    filter(symbol == two_genes[1]) %>% 
    .$ensgene %>% 
    unique() -> goi1 
  
  grcm38 %>% 
    filter(symbol == two_genes[2]) %>% 
    .$ensgene %>% 
    unique() -> goi2
  
  
  suppressMessages(y_list[[two_regions[1]]]$E[goi1,] %>%
                     as.data.frame %>% 
                     rename(Expression1 = ".") %>%  
                     rownames_to_column("sampleID") %>% 
                     left_join(sample_data %>% select(status, sampleID,subjectID)) %>% 
                     select(-sampleID) -> plotdf1)
  
  suppressMessages(y_list[[two_regions[2]]]$E[goi2,] %>%
                     as.data.frame %>% 
                     rename(Expression2 = ".") %>%  
                     rownames_to_column("sampleID") %>% 
                     left_join(sample_data %>% select(status, sampleID,subjectID)) %>% 
                     select(-sampleID) -> plotdf2)
  
  suppressMessages(full_join(plotdf1,plotdf2) -> plotdf)

  print(summary(lm(plotdf$Expression2 ~ plotdf$Expression1))  )
  print(cor.test(plotdf$Expression2, plotdf$Expression1,method =  "spearman"))
}



# despotism, alpha only ========================================


alpha_desp_plot_any_gene <- function(gene = "Cartpt", region = "vHYP",
                                     base_size = 11){
  
  
  limma_alpha_desp <- readRDS(glue("results_RNAseqRDS/limma_alpha_desp_{region}.RDS")) %>% 
    filter(symbol == gene)
  v.dl_alpha_desp <- readRDS(glue("results_RNAseqRDS/limma_vdl_alpha_desp_{region}"))
  
  limma_alpha_desp$logFC -> logfc
  limma_alpha_desp$P.Value -> pval
  
  grcm38 %>% 
    filter(symbol == gene) %>% 
    .$ensgene %>% 
    unique() -> goi
  
  v.dl_alpha_desp$E[goi,] %>% 
    as.data.frame() %>% 
    rename(Expression  = ".") %>% 
    rownames_to_column("sampleID") %>% 
    left_join(sample_data %>% select(sampleID, subjectID)) %>% 
    left_join(sample_data) %>% 
    ggplot(aes(despotism, Expression))+
    geom_point(size = 2.5, shape = 21)+
    geom_smooth(method = 'lm', alpha = 0.1)+
    theme_bw(base_size = base_size)+
    
    labs(x = "Despotism",
         y = "Normalized gene expression",
      title = glue("{region}: {gene}, logFC= {round(logfc,2)}, eFDR= {round(pval,3)}")) -> pp
  
  print(pp)
  print(limma_alpha_desp)  
  
}


# by bwchangerate  ===============================================================


bwchangerate_plot_any_gene <- function(gene = "Cartpt", region = "vHYP",
                                     base_size = 11){
  
  
  limma_bwchangerate <- readRDS(glue("results_RNAseqRDS/limma_bwchangerate_{region}.RDS")) %>% 
    filter(symbol == gene)
  v.dl_bwchangerate <- readRDS(glue("results_RNAseqRDS/limma_vdl_bwchangerate_{region}"))
  
  limma_bwchangerate$logFC -> logfc
  limma_bwchangerate$P.Value -> pval
  
  grcm38 %>% 
    filter(symbol == gene) %>% 
    .$ensgene %>% 
    unique() -> goi
  
  v.dl_bwchangerate$E[goi,] %>% 
    as.data.frame() %>% 
    rename(Expression  = ".") %>% 
    rownames_to_column("sampleID") %>% 
    left_join(sample_data %>% select(sampleID, subjectID)) %>% 
    left_join(sample_data) %>% 
    ggplot(aes(despotism, Expression))+
    geom_point(size = 2.5, shape = 21)+
    geom_smooth(method = 'lm', alpha = 0.1)+
    theme_bw(base_size = base_size)+
    
    labs(x = "Body weight change (%)",
         y = "Normalized gene expression",
         title = glue("{region}: {gene}, logFC= {round(logfc,2)}, eFDR= {round(pval,3)}")) -> pp
  
  print(pp)
  print(limma_bwchangerate)  
  
}
