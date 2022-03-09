# feeding genes in vHYP x cort level #

#make expression dataframe for each gene ###
library(janitor)
region <- "vHYP"
vhyp_df <-  as.data.frame(y_list[[region]]$E)

goi <- grcm38[grcm38$symbol=="Npy",1]
goi <- as.vector(goi$ensgene)
npy <-  vhyp_df %>%
  filter(rownames(vhyp_df) %in% goi) %>%
        #as.data.frame(col.names = c("Expression"))%>% 
        #rename(Expression = ".")%>%  
        rownames_to_column("sampleID") %>%
    t() %>%
    row_to_names(row_number = 1)  %>% 
    as.data.frame() %>%
    rownames_to_column("sampleID")  %>% 
    left_join(y_list[[region]]$target %>% 
                    rownames_to_column("sampleID") %>% 
                    select(sampleID, group)) %>%
    left_join(sample_data %>% select(sampleID, cohort, cort_post) ) 
colnames(npy)[2] <- "Expression"
npy$Expression <- as.numeric(npy$Expression)


goi <- grcm38[grcm38$symbol=="Agrp",1]
goi <- as.vector(goi$ensgene)
agrp <-  vhyp_df %>%
  filter(rownames(vhyp_df) %in% goi) %>%
  #as.data.frame(col.names = c("Expression"))%>% 
  #rename(Expression = ".")%>%  
  rownames_to_column("sampleID") %>%
  t() %>%
  row_to_names(row_number = 1)  %>% 
  as.data.frame() %>%
  rownames_to_column("sampleID")  %>% 
  left_join(y_list[[region]]$target %>% 
              rownames_to_column("sampleID") %>% 
              select(sampleID, group)) %>%
  left_join(sample_data %>% select(sampleID, cohort, cort_post) ) 
colnames(agrp)[2] <- "Expression"
agrp$Expression <- as.numeric(agrp$Expression)

goi <- grcm38[grcm38$symbol=="Pomc",1]
goi <- as.vector(goi$ensgene)
pomc <-  vhyp_df %>%
  filter(rownames(vhyp_df) %in% goi) %>%
  #as.data.frame(col.names = c("Expression"))%>% 
  #rename(Expression = ".")%>%  
  rownames_to_column("sampleID") %>%
  t() %>%
  row_to_names(row_number = 1)  %>% 
  as.data.frame() %>%
  rownames_to_column("sampleID")  %>% 
  left_join(y_list[[region]]$target %>% 
              rownames_to_column("sampleID") %>% 
              select(sampleID, group)) %>%
  left_join(sample_data %>% select(sampleID, cohort, cort_post) ) 
colnames(pomc)[2] <- "Expression"
pomc$Expression <- as.numeric(pomc$Expression)


goi <- grcm38[grcm38$symbol=="Cartpt",1]
goi <- as.vector(goi$ensgene)
cartpt <-  vhyp_df %>%
  filter(rownames(vhyp_df) %in% goi) %>%
  #as.data.frame(col.names = c("Expression"))%>% 
  #rename(Expression = ".")%>%  
  rownames_to_column("sampleID") %>%
  t() %>%
  row_to_names(row_number = 1)  %>% 
  as.data.frame() %>%
  rownames_to_column("sampleID")  %>% 
  left_join(y_list[[region]]$target %>% 
              rownames_to_column("sampleID") %>% 
              select(sampleID, group)) %>%
  left_join(sample_data %>% select(sampleID, cohort, cort_post) ) 
colnames(cartpt)[2] <- "Expression"
cartpt$Expression <- as.numeric(cartpt$Expression)

goi <- grcm38[grcm38$symbol=="Gal",1]
goi <- as.vector(goi$ensgene)
gal <-  vhyp_df %>%
  filter(rownames(vhyp_df) %in% goi) %>%
  #as.data.frame(col.names = c("Expression"))%>% 
  #rename(Expression = ".")%>%  
  rownames_to_column("sampleID") %>%
  t() %>%
  row_to_names(row_number = 1)  %>% 
  as.data.frame() %>%
  rownames_to_column("sampleID")  %>% 
  left_join(y_list[[region]]$target %>% 
              rownames_to_column("sampleID") %>% 
              dplyr::select(sampleID, group)) %>%
  left_join(sample_data %>% select(sampleID, cohort, cort_post) ) 
colnames(gal)[2] <- "Expression"
gal$Expression <- as.numeric(gal$Expression)


saveRDS(npy, "results_RDS/npy.RDS")
saveRDS(gal, "results_RDS/gal.RDS")
saveRDS(agrp, "results_RDS/agrp.RDS")
saveRDS(cartpt, "results_RDS/cartpt.RDS")
saveRDS(pomc, "results_RDS/pomc.RDS")




# brms ####

library(brms)

genelist <- c("agrp","npy","gal","cartpt","pomc")


for(i in 1:5){
  
  my_gene <- genelist[i]
  
  
  generds <- readRDS(glue("results_RDS/{my_gene}.RDS")) #make sure your RDS file is in the noted folder.
  
  brm(Expression ~ cort_post + group + (1|cohort), 
      family= "student",
      chains = 3,
      file= glue("results_RNAseqRDS/brms_result_{my_gene}_status_fixed.RDS"),
      control = list(adapt_delta = 0.95),
      data = generds) -> bmod1
  
  
  dat<- conditional_effects(bmod1)$cort_post %>% 
    mutate(xvar=cort_post, 
           yvar=estimate__,
           lower=lower__,
           upper=upper__) %>% 
    select(xvar,yvar,lower,upper) %>%
    as.data.frame()
  
  dat
  
  geneplot <- ggplot()+
    geom_point(data=generds,
               aes(cort_post,Expression, fill = group, color = group),
               shape = 21, size = 1.5)+
    geom_point(data=generds,
               aes(cort_post,Expression),
               shape = 21, size = 1.5,  fill = "grey", alpha = 0.5)+
    geom_ribbon(data = dat,
                aes(xvar,ymax=upper,ymin=lower),
                alpha=0.1,fill="navy")+
    geom_line(data = dat,
              aes(xvar,yvar),
              size=1.5,alpha=0.5,color="navy")+
    labs(x="Corticosterone lelvel (GD14) (ng/ml)",
         y="Normalized gene expression",
         title = glue("vHYP: {str_to_title(my_gene)}"))+
    theme_bw()+
    # scale_color_brewer(palette = "Dark2")+
    # scale_fill_brewer(palette = "Dark2")+
    scale_color_manual(values = c("purple4","#21908CFF","orange"), name = "Status")+
    scale_fill_manual(values = c("purple4","#21908CFF","orange"), name = "Status")+
    theme(legend.position ="bottom")+
    #facet_wrap(vars(dat), scales = "free_y", nrow=2)+
    theme_bw(base_size = 6) #+
    #scale_fill_discrete(name= "Status")
  
  
  png(filename = glue("results_figures/{my_gene}_cort_expression_brms_new.png"),
      width = 8, height = 5, units = "cm", res = 600)
  print(geneplot)
  invisible(dev.off())
  
}


# look at CI
summary(brms_result_cartpt_new)
brms_result_cartpt_new %>%
  gather_draws(b_cort_post) %>%
  median_qi() -> brms_result_cartpt_new
brms_result_cartpt_new

