

library(brms)
library(tidybayes)


my_power = 4 #(finalized 061721, both for mPOA and VMH) 

my_region = "mPOA"
my_height = 7.56
my_width = 25.5

my_region = "vHYP"
my_height = 7.56
my_width = 20


lnames = load(glue("results_WGCNA/{my_region}-networkConstruction-auto_power{my_power}_noIEGs.RData")) #this is also in 12_WGCNA ~ line 435
#The variable lnames contains the names of loaded variables.
lnames
nGenes = ncol(datExpr)  
nSamples = nrow(datExpr)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes 
orderMEs(MEs0) %>% 
  rownames_to_column("sampleID") %>%
  left_join(sample_data %>% #left_join(sample_data) %>% 
              select(sampleID,status, cohort, cort_post, bw_change_rate) %>% 
              mutate(cohort = as.factor(cohort)) %>% 
              mutate_if(is.numeric,scale)) %>% 
  relocate(status, cohort, cort_post, bw_change_rate) %>% 
  mutate(statusx = factor(status, levels = c("Subdominant","Alpha","Subordinate"))) %>% 
  relocate(status, statusx) %>% 
  select(-sampleID)-> ME_df

colnames(ME_df)

lm_result_list <- list()

library(lme4)
library(lmerTest)


for(x in 1:length(MEs0)){
  k = x + 5
  ME_df[,c(1:5,k)] -> df
  md <- gsub("ME","",colnames(df)[6])
  colnames(df)[6] <- "module"
  # lmer(module ~ status +(1|cohort) , data = df) -> mod1
  # lmer(module ~ statusx +(1|cohort) , data = df) -> mod2
  # summary(mod1)
  # summary(mod2)
  brm(module ~ status +(1|cohort) , data = df) -> bmod1
  brm(module ~ statusx +(1|cohort) , data = df) -> bmod2
  # summary(bmod1)
  # plot(bmod1)
  bmod1 %>% 
    gather_draws(b_statusSubdominant, b_statusSubordinate) %>% 
    median_qi() -> sum1
  
  bmod2 %>% 
    gather_draws(b_statusxSubordinate) %>% 
    median_qi() -> sum2
  
  rbind(sum1, sum2) %>% 
    as.data.frame() %>% 
    cbind(key = c("Alpha-subdom","Alpha-Sub","Subdom-Sub")) %>% 
    mutate(module = md) -> lm_result_list[[x]] 
  

}

lm_result_list %>% 
  rbindlist -> lm_result_all

lm_result_all %>% 
  filter(module != "grey") -> lm_result_all
  

my_region <- "vHYP"
saveRDS(lm_result_all,glue("results_rnaseqRDS/lm_result_all_{my_region}_noIEGs.RDS"))


my_region <- "vHYP"
my_region <- "mPOA"


lm_result_all<-readRDS(glue("results_rnaseqRDS/lm_result_all_{my_region}_noIEGs.RDS"))


lm_result_all %>% 
  rename(Estimate = .value) %>% 
  mutate(my_alpha = ifelse((.lower)*(.upper) >0, 1, 0)) %>% 
  mutate(my_alpha2 = Estimate ) -> heatmap_df



moduleColors %>% 
  table() %>% 
  as.data.frame() %>% arrange(Freq)  -> modnum 
colnames(modnum) <- c("module","count")

heatmap_df %>% as_tibble() %>% 
  left_join(modnum ) -> heatmap_dfx



heatmap_df %>% 
  filter(my_alpha > 0.95) %>% 
  .$Estimate %>% 
  abs() %>% 
  max() -> my_limit
rev(levels(as.factor(heatmap_df$module))) -> xx
xx[xx != "grey"] -> y_limit


heatmap_dfx %>% 
  as_tibble() %>% 
  select(module, count) %>% 
  distinct() %>% 
  # mutate(module = factor(module)) %>% 
  arrange(desc(count)) %>% 
  filter(module!="grey") %>% 
  .$module -> my_module_level


key_level <- c("Alpha-subdom","Alpha-Sub","Subdom-Sub")


vv = length(my_module_level)-1
vlines = c(1:vv)+0.495
hlines = c(1:2)+0.51
ggplot(data = heatmap_df %>%  filter(my_alpha >0),
       aes(x = module, y =key, fill = Estimate))+
  geom_tile(color = "white",size =1 )+
  geom_vline(xintercept = vlines, color = "grey",size = 0.2)+
  geom_hline(yintercept = hlines, color = "grey",size = 0.2)+
  scale_x_discrete(limits = my_module_level, position = 'top')+
  scale_y_discrete(limits = rev(key_level))+
  scale_fill_distiller(palette = "PiYG",
                       limit = c(-my_limit,my_limit),
                       breaks = c(-0.1,0,0.1))+
  labs(x = "", y = "")+
  theme_bw(base_size = 15)+
  theme(legend.title = element_text(size =rel(1.2)),
        legend.text = element_text(size =rel(1.)),
        legend.position = "none",
        legend.key.size = unit(.9, 'cm'),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "grey"),
        # panel.background = element_blank(),
        # panel.grid.minor = element_line(color = 'grey50'),
        panel.grid = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm"),
        axis.text.y = element_text(hjust=0.95,vjust=0.2,size = rel(1.2)),
        axis.text.x = element_text(angle = 45,hjust=0.08,vjust=0.15,size = rel(1))) -> traitmodule2


traitmodule2



png(filename = glue("results_figures/traitmodule_{my_region}_noIEGs.png"),
    width = my_width, height = my_height, units = "cm", res = 600)

traitmodule2

invisible(dev.off())



