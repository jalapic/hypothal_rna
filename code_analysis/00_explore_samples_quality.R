# total gene counts per sample ========================================

counts %>% 
  summarize_if(is.numeric,sum,na.rm = T) %>% 
  t() %>% 
  as.data.frame %>% 
  rename(genecounts = V1) %>% 
  rownames_to_column(var = 'sampleID') -> genecounts

genecounts %>% 
  ggplot(aes(genecounts)) +
  geom_histogram(bins = 40,alpha =0.5,color = 'grey') +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  theme_base() -> p_genecounts

print(p_genecounts)
#ggsave(p_genecounts,filename = 'results_figures/p_genecounts.png',width = 8, height = 5)



