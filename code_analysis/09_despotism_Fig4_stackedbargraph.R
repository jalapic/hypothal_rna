
readRDS("results_RDS/sample_data.RDS") -> df

df %>% 
  dplyr::select(cohort, subjectID, Win, Loss,glicko_rank) -> dfx

df %>% 
  group_by(cohort) %>% 
  summarize(totalwins = sum(Win)) -> total

df %>% 
  filter(glicko_rank == 1) %>% 
  dplyr::select(cohort, Win) -> alpha

total
alpha %>% 
  left_join(total) %>% 
  mutate(despot = (Win/totalwins)) %>% 
  mutate(cohort = as.numeric(factor(cohort))+104) %>% 
  filter(cohort !=111) %>% 
  arrange(desc(despot)) %>% 
  mutate(alphaname = glue("Alpha {row.names(.)}")) -> plotdf

summary(plotdf)

i = 3
x = 3.5
png(filename = glue("results_figures/despotism_stacked_barx.png"),
    width = 14, height = 16, units = "cm", res = 600)

plotdf %>% 
  mutate(rest = totalwins - Win) %>% 
  dplyr::select(-totalwins) %>% 
  relocate(Win, rest) %>% 
  gather(key, value, 1:2) %>% 
  ggplot(aes(x = reorder(alphaname,desc(despot)), y = value, fill = key, color = key))+
  geom_col(position = "fill")+
  labs(x = "",
       y = "Win proportion",
       color = "",
       fill = "",
       title = "Despotism across 11 alpha males")+
  scale_color_manual(values = c("#deebf7","#3182bd"))+
  scale_fill_manual(values = c("#deebf7","#3182bd"))+
  theme_classic(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1.5),
        legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())+
  annotate("text", x = 1, y = plotdf$despot[1]-0.05, label = round(plotdf$despot[1],i),size = x)+
  annotate("text", x = 2, y = plotdf$despot[2]-0.05, label = round(plotdf$despot[2],i),size = x)+
  annotate("text", x = 3, y = plotdf$despot[3]-0.05, label = round(plotdf$despot[3],i),size = x)+
  annotate("text", x = 4, y = plotdf$despot[4]-0.05, label = round(plotdf$despot[4],i),size = x)+
  annotate("text", x = 5, y = plotdf$despot[5]-0.05, label = round(plotdf$despot[5],i),size = x)+
  annotate("text", x = 6, y = plotdf$despot[6]-0.05, label = round(plotdf$despot[6],i),size = x)+
  annotate("text", x = 7, y = plotdf$despot[7]-0.05, label = round(plotdf$despot[7],i),size = x)+
  annotate("text", x = 8, y = plotdf$despot[8]-0.05, label = round(plotdf$despot[8],i),size = x)+
  annotate("text", x = 9, y = plotdf$despot[9]-0.05, label = round(plotdf$despot[9],i),size = x)+
  annotate("text", x = 10, y = plotdf$despot[10]-0.05, label = round(plotdf$despot[10],i),size = x)+
  annotate("text", x = 11, y = plotdf$despot[11]-0.05, label = round(plotdf$despot[11],i),size = x)



dev.off()
