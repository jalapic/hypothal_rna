theme_set(theme_bw())

# scale_color_manual(values = c('#d8b365','#5ab4ac'))
 #   scale_fill_manual(values = c('#d8b365','#5ab4ac'))+
 #   scale_color_manual(values = c('#e9a3c9','#a1d76a'))+
 #   scale_color_manual(values = c('#ef8a62','#67a9cf'))+

 options(ggplot2.continuous.colour="viridis")
 options(ggplot2.continuous.fill="viridis")
 scale_colour_discrete <- scale_colour_viridis_d
 scale_fill_discrete <- function(...) {
   scale_fill_manual(..., values = viridis_qualitative_pal7)
 } 

PTO <- c("purple4","#21908CFF","orange") # purple, turquoise (or whater you call it), orange 
BTS <- c("#4169E1","#d8b365","#5ab4ac") #blue, teal, soil 
JCO2 <- c("#FC4E07","#d8b365","#5ab4ac") #similar to JCO https://rpkgs.datanovia.com/ggpubr/index.html
JCO <- c("#FC4E07","#E7B800","#00AFBB")

Won_preset <-theme(legend.position = c(0.8,0.2),
      legend.key.height = unit(0,"mm"),
      legend.title = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 3.3)))


