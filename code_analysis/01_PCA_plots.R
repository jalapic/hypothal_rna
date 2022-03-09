i = 2
filter_counts = 50

library(DESeq2)

for (i in c(1:2)){
  

regions[[i]] -> my_region

df <- rawcount_list[[i]] %>% 
  dplyr::select(-ensgene)

coldata = sample_data %>% 
  filter(sampleID %in% colnames(df)) %>% 
  dplyr::select(sampleID, subjectID, status) #%>% 
  # filter(status == "Alpha") %>% # # despotism - alpha only 

# # despotism - alpha only 
# df <- df %>% 
#   select_if(colnames(.) %in% coldata$sampleID)

countData <- as.matrix(df) 
rownames(countData) <- rawcount_list[[i]]$ensgene 

countData <- countData[!is.na(rowSums(countData)),]
countData %>% is.na() %>% sum

dim(countData)

countData <- countData[rowSums(countData > filter_counts) > round((length(df))*0.9), ]

dim(countData)

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = coldata,
                              design = ~status)

# # despotism - alpha only 
# dds <- DESeqDataSetFromMatrix(countData = countData,
#                               colData = coldata,
#                               design = ~despotism)


dds <- DESeq(dds)

# saveRDS(dds, glue("results_RNAseqRDS/dds2_{my_region}.RDS"))


vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup = c("status"), returnData = T) -> d
attributes(d)
ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) +
  geom_point(size = 3, alpha = 0.5) +
  # geom_text(aes(label = subjectID), vjust = -0.5)+
  scale_color_manual(values = c("purple4", "#21908CFF","orange")) +
  # scale_fill_manual(values = c("purple4", "#21908CFF","orange")) +
  labs(x = paste0("PC1: ",
                  round(attributes(d)$percentVar[1] * 100),
                  "% variance"),
       y = paste0("PC2: ",
                  round(attributes(d)$percentVar[2] * 100),
                  "% variance"),
       color = "Social status",
       fill = "Social status") +
  theme(legend.position = "right",
    # legend.position = c(0.11,0.87),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))+
  ggtitle(glue('{my_region} PCA plot')) -> p



theme_set(theme_bw(base_size = 10))

png(filename = glue("results_figures/PCA_{my_region}_filtercounts{filter_counts}.png"),
    width = 10.5, height = 7.3, units = "cm", res = 600)
print(p)
invisible(dev.off())

}

