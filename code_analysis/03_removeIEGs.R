
ARG_list$rapid_PRGs %>% 
  rename(symbol = `Gene ID`) %>% 
  filter(!is.na(symbol)) %>% 
  left_join(grcm38 %>% select(ensgene, symbol)) %>% 
  .$symbol %>% 
  unique -> IEG_symbol

my_region <- "vHYP"
my_region <- "mPOA"

limma_genes <- readRDS(glue("results_RNAseqRDS/limma_{my_region}.RDS"))
limma_genes$alphasubdom <- limma_genes$alphasubdom[ ! limma_genes$alphasubdom$symbol %in% IEG_symbol, ]
limma_genes$alphasub <- limma_genes$alphasub[ ! limma_genes$alphasub$symbol %in% IEG_symbol, ]
limma_genes$subdomsub <- limma_genes$subdomsub[ ! limma_genes$subdomsub$symbol %in% IEG_symbol, ]

saveRDS(limma_genes,glue("results_RNAseqRDS/limma_{my_region}_noIEGs.RDS"))



E_limma <- readRDS(glue("results_RNAseqRDS/limma_vdl_{my_region}"))$E 

# MFD: sanity check: Both mine and Won's limma_vdl_{my_region} files are the same

dlNorm <- rawcount_list[[my_region]] %>% 
  column_to_rownames('ensgene')

datExpr <- t(E_limma)
datExpr <- as.data.frame(datExpr)   
datExpr <- datExpr[,both_wo_IEG] #need to get rid of IEGs
datExpr <- datExpr[,!(names(datExpr) %in% IEG_ensgene)] 

## check rows/cols
nrow(datExpr)
ncol(datExpr)
rownames(datExpr)
Samples = rownames(datExpr);
collectGarbage()

saveRDS(datExpr, glue("results_WGCNA/datExpr_{my_region}_limma_vdl_noIEGS.RDS"))
