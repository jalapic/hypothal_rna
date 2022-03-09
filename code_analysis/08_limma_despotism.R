# How many random sampling
R = 5000
i = 4


for (i in length(regions)){
  
  regions[[i]] -> my_region
  
  dlNorm <- rawcount_list[[i]] %>% 
    column_to_rownames('ensgene')
  
  var_info <- sample_data %>% 
    filter(sampleID %in% colnames(dlNorm)) %>% 
    select(sampleID, subjectID, status, despotism) %>% 
    filter(status == "Alpha")
  
  var_info$despotism -> group.dl  #MFD: added this
  
  dlNorm %>% 
    select_if(colnames(.) %in% var_info$sampleID) -> dlNorm  
  
  colnames(dlNorm)
  
  dlNorm <- dlNorm[!is.na(rowSums(dlNorm)),]
  
  d = apply(dlNorm, 2, as.numeric)
  dim(d)
  
  d0= DGEList(d, group = group.dl)
  dim(d0)
  rownames(d0) <- rownames(dlNorm)
  d0 <- calcNormFactors(d0)
  
  cutoff <- 5
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  dge.dl <- d0[-drop,]
  dim(dge.dl)
  
  var_info$despotism -> desp.dl
  
  design.dl <- model.matrix(~ desp.dl)
  colnames(design.dl) -> mycolnames
  
  v.dl = voom(dge.dl, design.dl, plot = F)
  vfit.dl = lmFit(v.dl, design.dl)
  efit.dl = eBayes(vfit.dl)
  
  

  p.dl.limma = efit.dl[["p.value"]]



saveRDS(v.dl, glue("results_RNAseqRDS/limma_vdl_alpha_desp_{my_region}"))

 p.dl.rand = vector('list',length = R)

  for(g in 1 : R){
    print(paste("Starting on Permutation", g))

    # Randomize the traits

    desp.dl.rand = sample(desp.dl)

    # Model
    design.dl.rand = model.matrix(~desp.dl.rand)
    colnames(design.dl.rand) <- mycolnames

    # Calculate p-values based on randomized traits
    v.dl.rand = voom(dge.dl, design.dl.rand, plot = F)
    vfit.dl.rand = lmFit(v.dl.rand, design.dl.rand)

    efit.dl.rand = eBayes(vfit.dl.rand)

    p.dl.rand[[g]] = efit.dl.rand[["p.value"]]
    head(p.dl.rand[[g]])
  }

  q.dl = matrix(0, nrow = nrow(p.dl.limma), ncol = ncol(p.dl.limma))

  for(h in 1 : R){
    print(paste("Calculating Permutation", h))

    temp = p.dl.rand[[h]]

    for(c in 1 : 2){
      for(r in 1 : nrow(p.dl.limma)){
        if(temp[r, c] <= p.dl.limma[r, c]){
          q.dl[r, c] = q.dl[r, c] + 1
        }
      }
    }
  }

  q.dl = q.dl / R
  colnames(q.dl) <- mycolnames
  q.dl = as.data.frame(q.dl)
  row.names(q.dl) <- rownames(dge.dl)

  saveRDS(q.dl,glue("results_RNAseqRDS/limma_eFDR_{my_region}_cutoff5_R{R}_alpha_desp.RDS"))

  q.dl <- readRDS(glue("results_RNAseqRDS/limma_eFDR_{my_region}_cutoff5_R{R}_alpha_desp.RDS"))
  
  
  png(filename = glue("results_figures/eFDR_hist_{my_region}R{R}_alpha_desp.png"),
      width = 18, height = 17, units = "cm", res = 600)
  hist(q.dl-p.dl.limma, main = glue("{my_region}"))
  invisible(dev.off())

  efit.dl[["p.value"]] <- q.dl
  row.names(q.dl) <- NULL
  
  topTable(efit.dl, sort.by = "P", n = Inf) %>% 
    rownames_to_column('ensgene') %>% 
    left_join(grcm38) %>%
    filter(!is.na(symbol)) %>% 
    select(symbol,logFC,P.Value,adj.P.Val) -> limma_alpha_despot 
  
  saveRDS(limma_alpha_despot,glue("results_RNAseqRDS/limma_alpha_desp_{my_region}.RDS"))
  
}

# remove IEGs ####
limma_df<- readRDS(glue("results_RNAseqRDS/limma_alpha_desp_{my_region}.RDS")) 
  
ARG_list$rapid_PRGs %>% 
  rename(symbol = `Gene ID`) %>% 
  filter(!is.na(symbol)) %>% 
  left_join(grcm38 %>% select(ensgene, symbol)) %>% 
  .$symbol %>% 
  unique -> IEG_symbol

limma_df <- limma_df[ !limma_df$symbol %in% IEG_symbol, ]

saveRDS(limma_df,glue("results_RNAseqRDS/limma_alpha_desp_{my_region}_noIEGs.RDS"))


# GO analysis ===============================================================

my_ont = "BP"
my_showCategory = 10
my_logFC_threshold = 0.15

my_region = "mPOA"
my_region = "vHYP"


limma_df<- readRDS(glue("results_RNAseqRDS/limma_alpha_desp_{my_region}_noIEGs.RDS")) %>% 
 filter(abs(logFC) >= my_logFC_threshold) %>%
  filter(.,P.Value <0.05) 

limma_df %>% 
  mutate(updown = ifelse(logFC>0, "Up","Down")) %>% 
  group_by(updown) %>% 
  distinct %>% 
  count

gettopGO <- function(limma_df,my_showCategory = 5){
  go_df_up <- limma_df %>% 
    filter(logFC>0) %>% 
    left_join(grcm38 %>% dplyr::select(symbol, entrez)) %>% 
    mutate(entrez = as.character(entrez)) %>% 
    dplyr::select(entrez, logFC) %>% 
    filter(!is.na(entrez)) %>% 
    arrange(desc(logFC))
  
  
  go_df_down <- limma_df %>% 
    filter(logFC<0) %>% 
    left_join(grcm38 %>% dplyr::select(symbol, entrez)) %>% 
    mutate(entrez = as.character(entrez)) %>% 
    dplyr::select(entrez, logFC) %>% 
    filter(!is.na(entrez)) %>% 
    arrange(desc(logFC))
  
  ggo_up <- enrichGO(gene = go_df_up$entrez %>% unique(),
                     OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                     keyType = "ENTREZID",
                     ont = my_ont,
                     readable = T,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.1,
                     qvalueCutoff  = 0.50)
  
  ggo_down <- enrichGO(gene = go_df_down$entrez %>% unique() ,
                       OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                       keyType = "ENTREZID",
                       ont = my_ont,
                       readable = T,
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.1, # change to 0.2 for Liver-status only
                       qvalueCutoff  = 0.50) 
  
  fortify(
    ggo_up,
    showCategory = my_showCategory,
    by = "Count",
    split = NULL,
    includeAll = TRUE
  ) %>% 
    arrange(desc(GeneRatio)) %>% 
    mutate(direction = "Up") -> temp1
  
  
  fortify(
    ggo_down,
    showCategory = my_showCategory,
    by = "Count",
    split = NULL,
    includeAll = TRUE
  ) %>% 
    arrange(desc(GeneRatio)) %>% 
    mutate(direction = "Down") -> temp2
  return(rbind(temp1,temp2))
  
}


limma_df %>% 
  gettopGO -> GO_result


write.csv(GO_result, glue("results_tables/GO_results_despot_{my_region}_noIEGs.csv"), row.names = F)

  

# Genes with greatest LFC ==========================================

my_region = "mPOA"
#my_region = "vHYP"
limma_alpha_despot <- readRDS(glue("results_RNAseqRDS/limma_alpha_desp_{my_region}_noIEGs.RDS"))


## How many gene significantly high/low expressed in despotic alphas?
limma_alpha_despot %>%
  filter(P.Value <0.05) %>%
  filter(logFC >0) -> limma_alpha_despot_up
print(nrow(limma_alpha_despot_up))

limma_alpha_despot %>%
  filter(P.Value <0.05) %>%
  filter(logFC <0) -> limma_alpha_despot_down
print(nrow(limma_alpha_despot_down))


# filter genes from aggression gene set
grcm38 %>%
  filter(symbol %in% both) %>%
  arrange(symbol) %>%
  select(symbol, description) %>%
  as.data.frame


limma_alpha_despot %>%
  filter(symbol %in% mouse_agg_geneset) 



# Figure 4B =====
limma_alpha_despot %>%
  filter(P.Value <0.05) %>%
  arrange(logFC) %>%
  head(14) %>%
  .$symbol


limma_alpha_despot %>%
  filter(P.Value <0.05) %>%
  arrange(desc(logFC)) %>%
  head(20) %>%
  .$symbol


alpha_desp_plot_any_gene("Rpl21","mPOA")

alpha_desp_plot_any_gene("Rpl21","vHYP")
alpha_desp_plot_any_gene("Ache","vHYP")
alpha_desp_plot_any_gene("Cartpt","vHYP")


my_region = "mPOA"
gene = "Rpl21"


my_region = "vHYP"
gene = "Rpl21"
gene = "Ache"
gene = "Cartpt"


png(filename = glue("results_figures/alpha_desp_{gene}_{my_region}.png"),
    width = 7, height = 7, units = "cm", res = 600)

alpha_desp_plot_any_gene(gene, my_region, base_size = 7)

invisible(dev.off())

