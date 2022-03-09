
WGCNA_part_1 <- function(my_region){
  

  datExpr<- readRDS(glue("results_WGCNA/datExpr_{my_region}_limma_vdl_noIEGS.RDS"))
  
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, )
  # Plot the results:
  print(sft)
  sft$powerEstimate -> mynumber
  
  png(filename = glue("results_figures/WGCNA_{my_region}_softpower{mynumber}_limma_vdl_noIEGS.png"),
      width = 20, height = 12, units = "cm", res = 600)
  
  # sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = glue("{my_region}: Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  
  abline(h=0.85,col="red")  #changed to 0.8 (originally 0.9 in the tutorial)
  abline(h=0.90,col="blue")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = glue("{my_region}:Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  
  invisible(dev.off())
  
}


WGCNA_part_1("mPOA")  
WGCNA_part_1("vHYP") 

my_region = "vHYP"
my_power =  4
my_TOMType = "signed"
my_networkType = "signed hybrid"

WGCNA_get_net <- function(my_region = "mPOA",
                          my_power =  4, 
                          my_TOMType = "unsigned", 
                          my_networkType = "unsigned"){
  
  x <- readRDS(glue("results_WGCNA/datExpr_{my_region}_limma_vdl_noIEGS.RDS"))  
  net = blockwiseModules(x, 
                         power = my_power,
                         TOMType = my_TOMType, 
                         networkType = my_networkType,
                         minModuleSize = 50,
                         reassignThreshold = 0, 
                         mergeCutHeight = 0.25,
                         numericLabels = TRUE, 
                         pamRespectsDendro = FALSE,
                         saveTOMs = FALSE,
                         verbose = 3)
  
  
  saveRDS(net, glue("results_WGCNA/net_limma_vdl_{my_networkType}_{my_region}_{my_power}_noIEGS.RDS"))
  
}


WGCNA_get_net("mPOA", 4, "signed", "signed hybrid")

WGCNA_get_net("vHYP", 4, "signed", "signed hybrid") 


# ========================================================================================

files <- list.files(path = "results_WGCNA", pattern = "net_limma_vdl_*")
files <- files[-c(1:10)] #make sure im working with files i want
files

for (i in 1:length(files)){
  filename <- files[i]  #changed 1 to i
  
  net <- readRDS(glue('results_WGCNA/{filename}')) #its reading in vHYP here first even though its 2nd in list
  
  # open a graphics window
  sizeGrWindow(12, 9)
  # Convert labels to colors for plotting
  mergedColors = labels2colors(net$colors)
  
  # Plot the dendrogram and the module colors underneath
  
  dev.off() # make sure you do this before AND after 
  png(file = glue("results_figures/cluster_dendo_{str_sub(filename,5,-5)}.png"),
      width=600, height=350)
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05, 
                      main = str_sub(filename, 5, -5))
  
  dev.off()
}
  
# save WGCNA =============================================================================
mynrow = 5
myheight = 28

my_region = "mPOA" 

files <- list.files(path = "results_WGCNA", pattern = "net_limma_vdl_*")

files

 for (i in 1:length(files)){ 
my_region <- str_sub(gsub("net_limma_vdl_signed hybrid_","",files[i]),24,-14) #changed from -7 to -20 temporarily
my_power <- as.numeric(str_extract(files[i], "\\-*\\d+\\.*\\d*"))

net <- readRDS(glue('results_WGCNA/{files[i]}')) 
datExpr <- readRDS(glue("results_WGCNA/datExpr_{my_region}_limma_vdl_noIEGS.RDS"))  



sample_data %>% 
  filter(sampleID %in% rownames(datExpr)) %>% 
  mutate(statusx = ifelse(status == "Alpha", 2,ifelse(status == "Subdominant", 1,0))) %>%
  left_join(sample_data ) %>% 
  column_to_rownames("sampleID") %>% 
  dplyr::select(ds, status, statusx, cort_post, GD14_bw, bw_change_rate)-> datTraits 


datTraits 


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

moduleNumber = length(unique(moduleColors))

table(moduleColors)

save(MEs, moduleLabels, moduleColors, geneTree, datExpr,
     file = glue("results_WGCNA/{my_region}-networkConstruction-auto_power{my_power}_noIEGS.RData"))



# Module trait relationships ================================================================

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes   #make sure WGCNA package is loaded. kept getting error here
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits %>% select(-status), use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


png(file = glue("results_figures/module_trait_{my_region}_power{my_power}_noIEGS.png"),
    width = 12, height = 27, units = "cm", res = 600)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits %>% select(-status)),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


dev.off()


# look at regression ==================================================================

MEs %>%
  rownames_to_column("sampleID") %>%
  left_join(datTraits %>%
              rownames_to_column("sampleID")) -> ME_df



ME_df %>%
  gather(modules,value,2:(1+moduleNumber) ) %>%
  ggplot(aes(cort_post, value))+
  geom_point(alpha = 0.5)+
  geom_smooth(method = "lm", alpha = 0.2)+
  facet_wrap(~modules, scales = "free_y", nrow = mynrow)+
  theme_bw()+
  labs(x = "Corticosterone (ng/ml)",
       y = "Module eigengene",
       title = glue("{my_region}: Module eigengene ~ plasma corticosterone")) -> p2

png(filename = glue("results_figures/ME_cort_{my_region}_power{my_power}_noIEGS.png"),
    width = 30, height = myheight, units = "cm", res = 600)

print(p2)
invisible(dev.off())


ME_df %>%
  gather(modules,value,2:(1+moduleNumber) ) %>%
  # filter(modules %in% c("MEmagenta","MEcyan","MEgreen")) %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~modules, scales = "free_y", nrow = mynrow)+
  theme_bw()+
  theme(legend.position = "none") +
  scale_x_discrete(labels=c("Alpha" = "Alpha", "Subdominant" = "Subdom",
                            "Subordinate" = "Sub"))+
  labs(x = "Social status",
       y = "Module eigengene",
       title = glue("{my_region}: Module eigengene across social status")) -> p3

png(filename = glue("results_figures/ME_status_{my_region}_power{my_power}_noIEGS.png"),
    width = 35, height = myheight, units = "cm", res = 600)

print(p3)
invisible(dev.off())


# Module gene numbers ==============

moduleColors %>%
  table() %>%
  as.data.frame() %>% arrange(Freq)  -> modnum
colnames(modnum) <- c("module","count")
str(modnum)

modnum$module <- factor(modnum$module,levels=unique(modnum$module[order(modnum$count)]))
mycol <- modnum$module %>% as.character() # for some reason this has to be character again to be ordered properly in the figure...!!


modnum %>%
  filter(module != "grey") %>%
  ggplot(aes(y = module, x = count, fill = module))+
  geom_bar(stat = 'identity', color = NA)+
  geom_text(aes(label = count),hjust = -0.05, size = 1.4)+
  xlim(0, max(modnum$count)+75)+
  scale_fill_manual(values = mycol)+
  theme_minimal(base_size = 6)+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(margin = margin(r=-6)),
        panel.grid = element_blank())+
  labs(x = "",
       y = "",
       title = glue("{my_region}: Module size")) -> temp_p
temp_p

png(filename = glue("results_figures/modulesize_{my_region}_power{my_power}_noIEGS.png"),
    width = 6, height = 6, units = "cm", res = 600)
print(temp_p)
invisible(dev.off())



# Calculate kIN ====================================

ADJ=abs(cor(datExpr,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ, moduleColors) 

Alldegrees1 %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38 %>% select(ensgene,symbol)) -> kIN_df

xx <- cbind(ensgene = colnames(datExpr), module = moduleColors) %>% 
  as_tibble()


kIN_df %>% 
  left_join(xx) -> kIN_df

saveRDS(kIN_df,glue("results_RNAseqRDS/kIN_dataframe_{my_region}_noIEGS.RDS"))


# Make MM and GS dataframe =========================
# Define variable David's score containing the David's score column of datTrait
ds = as.data.frame(datTraits$ds);
names(ds) = "Davids_score"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, ds, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(ds), sep="");
names(GSPvalue) = paste("p.GS.", names(ds), sep="");


geneModuleMembership %>% 
  rownames_to_column("ensgene") %>% 
  left_join(geneTraitSignificance %>% rownames_to_column("ensgene")) -> gene_MM_TS_all

gene_MM_TS_all$module <- moduleColors
length(moduleColors)

get_mm <- function(x){
  x$moduleMembership <- x[colnames(x) == paste("MM",x$module,sep = "")] %>%
  unlist %>%
   as.numeric
  xx <- x %>%
    dplyr::select(ensgene,module,moduleMembership,GS.Davids_score)
  return(xx)
}



wgcna_whole <- get_mm(gene_MM_TS_all[1,])

for(i in 2:nrow(gene_MM_TS_all)){
  wgcna_whole <- rbind(wgcna_whole,get_mm(gene_MM_TS_all[i,]))
}

cort_post = as.data.frame(datTraits$cort_post);
names(cort_post) = "cort"
geneTraitSignificance = as.data.frame(cor(datExpr, cort_post, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(cort_post), sep="");
names(GSPvalue) = paste("p.GS.", names(cort_post), sep="");

wgcna_whole %>%
  left_join(geneTraitSignificance %>% rownames_to_column("ensgene")) %>% 
  left_join(grcm38 %>% select(ensgene, symbol)) %>% 
  select(-ensgene) %>% 
  as_tibble() %>% 
  relocate(symbol)-> wgcna_all 

saveRDS(wgcna_all, glue("results_RNAseqRDS/WGCNA_MM_GS_all_{my_region}_noIEGS.RDS"))
}




# TOM plot ===========================================================================
my_power = 4 #(finalized 061721, both for mPOA and vHYP)

files <- list.files(path = "results_WGCNA", pattern = "net_limma_vdl_*")

net <- readRDS(glue('results_WGCNA/{files[1]}')) #net_limma_vdl_signed hybrid_ files have 15766 (mPOA) and 15258 (vHYP)
#net <- readRDS(glue('results_WGCNA/{files[2]}')) 

my_region = "mPOA"
#my_region = "vHYP"

lnames = load(glue("results_WGCNA/{my_region}-networkConstruction-auto_power{my_power}_noIEGS.RData"))
#The variable lnames contains the names of loaded variables.
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


moduleLabels = net$colors  # make sure net is correct region 
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

moduleNumber = length(unique(moduleColors))

restGenes= (moduleColors != "grey") # MFD: why get rid of grey module?
#restGenes= moduleColors # MFD: why get rid of grey module?


dissTOM = 1-TOMsimilarityFromExpr(datExpr[,restGenes],
                                  power = my_power,
                                  networkType = "signed hybrid")
plotTOM = dissTOM^9; 
# look at the network w/out "grey" modules
colnames(dissTOM) = rownames(dissTOM) = colnames(datExpr[restGenes])
hier1=flashClust(as.dist(dissTOM), method = "average" )  


plotDendroAndColors(hier1, colors = data.frame(moduleColors[restGenes]),
                    c("Module Tree Cut"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
                    guideHang = 0.05, main = "Gene dendrogram and module colors")

# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
diag(dissTOM) = NA


png(glue("results_figures/TOMplot_{my_region}_{my_power}_noIEGS.png"),   width = 12, height = 12, units = "cm", res = 600)

TOMplot(dissim = 1-dissTOM^7, 
        hier1, 
        as.character(moduleColors[restGenes]), 
        min = "Network heatmap plot, removed unassigned genes")
dev.off()



# get MM ####

mpoa_mm<- as.data.frame(readRDS("results_RNAseqRDS/WGCNA_MM_GS_all_mPOA_noIEGs.RDS"))
vhyp_mm<- as.data.frame(readRDS("results_RNAseqRDS/WGCNA_MM_GS_all_vHYP_noIEGs.RDS"))

mpoa_mm <- mpoa_mm[order(mpoa_mm$module,mpoa_mm$moduleMembership), ]
vhyp_mm <- vhyp_mm[order(vhyp_mm$module,vhyp_mm$moduleMembership), ]

# mpoa
mpoa_mm %>% 
  group_by(module) %>% 
  arrange(module, desc(abs(moduleMembership))) -> mpoa_mm_top

mpoa_mm_top %>% filter(module == "magenta") %>% slice(1:20) -> mpoa_mm_magenta
as.vector(mpoa_mm_magenta$symbol) -> mpoa_mm_magenta

mpoa_mm_top %>% filter(module == "cyan") %>% slice(1:20) -> mpoa_mm_cyan
as.vector(mpoa_mm_cyan$symbol) -> mpoa_mm_cyan

mpoa_mm_top %>% filter(module == "black") %>% slice(1:20) -> mpoa_mm_black
as.vector(mpoa_mm_black$symbol) -> mpoa_mm_black

mpoa_mm_top %>% filter(module == "lightgreen") %>% slice(1:20) -> mpoa_mm_lightgreen
as.vector(mpoa_mm_lightgreen$symbol) -> mpoa_mm_lightgreen

mpoa_mm_top %>% filter(module == "darkgrey") %>% slice(1:20) -> mpoa_mm_darkgrey
as.vector(mpoa_mm_darkgrey$symbol) -> mpoa_mm_darkgrey

mpoa_mm_top %>% filter(module == "grey") %>% slice(1:20) -> mpoa_mm_grey
as.vector(mpoa_mm_grey$symbol) -> mpoa_mm_grey

mpoa_mm_top %>% filter(module == "brown") %>% slice(1:20) -> mpoa_mm_brown
as.vector(mpoa_mm_brown$symbol) -> mpoa_mm_brown

mpoa_mm_top %>% filter(module == "tan") %>% slice(1:20) -> mpoa_mm_tan
as.vector(mpoa_mm_tan$symbol) -> mpoa_mm_tan

mpoa_mm_top %>% filter(module == "lightyellow") %>% slice(1:20) -> mpoa_mm_lightyellow
as.vector(mpoa_mm_lightyellow$symbol) -> mpoa_mm_lightyellow

mpoa_mm_top %>% filter(module == "darkred") %>% slice(1:20) -> mpoa_mm_darkred
as.vector(mpoa_mm_darkred$symbol) -> mpoa_mm_darkred

mpoa_mm_top %>% filter(module == "blue") %>% slice(1:20) -> mpoa_mm_blue
as.vector(mpoa_mm_blue$symbol) -> mpoa_mm_blue

mpoa_mm_top %>% filter(module == "purple") %>% slice(1:20) -> mpoa_mm_purple
as.vector(mpoa_mm_purple$symbol) -> mpoa_mm_purple

mpoa_mm_top %>% filter(module == "red") %>% slice(1:20) -> mpoa_mm_red
as.vector(mpoa_mm_red$symbol) -> mpoa_mm_red

mpoa_mm_top %>% filter(module == "darkturquoise") %>% slice(1:20) -> mpoa_mm_darkturquoise
as.vector(mpoa_mm_darkturquoise$symbol) -> mpoa_mm_darkturquoise

mpoa_mm_top %>% filter(module == "turquoise") %>% slice(1:20) -> mpoa_mm_turquoise
as.vector(mpoa_mm_turquoise$symbol) -> mpoa_mm_turquoise

mpoa_mm_top %>% filter(module == "midnightblue") %>% slice(1:20) -> mpoa_mm_midnightblue
as.vector(mpoa_mm_midnightblue$symbol) -> mpoa_mm_midnightblue

mpoa_mm_top %>% filter(module == "royalblue") %>% slice(1:20) -> mpoa_mm_royalblue
as.vector(mpoa_mm_royalblue$symbol) -> mpoa_mm_royalblue

mpoa_mm_top %>% filter(module == "grey60") %>% slice(1:20) -> mpoa_mm_grey60
as.vector(mpoa_mm_grey60$symbol) -> mpoa_mm_grey60

mpoa_mm_top %>% filter(module == "green") %>% slice(1:20) -> mpoa_mm_green
as.vector(mpoa_mm_green$symbol) -> mpoa_mm_green

mpoa_mm_top %>% filter(module == "greenyellow") %>% slice(1:20) -> mpoa_mm_greenyellow
as.vector(mpoa_mm_greenyellow$symbol) -> mpoa_mm_greenyellow

mpoa_mm_top %>% filter(module == "lightcyan") %>% slice(1:20) -> mpoa_mm_lightcyan
as.vector(mpoa_mm_lightcyan$symbol) -> mpoa_mm_lightcyan

mpoa_mm_top %>% filter(module == "salmon") %>% slice(1:20) -> mpoa_mm_salmon
as.vector(mpoa_mm_salmon$symbol) -> mpoa_mm_salmon

mpoa_mm_top %>% filter(module == "pink") %>% slice(1:20) -> mpoa_mm_pink
as.vector(mpoa_mm_pink$symbol) -> mpoa_mm_pink

mpoa_mm_top_list <- list(mpoa_mm_salmon, mpoa_mm_lightcyan, mpoa_mm_greenyellow, mpoa_mm_green, mpoa_mm_grey60,
     mpoa_mm_royalblue, mpoa_mm_midnightblue, mpoa_mm_turquoise, mpoa_mm_darkturquoise,
     mpoa_mm_red, mpoa_mm_purple, mpoa_mm_blue, mpoa_mm_darkred, mpoa_mm_lightyellow,
     mpoa_mm_tan, mpoa_mm_darkgrey,mpoa_mm_brown, mpoa_mm_grey, mpoa_mm_magenta,
     mpoa_mm_lightgreen, mpoa_mm_black, mpoa_mm_cyan, mpoa_mm_pink)

names(mpoa_mm_top_list) <- list("mpoa_mm_salmon", "mpoa_mm_lightcyan", "mpoa_mm_greenyellow", "mpoa_mm_green", "mpoa_mm_grey60",
                         "mpoa_mm_royalblue", "mpoa_mm_midnightblue", "mpoa_mm_turquoise", "mpoa_mm_darkturquoise",
                         "mpoa_mm_red", "mpoa_mm_purple", "mpoa_mm_blue", "mpoa_mm_darkred", "mpoa_mm_lightyellow",
                         "mpoa_mm_tan", "mpoa_mm_darkgrey","mpoa_mm_brown", "mpoa_mm_grey", "mpoa_mm_magenta",
                         "mpoa_mm_lightgreen", "mpoa_mm_black", "mpoa_mm_cyan", "mpoa_mm_pink")

# vhyp

vhyp_mm %>% 
  group_by(module) %>% 
  arrange(module, desc(abs(moduleMembership))) -> vhyp_mm_top

vhyp_mm_top %>% filter(module == "magenta") %>% slice(1:20) -> vhyp_mm_magenta
as.vector(vhyp_mm_magenta$symbol) -> vhyp_mm_magenta

vhyp_mm_top %>% filter(module == "cyan") %>% slice(1:20) -> vhyp_mm_cyan
as.vector(vhyp_mm_cyan$symbol) -> vhyp_mm_cyan

vhyp_mm_top %>% filter(module == "black") %>% slice(1:20) -> vhyp_mm_black
as.vector(vhyp_mm_black$symbol) -> vhyp_mm_black

vhyp_mm_top %>% filter(module == "lightgreen") %>% slice(1:20) -> vhyp_mm_lightgreen
as.vector(vhyp_mm_lightgreen$symbol) -> vhyp_mm_lightgreen

vhyp_mm_top %>% filter(module == "darkgrey") %>% slice(1:20) -> vhyp_mm_darkgrey
as.vector(vhyp_mm_darkgrey$symbol) -> vhyp_mm_darkgrey

vhyp_mm_top %>% filter(module == "grey") %>% slice(1:20) -> vhyp_mm_grey
as.vector(vhyp_mm_grey$symbol) -> vhyp_mm_grey

vhyp_mm_top %>% filter(module == "brown") %>% slice(1:20) -> vhyp_mm_brown
as.vector(vhyp_mm_brown$symbol) -> vhyp_mm_brown

vhyp_mm_top %>% filter(module == "tan") %>% slice(1:20) -> vhyp_mm_tan
as.vector(vhyp_mm_tan$symbol) -> vhyp_mm_tan

vhyp_mm_top %>% filter(module == "lightyellow") %>% slice(1:20) -> vhyp_mm_lightyellow
as.vector(vhyp_mm_lightyellow$symbol) -> vhyp_mm_lightyellow

vhyp_mm_top %>% filter(module == "darkred") %>% slice(1:20) -> vhyp_mm_darkred
as.vector(vhyp_mm_darkred$symbol) -> vhyp_mm_darkred

vhyp_mm_top %>% filter(module == "blue") %>% slice(1:20) -> vhyp_mm_blue
as.vector(vhyp_mm_blue$symbol) -> vhyp_mm_blue

vhyp_mm_top %>% filter(module == "purple") %>% slice(1:20) -> vhyp_mm_purple
as.vector(vhyp_mm_purple$symbol) -> vhyp_mm_purple

vhyp_mm_top %>% filter(module == "red") %>% slice(1:20) -> vhyp_mm_red
as.vector(vhyp_mm_red$symbol) -> vhyp_mm_red

vhyp_mm_top %>% filter(module == "darkturquoise") %>% slice(1:20) -> vhyp_mm_darkturquoise
as.vector(vhyp_mm_darkturquoise$symbol) -> vhyp_mm_darkturquoise

vhyp_mm_top %>% filter(module == "turquoise") %>% slice(1:20) -> vhyp_mm_turquoise
as.vector(vhyp_mm_turquoise$symbol) -> vhyp_mm_turquoise

vhyp_mm_top %>% filter(module == "midnightblue") %>% slice(1:20) -> vhyp_mm_midnightblue
as.vector(vhyp_mm_midnightblue$symbol) -> vhyp_mm_midnightblue

vhyp_mm_top %>% filter(module == "royalblue") %>% slice(1:20) -> vhyp_mm_royalblue
as.vector(vhyp_mm_royalblue$symbol) -> vhyp_mm_royalblue

vhyp_mm_top %>% filter(module == "grey60") %>% slice(1:20) -> vhyp_mm_grey60
as.vector(vhyp_mm_grey60$symbol) -> vhyp_mm_grey60

vhyp_mm_top %>% filter(module == "green") %>% slice(1:20) -> vhyp_mm_green
as.vector(vhyp_mm_green$symbol) -> vhyp_mm_green

vhyp_mm_top %>% filter(module == "greenyellow") %>% slice(1:20) -> vhyp_mm_greenyellow
as.vector(vhyp_mm_greenyellow$symbol) -> vhyp_mm_greenyellow

vhyp_mm_top %>% filter(module == "lightcyan") %>% slice(1:20) -> vhyp_mm_lightcyan
as.vector(vhyp_mm_lightcyan$symbol) -> vhyp_mm_lightcyan

vhyp_mm_top %>% filter(module == "salmon") %>% slice(1:20) -> vhyp_mm_salmon
as.vector(vhyp_mm_salmon$symbol) -> vhyp_mm_salmon

vhyp_mm_top %>% filter(module == "pink") %>% slice(1:20) -> vhyp_mm_pink
as.vector(vhyp_mm_pink$symbol) -> vhyp_mm_pink

vhyp_mm_top_list <- list(vhyp_mm_salmon, vhyp_mm_lightcyan, vhyp_mm_greenyellow, vhyp_mm_green, vhyp_mm_grey60,
                         vhyp_mm_royalblue, vhyp_mm_midnightblue, vhyp_mm_turquoise, vhyp_mm_darkturquoise,
                         vhyp_mm_red, vhyp_mm_purple, vhyp_mm_blue, vhyp_mm_darkred, vhyp_mm_lightyellow,
                         vhyp_mm_tan, vhyp_mm_darkgrey,vhyp_mm_brown, vhyp_mm_grey, vhyp_mm_magenta,
                         vhyp_mm_lightgreen, vhyp_mm_black, vhyp_mm_cyan, vhyp_mm_pink)

names(vhyp_mm_top_list) <- list("vhyp_mm_salmon", "vhyp_mm_lightcyan", "vhyp_mm_greenyellow", "vhyp_mm_green", "vhyp_mm_grey60",
                                "vhyp_mm_royalblue", "vhyp_mm_midnightblue", "vhyp_mm_turquoise", "vhyp_mm_darkturquoise",
                                "vhyp_mm_red", "vhyp_mm_purple", "vhyp_mm_blue", "vhyp_mm_darkred", "vhyp_mm_lightyellow",
                                "vhyp_mm_tan", "vhyp_mm_darkgrey","vhyp_mm_brown", "vhyp_mm_grey", "vhyp_mm_magenta",
                                "vhyp_mm_lightgreen", "vhyp_mm_black", "vhyp_mm_cyan", "vhyp_mm_pink")



