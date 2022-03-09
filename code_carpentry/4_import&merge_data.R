
# Import Tag-seq data
rawcount_list <- list()
temp = list.files(path = "data_raw",pattern="*_counts.csv")
temp0 <- paste0("data_clean/",temp)
rawcount_list = lapply(temp0, read_csv)
regions <- gsub("_counts.csv","",temp)

# Import count sample ID info 
sample_id <- read_csv('data_raw/sample_id.csv') %>% 
  rename(sampleID = sample_id)

# Import behavior data 
names(rawcount_list) <-regions

behavior <- read_csv('data_clean/hbdata.csv') %>% 
  mutate(cohort = as.numeric(as.factor(cohort))+104) %>% 
  mutate(subjectID = paste(cohort, mouseID, sep = '-'))  %>% 
  mutate(bw_change_rate = (GD14_bw-GD01_BW)/GD01_BW)


# Import data on time out of group housing prior to euthanasia 
time_out <- read_csv("data_raw/end_bodyweight.csv") %>% 
  mutate(subjectID = glue("{cohort}-{mouseID}")) %>% 
  select(subjectID,time_out) %>% 
  filter(subjectID %in% unique(sample_id$subjectID)) %>% 
  mutate(timeout_minute = lubridate::minute(time_out)) 

  # median time out of vivaria:
median(time_out$timeout_minute)
sd(time_out$timeout_minute)


# Import data on activity-regulated genes (ARGs) from Tyssowski et al. 2018 Neuron ========
# https://doi.org/10.1016/j.neuron.2018.04.001

#mmc3: Table S2. ARG-Seq and eRNA-Seq Data
readxl::read_excel("data_clean/Tyssowski_et_al_2018/mmc3.xlsx",
                   sheet = 2) %>% 
  rename(symbol = `gene name`,
         gene_class = `Gene Class`) %>% 
  select(symbol, gene_class) -> ARG_df

#mmc6: Table S5. Functional Annotation of ARGs by Class
readxl::read_excel("data_clean/Tyssowski_et_al_2018/mmc6.xlsx",
                   sheet = 1) 


"data_clean/Tyssowski_et_al_2018/mmc6.xlsx" -> xl_data
tab_names <- readxl::excel_sheets(path = xl_data)

ARG_list <- lapply(tab_names, function(x) 
  readxl::read_excel(xl_data, sheet = x)) 
names(ARG_list) <- gsub(" ","_",tab_names)


ARG_df2 <- ARG_list %>%
  map2_df(.,names(.), ~mutate(.x, gene_class = .y)) %>%
  rename(symbol = `Gene ID`) %>%
  select(symbol, gene_class)


# Import aggression gene set ======
my_geneset <- read_csv("data_clean/top40_aggressive_genes_ZhnagJames2018.csv")
mouse_agg_geneset <- convertHumanGeneList(my_geneset$human_symbol)

# Import neurondocrine gene set ======
source("code_functions/functions_convert_human_mouse_biomaRt.R")
my_geneset <- read_csv("data_clean/neuroendocrine_gene_set_Donkelaar2020.csv")
neuroendo_mouse_geneset <- rbind(
  cbind(location = "sex_chr", symbol = convertHumanGeneList(my_geneset$sex_chromosome)),
  cbind(location = "auto_chr", symbol = convertHumanGeneList(my_geneset$auto_chromosome))
) 

saveRDS(neuroendo_mouse_geneset, "results_RNAseqRDS/neuroendo_mouse_geneset.RDS")


# Merge behavior and sample ID ####
sample_data <- sample_id %>% 
  left_join(behavior) %>% 
  filter(!is.na(status))

saveRDS(sample_data, file = "results_RDS/sample_data.rds")

counts <- rawcount_list %>%
  reduce(full_join, by = 'ensgene')
