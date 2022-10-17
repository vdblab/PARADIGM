library(ggplot2)
library(ggrepel)
library(yarr)
library(tidyverse)
library(ggpubr)

# download publicly available data from Maier et al., 2018 Nature (DOI: 10.1038/nature25979)
# download link: https://figshare.com/articles/dataset/Extensive_impact_of_non-antibiotic_drugs_on_human_gut_bacteria/4813882

combined_pv = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Aug2021/combined_pv.csv")
prestwick_atc = read_tsv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Aug2021/prestwick_atc.tsv")
species_overview = read_tsv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Aug2021/species_overview.tsv")
maier_data = combined_pv %>% 
  left_join(prestwick_atc %>% select(prestwick_ID, chemical_name), by = "prestwick_ID") %>% 
  left_join(species_overview %>% select(NT_code, species), by = "NT_code") %>% 
  filter(!duplicated(.)) %>% 
  mutate(chemical_name = as.character(chemical_name)) %>% 
  mutate(chemical_name2 = tolower(gsub(' [A-z ]*', '' , chemical_name)) ) %>% 
  mutate(chemical_name2 = case_when(chemical_name == "Baclofen (R,S)" ~ "baclofen",
                                    chemical_name == "Methylprednisolone, 6-alpha" ~ "methylprednisolone", 
                                    chemical_name == "Metoprolol-(+,-) (+)-tartrate salt" ~ "metoprolol", 
                                    chemical_name == "Ethinylestradiol" ~ "estradiol",
                                    chemical_name == "Estradiol Valerate" ~ "estradiol valerate",
                                    chemical_name == "Cyclosporin A" ~ "cyclosporine", 
                                    chemical_name == "Amlodipine" ~ "amlodipine_besylate",
                                    TRUE ~ chemical_name2) )
                  
# two antibiotics, piperacillin/tazobactam and sulfamethoxazole/trimethoprim, were administered in human as combination
# but were screened in vitro individually 
# we assume that drugs given in combination in the clinical setting will have synergistic effect in vitro

zosyn = maier_data %>% 
  filter(chemical_name2 == "piperacillin") %>% 
  select(NT_code, species, AUC, hit) %>% 
  left_join(maier_data %>% filter(chemical_name2 == "tazobactam") %>% select(NT_code, species, AUC, hit), by = "NT_code") %>% 
  mutate(AUC = AUC.x * AUC.y) %>% 
  mutate(hit.x = as.numeric(hit.x) ) %>% 
  mutate(hit.y = as.numeric(hit.y) ) %>% 
  mutate(hit = hit.x + hit.y) %>% 
  mutate(hit = case_when(hit >= 1 ~ T, T ~ F))
zosyn2 = maier_data %>% 
  filter(chemical_name2 == "piperacillin") %>% 
  mutate(chemical_name2 = "piperacillin_tazobactam") %>% 
  mutate(AUC = zosyn$AUC) %>% 
  mutate(hit = zosyn$hit)

bactrim = maier_data %>% 
  filter(chemical_name2 == "sulfamethoxazole") %>% 
  select(NT_code, species, AUC, hit) %>% 
  left_join(maier_data %>% filter(chemical_name2 == "trimethoprim") %>% select(NT_code, species, AUC, hit), by = "NT_code") %>% 
  mutate(AUC = AUC.x * AUC.y) %>% 
  mutate(hit.x = as.numeric(hit.x) ) %>% 
  mutate(hit.y = as.numeric(hit.y) ) %>% 
  mutate(hit = hit.x + hit.y)  %>% 
  mutate(hit = case_when(hit >= 1 ~ T, T ~ F))
bactrim2 = maier_data %>% 
  filter(chemical_name2 == "sulfamethoxazole") %>% 
  mutate(chemical_name2 = "sulfamethoxazole_trimethoprim") %>% 
  mutate(AUC = bactrim$AUC)%>% 
  mutate(hit = bactrim$hit)

maier_data = rbind.data.frame(maier_data, zosyn2, bactrim2)
maier_data$species = gsub(" ", ".\\", maier_data$species)

# remove one E.coli strain and one B.fragilis since there are two strains for each species screened in vitro
# we filter out the second strain, keeping only the first in terms of NT_code number 

maier_data = maier_data %>% filter(NT_code != "NT5078") %>% filter(NT_code != "NT5033")

# calculate response scores for species 

tblcounts = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/tblcounts_master_table_deid_MSKCC_9674_Duke_500.csv")
tbleuclidean_distance = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tbleuclidean_distance_10clusters.csv")
tblsample = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/tblsample_cohort_master_table_deid_MSKCC_9674_Duke_500.csv")
tblsample_discovery = tblsample %>% 
  filter(set == "discovery")

abundance_by_species = tblcounts %>% 
  filter(oligos_id %in% tblsample_discovery$oligos_id) %>% 
  group_by(oligos_id, species) %>% 
  summarise(.,total_frequency = sum(count/count_total)) %>% 
  spread(species, total_frequency, fill=0) 

colnames(abundance_by_species) = gsub("\\]", "", colnames(abundance_by_species))
colnames(abundance_by_species) = gsub("\\[", "", colnames(abundance_by_species))
colnames(abundance_by_species)[-1] = gsub("_", ".", colnames(abundance_by_species)[-1])

# select species with a minimum abundance threshold of 10^-4 and present in at least 10% of the samples

select_sp <- function(fc){
  threshold = fc >= 10^-4
  threshold = factor(threshold, level = c(TRUE, FALSE))
  t_th = table(threshold)/sum(table(threshold))
  if (t_th[["TRUE"]] >= 0.1){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

abundance_by_species_filtered = abundance_by_species %>% 
  left_join(tblsample_discovery %>% select(oligos_id, cluster_assignment), by = "oligos_id") %>% 
  select(oligos_id, match(unique(maier_data$species), colnames(abundance_by_species), nomatch=0), cluster_assignment) 

abundance_by_species_filtered = abundance_by_species_filtered[,c(1,which(apply(abundance_by_species_filtered[,-c(1,36)],2,select_sp) == TRUE) + 1, 36)]

species_profile_by_cluster = abundance_by_species_filtered [,-1] %>% 
  group_by(cluster_assignment) %>%
  summarise_all(mean) %>% 
  data.frame()

# source("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Script/github_deposit/paradigm_example_Fig2_bc.R")
## will replace source with the table containing bacteria response score 

microbiomeFeatureList_invivo = colnames(species_profile_by_cluster)[-1]
drugList_invivo = variables[which(variables %in% maier_data$chemical_name2)]
response_score_invivo_comparison = data.frame(exposure_name = drugList_invivo) 

for (i in 1:length(microbiomeFeatureList_invivo)){
  
  feature_profile_by_cluster_i = species_profile_by_cluster[[i]]
  feature_profile_name = microbiomeFeatureList_invivo[i]
  
  response_score_aggregate = c() 
  
  for (drug in drugList_invivo){
    
    response_score = bacteria_response_score(self_coef = self_coefficient_matrix, 
                                             attractor_coef = attractor_coefficient_matrix, 
                                             exposure_name = drug,
                                             feature_profile_by_cluster = feature_profile_by_cluster_i, 
                                             euclidean_distance = tbleuclidean_distance, 
                                             d = 0, 
                                             numCluster = 10)
    
    
    response_score_aggregate = c(response_score_aggregate, response_score)
    
  }
  response_score_invivo_comparison[feature_profile_name] = response_score_aggregate
}

response_score_invivo_comparison_long = melt(response_score_invivo_comparison, id = "exposure_name")

antibiotic_names = c("aztreonam", "cefepime", "ciprofloxacin", "meropenem", "metronidazole", "piperacillin_tazobactam",
                     "sulfamethoxazole_trimethoprim")

response_score_invivo_comparison_long = response_score_invivo_comparison_long %>%
  select(exposure_name, species = variable, response_score = value) %>% 
  left_join(maier_data %>% select(exposure_name = chemical_name2,
                                  species = species,
                                  AUC = AUC, 
                                  hit = hit), by = c("exposure_name", "species") ) %>% 
  mutate(drug_class = case_when(exposure_name %in% antibiotic_names ~ "antibiotics", 
                                TRUE ~ "others")) %>% 
  mutate(significance_char = case_when(hit == TRUE ~ "Significance", 
                                       TRUE ~ "Non-sign.")) %>% 
  mutate(negative_score = case_when(response_score >= 0 ~ "positive", 
                                    TRUE ~ "negative"))

# perform Fisher's Exact Test 

table(response_score_invivo_comparison_long$hit, response_score_invivo_comparison_long$negative_score)
fisher.test( table(response_score_invivo_comparison_long$hit, response_score_invivo_comparison_long$negative_score) )

# plot Fig. 3e 

ggplot(response_score_invivo_comparison_long, aes(y = AUC, x = response_score)) + 
  geom_point(aes(col = factor(hit), shape = factor(drug_class))) + 
  geom_smooth(data = response_score_invivo_comparison_long[response_score_invivo_comparison_long$hit == T,], method = "lm") +
  scale_color_manual(values = c("black", "red")) +
  scale_shape_manual(values = c("circle", "square")) +
  theme_classic() + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  scale_y_continuous(limits = c(0,1.5)) +
  xlab("Bacteria response scores \n (present study)") +
  ylab("in vitro growth (AUC) \n (Maier et al. 2018")

# plot Fig. 3f

ggplot(response_score_invivo_comparison_long, 
       aes(x = significance_char, y = response_score, group = significance_char, 
           col = significance_char, shape = factor(drug_class))) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  scale_color_manual(values = c("black", "red")) +
  scale_shape_manual(values = c("circle", "square")) +
  theme_classic() + stat_compare_means() +
  ylab("Bacteria response scores") +
  ylab("In vitro inhibition")

# perform t.test to see whether bacteria response scores are significantly different from 0 

t.test(response_score_invivo_comparison_long$response_score[response_score_invivo_comparison_long$hit == T])
t.test(response_score_invivo_comparison_long$response_score[response_score_invivo_comparison_long$hit == F])
