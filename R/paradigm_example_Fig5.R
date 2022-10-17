library(ggrepel)
library(ggpubr)
library(survival)

#' Source script calculates drug-taxon response scores. 

## source("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Script/github_deposit/paradigm_example_Fig2_bc.R")
## will replase source with the table containing bacteria response score 

'%notin%' = Negate('%in%')

tbldrugs = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tbldrugs_master_table_MSKCC_1227patients_day-14to14_101122.csv")

tblsample = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tblsample_cohort_master_table_deid_MSKCC_9674_Duke_500_101222.csv")

tblpatient = read.csv("tblpatient_cohort_characteristics_master_table_deid_MSKCC_1227_Duke_142.csv")
df_clinical_outcomes = read.csv("clinical_annotation_frozen_set_with_bacteria_response_scores_landmark_1227patients_set_stratification.csv")
tblpatient_validation = tblpatient %>% 
  filter(set == "validation") %>% 
  filter(simplesource != "PBSC T-cell Depleted") %>% 
  left_join(df_clinical_outcomes %>% select(PatientID, follow_up_days, follow_up_months, ci, event_type, age), by = "PatientID")

patient_level_response_score_calculation <- function(drug_exposure_table, 
                                                     bacteria_response_score_table, 
                                                     taxonomic_table, feature_name, 
                                                     drug_timepoint, microbiome_timepoint, cohort_table, 
                                                     categorical_yn, categorical_variable_names){
  
  patient_level_response_score = drug_exposure_table %>% 
    left_join(bacteria_response_score_table %>% select(exposure_name, feature_name ), by = "exposure_name") %>% 
    filter(exposure_day_relative_to_hct %in% drug_timepoint) %>%
    group_by(PatientID) %>%
    summarise(rc = as.numeric( across(all_of(feature_name), sum) ) ) %>% 
    mutate(rc_norm = as.numeric( scale(rc) ) ) %>% 
    filter(PatientID %in% cohort_table$PatientID) 
  
  if (categorical_yn == T){
    
    categorical_df = cohort_table %>% 
      select(PatientID, all_of(categorical_variable_names) ) %>% 
      pivot_longer(., cols = categorical_variable_names) %>% 
      mutate(value2 = gsub(" ", ".", value)) %>% 
      mutate(exposure_name = paste(name, value2, sep = "_")) %>% 
      left_join(bacteria_response_score_table %>% select(exposure_name, rc_feature = all_of(feature_name)), by = "exposure_name") %>% 
      mutate(days_on = length(drug_timepoint)) %>% 
      mutate(rc_categorical = rc_feature * days_on) %>% 
      group_by(PatientID) %>% 
      summarise(rc_categorical = sum(rc_categorical))
    
    patient_level_response_score = patient_level_response_score %>% 
      left_join(categorical_df, by = "PatientID") %>% 
      mutate(rc = rc + rc_categorical) %>% 
      mutate(rc_norm = scale(rc)) %>% 
      select(PatientID, rc, rc_norm)
    
  }
  
  genus_abundance_correlation = taxonomic_table %>% 
    filter(day_relative_to_hct %in% microbiome_timepoint) %>% 
    group_by(PatientID) %>% 
    summarise(genus_abundance = as.numeric( across(all_of(feature_name), median)) ) %>% 
    select(PatientID, genus_abundance) %>% 
    filter(PatientID %in% patient_level_response_score$PatientID) %>% 
    left_join(patient_level_response_score, by = "PatientID")
  
  return(list(patient_level_response_score, genus_abundance_correlation))
}

## survival analysis, will be excluded from the final script release 

list_of_microbiome_features = c("Blautia", "Enterococcus", "Erysipelatoclostridium", "simpson_reciprocal")
genus_abundance_aggregate = data.frame()
cox_single_results = data.frame()

for (m in list_of_microbiome_features){
  
  patient_level_response_score_result = patient_level_response_score_calculation(drug_exposure_table = tbldrugs, 
                                                                                 bacteria_response_score_table = response_score_feature_aggregate, 
                                                                                 taxonomic_table = tblsample , feature_name = m, 
                                                                                 drug_timepoint = -14:14 , microbiome_timepoint = 14:45,
                                                                                 cohort_table = tblpatient_validation, 
                                                                                 categorical_yn = F, 
                                                                                 categorical_variable_names = c("intensity", "disease_simple", "simplesource"))
  genus_abundance_by_response_score = patient_level_response_score_result[[2]]
  genus_abundance_by_response_score$feature_name = m
  genus_abundance_aggregate = rbind.data.frame(genus_abundance_aggregate, genus_abundance_by_response_score)
  
  patient_level_response_score_for_survival = patient_level_response_score_result[[1]]
  df_validation_rc = df_clinical_outcomes %>%
    left_join(patient_level_response_score_for_survival, by = "PatientID")
  
  cox_model = coxph(Surv(follow_up_months, event_type)~rc_norm+
                      simplesource+intensity+age+sex+disease_simple,
                    data = df_validation_rc)
  
  cox_single_results = rbind.data.frame(cox_single_results, cbind.data.frame(coef(summary(cox_model)),exp(confint(cox_model)), feature_name = m ))
  
}

## plot Fig. 5c (MSKCC)
ggplot(data = genus_abundance_aggregate, aes(x=rc_norm, y = log10(genus_abundance+0.0001), fill = feature_name, col = feature_name)) +
  geom_point() + theme_classic() + geom_smooth(method = "lm") + stat_cor() +
  scale_color_manual(values = c("#EC9B96", "#129246", "#7A6920", "black")) +
  scale_fill_manual(values = c("#EC9B96", "#129246", "#7A6920", "black")) + 
  scale_y_continuous(limits = c(-4,1))

cox_single_results$label = paste(rownames(cox_single_results ), cox_single_results$feature_name, sep = "_")
colnames(cox_single_results) = c("coef", "exp_coef", "se_coef", "z", "pvalue", "lower_95", "upper_95", "feature_name", "label")

## plot Fig. 5d (MSKCC) 

ggplot(data=cox_single_results[c(1,9,17,25),], aes(x=label, y=exp_coef, ymin=lower_95, ymax=upper_95, fill = label, col = label)) +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Label") + ylab("Mean (95% CI)") +
  scale_color_manual(values = c("#EC9B96", "#129246", "#7A6920", "black")) +
  scale_fill_manual(values = c("#EC9B96", "#129246", "#7A6920", "black")) + 
  theme_bw() +  # use a white background 
  theme(legend.position = "none") +
  geom_text(aes(label = round(pvalue, digits = 4) ), vjust = -2, ) + 
  geom_text(aes(label = round(exp_coef, digits = 4), vjust = 2)) 

p.adjust(cox_single_results$pvalue[c(1,9,17,25)], method = "BH")

## plot Fig. 5b

patient_level_response_score_result = patient_level_response_score_calculation(drug_exposure_table = tbldrugs, 
                                                                               bacteria_response_score_table = response_score_feature_aggregate, 
                                                                               taxonomic_table = tblsample , feature_name = "Enterococcus", 
                                                                               drug_timepoint = -14:14 , microbiome_timepoint = 14:45,
                                                                               cohort_table = tblpatient_validation, 
                                                                               categorical_yn = F, 
                                                                               categorical_variable_names = c("intensity", "disease_simple", "simplesource"))
patient_level_response_score_Enteroccus = patient_level_response_score_result[[1]]

gghistogram(data = patient_level_response_score_Enteroccus, x = "rc_norm") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  scale_y_continuous(limits = c(0,60), breaks = seq(0,60,by=10) )

## analysis for Duke cohort

clinical_annotation_duke = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/Duke_patient_cohort_142patients.csv")
tblantibiotics_duke = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tblantibiotics_master_table_Duke_142patients_day-14to14.csv")
duke_sample = read.csv("Duke_sample_cohort_142patients_513samples.csv")

unique_duke_mrn = unique(duke_sample$patient_uid)
unique_duke_deid = c()

for(mrn in unique_duke_mrn){
  
  samples_linked_to_mrn = duke_sample$sampleid[duke_sample$patient_uid == mrn]
  deid = paste("pt_with_samples", paste(samples_linked_to_mrn, collapse = "_"), sep = "_")
  unique_duke_deid = c(unique_duke_deid, deid)
}

unique_duke_deid_df = cbind.data.frame(patient_uid = unique_duke_mrn, PatientID = unique_duke_deid)

duke_sample = duke_sample %>% 
  left_join(unique_duke_deid_df, by = "patient_uid") %>% 
  select(patient_id = patient_uid, PatientID) %>% 
  filter(!duplicated(.))

clinical_annotation_duke = clinical_annotation_duke %>% 
  left_join(duke_sample, by = "patient_id")
  
genus_abundance_aggregate_duke = data.frame()
cox_single_results_duke = data.frame()

for (m in list_of_microbiome_features){
  
  patient_level_response_score_result = patient_level_response_score_calculation(drug_exposure_table = tblantibiotics_duke, 
                                                                                 bacteria_response_score_table = response_score_feature_aggregate, 
                                                                                 taxonomic_table = tblsample , feature_name = m, 
                                                                                 drug_timepoint = -14:14 , microbiome_timepoint = 14:45,
                                                                                 cohort_table = clinical_annotation_duke, 
                                                                                 categorical_yn = F, 
                                                                                 categorical_variable_names = c("intensity", "disease_simple", "simplesource"))
  genus_abundance_by_response_score = patient_level_response_score_result[[2]]
  genus_abundance_by_response_score$feature_name = m
  genus_abundance_aggregate_duke = rbind.data.frame(genus_abundance_aggregate_duke, genus_abundance_by_response_score)
  
  patient_level_response_score_for_survival = patient_level_response_score_result[[1]]
  df_validation_rc = clinical_annotation_duke %>%
    left_join(patient_level_response_score_for_survival, by = "PatientID")
  
  cox_model = coxph(Surv(follow_up_months, event_type)~rc_norm+
                      simplesource+conditioning_intensity+age+sex+disease_simple2,
                    data = df_validation_rc)
  
  cox_single_results_duke = rbind.data.frame(cox_single_results_duke, cbind.data.frame(coef(summary(cox_model)),exp(confint(cox_model)), feature_name = m ))
  
}

## plot Fig. 5c (Duke)
ggplot(data = genus_abundance_aggregate_duke, aes(x=rc_norm, y = log10(genus_abundance+0.0001), fill = feature_name, col = feature_name)) +
  geom_point() + theme_classic() + geom_smooth(method = "lm") + stat_cor() +
  scale_color_manual(values = c("#EC9B96", "#129246", "#7A6920", "black")) +
  scale_fill_manual(values = c("#EC9B96", "#129246", "#7A6920", "black")) + 
  scale_y_continuous(limits = c(-4,1))

cox_single_results_duke$label = paste(rownames(cox_single_results_duke ), cox_single_results_duke$feature_name, sep = "_")
colnames(cox_single_results_duke) = c("coef", "exp_coef", "se_coef", "z", "pvalue", "lower_95", "upper_95", "feature_name", "label")

## plot Fig. 5d (Duke) 

ggplot(data=cox_single_results_duke[c(1,9,17,25),], aes(x=label, y=exp_coef, ymin=lower_95, ymax=upper_95, fill = label, col = label)) +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Label") + ylab("Mean (95% CI)") +
  scale_color_manual(values = c("#EC9B96", "#129246", "#7A6920", "black")) +
  scale_fill_manual(values = c("#EC9B96", "#129246", "#7A6920", "black")) + 
  theme_bw() +  # use a white background 
  theme(legend.position = "none") +
  geom_text(aes(label = round(pvalue, digits = 4) ), vjust = -2, ) + 
  geom_text(aes(label = round(exp_coef, digits = 4), vjust = 2)) 

p.adjust(cox_single_results_duke$pvalue[c(1,9,17,25)], method = "BH")
