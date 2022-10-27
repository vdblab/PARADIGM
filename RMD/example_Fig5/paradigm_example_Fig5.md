---
title: 'PARADIGM example: Figure 5'
author: "Chi L. Nguyen"
date: "10/21/2022"
output:
  html_document:
    keep_md: yes
---




Table details:  
'tbldrugs' List of drug exposures and date of exposures for the entire MSKCC patient cohort (beween day -14 to 14 relative to HCT).  
'tblantibiotics_duke' List of antibiotics and date of exposures for the Duke patient cohort (between day -14 to 14 relative to HCT).  
'tblsample' Sample characteristics.  
'tblpatient' Patient charactestics, including discovery and validation set stratification.  
'tblresponse_score' Response scores for 4 microbiome features of interest. 


```r
tbldrugs = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tbldrugs_master_table_MSKCC_1198patients_day-14to14_102422.csv")
head(tbldrugs)
```

```
##   PatientID exposure_name exposure_day_relative_to_hct
## 1  FMT.0042   fludarabine                           -5
## 2  FMT.0042   fludarabine                           -4
## 3  FMT.0042   fludarabine                           -3
## 4  FMT.0042   fludarabine                           -2
## 5  FMT.0042     melphalan                           -2
## 6  FMT.0042      ursodiol                           -5
```

```r
length(unique(tbldrugs$PatientID))
```

```
## [1] 1198
```

```r
range(tbldrugs$exposure_day_relative_to_hct)
```

```
## [1] -14  14
```

```r
tblantibiotics_duke = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tblantibiotics_master_table_Duke_142patients_day-14to14_102422.csv")
head(tblantibiotics_duke)
```

```
##                                                           PatientID
## 1 pt_with_samples_13.aes.13_13.aes.30_13.aes.95_13.aes.96_13.tp.317
## 2 pt_with_samples_13.aes.13_13.aes.30_13.aes.95_13.aes.96_13.tp.317
## 3 pt_with_samples_13.aes.13_13.aes.30_13.aes.95_13.aes.96_13.tp.317
## 4 pt_with_samples_13.aes.13_13.aes.30_13.aes.95_13.aes.96_13.tp.317
## 5 pt_with_samples_13.aes.13_13.aes.30_13.aes.95_13.aes.96_13.tp.317
## 6 pt_with_samples_13.aes.13_13.aes.30_13.aes.95_13.aes.96_13.tp.317
##                   exposure_name exposure_day_relative_to_hct
## 1 sulfamethoxazole_trimethoprim                           -5
## 2 sulfamethoxazole_trimethoprim                           -4
## 3 sulfamethoxazole_trimethoprim                           -3
## 4 sulfamethoxazole_trimethoprim                           -2
## 5                 ciprofloxacin                           -5
## 6                 ciprofloxacin                           -4
```

```r
tblsample = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tblsample_cohort_master_table_deid_MSKCC_9167_Duke_473_post_filter_102422.csv")
dim(tblsample)
```

```
## [1] 9640   16
```

```r
tblpatient = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tblpatient_cohort_characteristics_master_table_deid_MSKCC_778_discovery_423_validation_142_Duke_102422.csv")
table(tblpatient$set)
```

```
## 
##       discovery      validation validation_Duke 
##             778             423             142
```

```r
head(tblpatient)
```

```
##   PatientID sex         intensity         simplesource disease_simple ci
## 1  FMT.0049   M          Ablative      PBSC unmodified            AML  5
## 2  FMT.0054   M Reduced Intensity      PBSC unmodified         Others  3
## 3       576   M       Nonablative      PBSC unmodified         Others  7
## 4       484   M          Ablative PBSC T-cell Depleted            AML  3
## 5      1010   M Reduced Intensity      PBSC unmodified            AML  1
## 6       799   M Reduced Intensity      PBSC unmodified         Others  4
##   age_range institution        set
## 1      >=65       MSKCC validation
## 2      >=65       MSKCC  discovery
## 3      >=65       MSKCC validation
## 4      >=65       MSKCC  discovery
## 5      >=65       MSKCC validation
## 6      >=65       MSKCC validation
```

```r
tblpatient_validation = tblpatient %>% 
  filter(set == "validation") 
length(unique(tblpatient_validation$PatientID))
```

```
## [1] 423
```

```r
tblresponse_score = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tblresponse_scores_4microbiome_features_2039samples_102422.csv")
head(tblresponse_score)
```

```
##                   exposure_name     Blautia Enterococcus Erysipelatoclostridium
## 1                   allopurinol  0.15873736   -0.3239883            -0.04685274
## 2           amlodipine_besylate  0.05942865    0.3280621             0.28421466
## 3 antithymocyte_globulin_rabbit  0.37822969   -1.0575307            -0.64055918
## 4                    aprepitant  1.65869058    0.5978437             0.07304849
## 5                     aztreonam -0.95657956    0.5924673             1.00549970
## 6                      baclofen -0.30682585    0.4586939            -0.29415818
##   simpson_reciprocal
## 1         0.05927848
## 2        -0.01740609
## 3         0.38638829
## 4         0.40101685
## 5        -0.53525022
## 6        -0.17278183
```

Pre-specified parameters. In this study, we investigated drug exposures between day -14 to 14 relative to HCT, and their associations with microbiome trajectories following drug exposures up until a month afterwards (day 14 to 45 relative to HCT). 


```r
drug_timepoint = -14:14 
microbiome_timepoint = 14:45 
```

We first explored *Enterococcus*-specific patient response scores. A negative score indicates that a given patient's drug exposure profiles between day -14 to 14 relative to HCT is associated with *Enterococcus* inhibition, and a positive score indicates that a given patient's drug exposure profiles between day -14 to 14 relative to HCT is associated with *Enterococcus* expansion. 
 

```r
patient_level_response_score_Enterococcus = tbldrugs %>% 
    left_join(tblresponse_score %>% select(exposure_name, all_of("Enterococcus") ), by = "exposure_name") %>% 
    filter(exposure_day_relative_to_hct %in% drug_timepoint) %>%
    group_by(PatientID) %>%
    summarise(feature_name_rc = across(all_of("Enterococcus"), as.numeric)) %>% 
    summarise(rc = sum(feature_name_rc)) %>% 
    mutate(rc_norm = as.numeric( scale(rc) ) ) %>% 
    filter(PatientID %in% tblpatient_validation$PatientID) 
```

```
## `summarise()` has grouped output by 'PatientID'. You can override using the
## `.groups` argument.
```

```r
gghistogram(data = patient_level_response_score_Enterococcus, x = "rc_norm") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  scale_y_continuous(limits = c(0,60), breaks = seq(0,60,by=10) ) +
  xlab("Patient-specific Enterococcus response scores") + 
  ylab("Number of patients") +
  geom_segment(aes(x = -0.5, y = 55, xend = -2, yend = 55), arrow = arrow(length = unit(0.2, "cm")), 
               color = "light blue") + 
  annotate(geom = "text", x = -1.3, y = 55.2, label = "Enterococcus-inhibiting\n association", 
           color = "light blue", size = 4) +
  geom_segment(aes(x = 1, y = 55, xend = 2.55, yend = 55), arrow = arrow(length = unit(0.2, "cm")), 
               color = "dark blue") +
  annotate(geom = "text", x = 1.7, y = 55.2, label = "Enterococcus-promoting\n association", 
           color = "dark blue", size = 4) 
```

```
## Warning: Using `bins = 30` by default. Pick better value with the argument
## `bins`.
```

![](https://github.com/ChiLNguyen/PARADIGM/blob/060c9074971ea8ff90d8dedb61d67d85e9aaed8f/RMD/example_Fig5/figures/Fig5b-1.png)<!-- -->

Let's explore whether *Enterococcus*-specific patient response scores (based solely on drug exposures between day -14 to 14 relative to HCT) are predictive of future microbiome trajectories (based on stool samples collected between day 14 to 45 relative to HCT). 


```r
genus_abundance_correlation_Enterococcus = tblsample %>% 
    filter(day_relative_to_hct %in% microbiome_timepoint) %>% 
    group_by(PatientID) %>% 
    summarise(genus_abundance = as.numeric( across(all_of("Enterococcus"), median)) ) %>% 
    select(PatientID, genus_abundance) %>% 
    filter(PatientID %in% patient_level_response_score_Enterococcus$PatientID) %>% 
    left_join(patient_level_response_score_Enterococcus, by = "PatientID") %>% 
    mutate(feature_name = "Enterococcus")

ggplot(data = genus_abundance_correlation_Enterococcus, aes(x=rc_norm, y = log10(genus_abundance+0.0001), fill = feature_name, col = feature_name)) +
  geom_point() + theme_classic() + geom_smooth(method = "lm") + stat_cor() +
  scale_color_manual(values = c("#129246"), name = "Feature name") +
  scale_fill_manual(values = c("#129246"), name = "Feature name") + 
  scale_y_continuous(limits = c(-4,1)) +
  xlab("Patient-specific\n bacteria response scores") +
  ylab("Genus abundance OR\n Alpha-diversity")  
```

```
## `geom_smooth()` using formula 'y ~ x'
```

![](https://github.com/ChiLNguyen/PARADIGM/blob/060c9074971ea8ff90d8dedb61d67d85e9aaed8f/RMD/example_Fig5/figures/genus_abundance_prediction_Enterococcus_example-1.png)<!-- -->

To reproduce Fig 5c, we looped the script through all four features of interest. 


```r
list_of_microbiome_features = c("Blautia", "Enterococcus", "Erysipelatoclostridium", "simpson_reciprocal")
genus_abundance_aggregate = data.frame()
patient_level_response_score_aggregate = data.frame(PatientID = tblpatient_validation$PatientID)

for (m in list_of_microbiome_features){
  
  patient_level_response_score = tbldrugs %>% 
    filter(PatientID %in% tblpatient_validation$PatientID) %>% 
    left_join(tblresponse_score %>% select(exposure_name, all_of(m) ), by = "exposure_name") %>% 
    filter(exposure_day_relative_to_hct %in% drug_timepoint) %>%
    group_by(PatientID) %>%
    summarise(feature_name_rc = across(all_of(m), as.numeric)) %>% 
    summarise(rc = sum(feature_name_rc)) %>% 
    mutate(rc_norm = as.numeric( scale(rc) ) ) 
  
  patient_level_response_score_aggregate = patient_level_response_score_aggregate %>% 
    left_join(patient_level_response_score, by = "PatientID")
  
  genus_abundance_correlation = tblsample %>%
    filter(PatientID %in% patient_level_response_score$PatientID) %>% 
    filter(day_relative_to_hct %in% microbiome_timepoint) %>% 
    group_by(PatientID) %>% 
    summarise(genus_abundance = as.numeric( across(all_of(m), median)) ) %>% 
    select(PatientID, genus_abundance) %>% 
    left_join(patient_level_response_score, by = "PatientID") %>% 
    mutate(feature_name = m)
  
  print(m)
  print(length(unique(patient_level_response_score$PatientID)))
  print(length(unique(genus_abundance_correlation$PatientID)))
  
  genus_abundance_aggregate = rbind.data.frame(genus_abundance_aggregate, genus_abundance_correlation)
  
}
```

```
## `summarise()` has grouped output by 'PatientID'. You can override using the
## `.groups` argument.
```

```
## [1] "Blautia"
## [1] 423
## [1] 237
```

```
## `summarise()` has grouped output by 'PatientID'. You can override using the
## `.groups` argument.
```

```
## [1] "Enterococcus"
## [1] 423
## [1] 237
```

```
## `summarise()` has grouped output by 'PatientID'. You can override using the
## `.groups` argument.
```

```
## [1] "Erysipelatoclostridium"
## [1] 423
## [1] 237
```

```
## `summarise()` has grouped output by 'PatientID'. You can override using the
## `.groups` argument.
```

```
## [1] "simpson_reciprocal"
## [1] 423
## [1] 237
```

```r
colnames(patient_level_response_score_aggregate) = c("PatientID", "Blautia_rc", "Blautia_rc_norm",
                                                     "Enterococcus_rc", "Enterococcus_rc_norm", 
                                                     "Erysipelatoclostridium_rc", "Erysipelatoclostridium_rc_norm",
                                                     "simpson_reciprocal_rc", "simpson_reciprocal_rc_norm")
```


```r
color_genus_abundance = c("#EC9B96", "#129246", "#7A6920", "black")
names(color_genus_abundance) = list_of_microbiome_features

ggplot(data = genus_abundance_aggregate, aes(x=rc_norm, y = log10(genus_abundance+0.0001), fill = feature_name, col = feature_name)) +
  geom_point() + theme_classic() + geom_smooth(method = "lm") + stat_cor() +
  scale_color_manual(values = color_genus_abundance, name = "Feature name") +
  scale_fill_manual(values = color_genus_abundance, name = "Feature name") + 
  scale_y_continuous(limits = c(-5,2)) +
  xlab("Patient-specific\n bacteria response scores") +
  ylab("Genus abundance OR\n Alpha-diversity")  
```

```
## `geom_smooth()` using formula 'y ~ x'
```

![](https://github.com/ChiLNguyen/PARADIGM/blob/060c9074971ea8ff90d8dedb61d67d85e9aaed8f/RMD/example_Fig5/figures/Fig5c_MSKCC-1.png)<!-- -->

Similarly, we could calculate the patient response scores for the Duke validation cohort. For this cohort, we only evaluated antibiotics exposures. 


```r
tblpatient_duke = tblpatient %>% 
  filter(set == "validation_Duke") 
length(unique(tblpatient_duke$PatientID))
```

```
## [1] 142
```

```r
genus_abundance_aggregate_duke = data.frame()
patient_level_response_score_aggregate_duke = data.frame(PatientID = tblpatient_duke$PatientID)

for (m in list_of_microbiome_features){
  
  patient_level_response_score = tblantibiotics_duke %>% 
    filter(PatientID %in% tblpatient_duke$PatientID) %>%
    left_join(tblresponse_score %>% select(exposure_name, all_of(m) ), by = "exposure_name") %>% 
    filter(exposure_day_relative_to_hct %in% drug_timepoint) %>%
    group_by(PatientID) %>%
    summarise(feature_name_rc = across(all_of(m), as.numeric)) %>% 
    summarise(rc = sum(feature_name_rc)) %>% 
    mutate(rc_norm = as.numeric( scale(rc) ) )
  
  patient_level_response_score_aggregate_duke = patient_level_response_score_aggregate_duke %>% 
    left_join(patient_level_response_score, by = "PatientID")
  
  genus_abundance_correlation = tblsample %>% 
    filter(PatientID %in% patient_level_response_score$PatientID) %>% 
    filter(day_relative_to_hct %in% microbiome_timepoint) %>% 
    group_by(PatientID) %>% 
    summarise(genus_abundance = as.numeric( across(all_of(m), median)) ) %>% 
    select(PatientID, genus_abundance) %>% 
    left_join(patient_level_response_score, by = "PatientID") %>% 
    mutate(feature_name = m)
  
  print(m)
  print(length(unique(patient_level_response_score$PatientID)))
  print(length(unique(genus_abundance_correlation$PatientID)))
  
  genus_abundance_aggregate_duke = rbind.data.frame(genus_abundance_aggregate_duke, genus_abundance_correlation)
  
}
```

```
## `summarise()` has grouped output by 'PatientID'. You can override using the
## `.groups` argument.
```

```
## [1] "Blautia"
## [1] 141
## [1] 73
```

```
## `summarise()` has grouped output by 'PatientID'. You can override using the
## `.groups` argument.
```

```
## [1] "Enterococcus"
## [1] 141
## [1] 73
```

```
## `summarise()` has grouped output by 'PatientID'. You can override using the
## `.groups` argument.
```

```
## [1] "Erysipelatoclostridium"
## [1] 141
## [1] 73
```

```
## `summarise()` has grouped output by 'PatientID'. You can override using the
## `.groups` argument.
```

```
## [1] "simpson_reciprocal"
## [1] 141
## [1] 73
```

```r
colnames(patient_level_response_score_aggregate_duke) = c("PatientID", "Blautia_rc", "Blautia_rc_norm",
                                                     "Enterococcus_rc", "Enterococcus_rc_norm", 
                                                     "Erysipelatoclostridium_rc", "Erysipelatoclostridium_rc_norm",
                                                     "simpson_reciprocal_rc", "simpson_reciprocal_rc_norm")

ggplot(data = genus_abundance_aggregate_duke, aes(x=rc_norm, y = log10(genus_abundance+0.0001), fill = feature_name, col = feature_name)) +
  geom_point() + theme_classic() + geom_smooth(method = "lm") + stat_cor() +
  scale_color_manual(values = color_genus_abundance, name = "Feature name") +
  scale_fill_manual(values = color_genus_abundance, name = "Feature name") + 
  scale_y_continuous(limits = c(-5,2)) +
  xlab("Patient-specific\n bacteria response scores") +
  ylab("Genus abundance OR\n Alpha-diversity")  
```

```
## `geom_smooth()` using formula 'y ~ x'
```

![]([paradigm_example_Fig5_files/figure-html/Fig5c_Duke-1.png](https://github.com/ChiLNguyen/PARADIGM/blob/060c9074971ea8ff90d8dedb61d67d85e9aaed8f/RMD/example_Fig5/figures/Fig5c_Duke-1.png))<!-- -->
