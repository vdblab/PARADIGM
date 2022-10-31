---
title: 'PARADIGM example: Fig4'
author: "Chi L. Nguyen"
date: "10/28/2022"
output:
  html_document:
    keep_md: yes
---



StrainPhlAn (https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-4) were used to identify strains of a species of interest reconstructed directly from the shotgun metagenomic samples, and to build phylogenetic tree showing strain evolution and strain phylogenetic distance among these ssamples. 


```r
strain_dynamics = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tblstrain_dynamics.csv")

tblresponse_score = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tblresponse_scores_4microbiome_features_2039samples_102422.csv")
```

We observed the strain-level divergence over time of sample collection. 


```r
species_scale = c("Blautia_coccoides" = "#EC9B96", 
                  "Blautia_producta" = "#c97671", 
                  "Enterococcus_faecalis" = "#129246", 
                  "Enterococcus_faecium" = "#0A572A", 
                  "Erysipelatoclostridium_ramosum" = "#7A6920")

ggplot(data = strain_dynamics, aes(x = day.x, y = log10(phylo_dist), 
                                        group = strain_name, fill = strain_name, color = strain_name )) + 
  geom_point() + geom_smooth(method = "loess") + 
  scale_fill_manual(values = species_scale) +
  scale_color_manual(values = species_scale) +
  theme_classic() +
  scale_y_continuous(limits = c(-7,1), breaks = seq(-6,1,by=2), name = "Phylogenetic distance (log10)") + 
  scale_x_continuous(limits = c(-15,15), breaks = seq(-14,14,by=7), name = "Day relative to HCT") +
  facet_wrap(~strain_name, nrow = 1) +
  theme(legend.position = "none")
```

```
## `geom_smooth()` using formula 'y ~ x'
```

![](https://github.com/vdblab/PARADIGM/blob/540fc16ef7c5ff808d52f6a25d32d03ca37746f2/RMD/example_Fig4/figures/Fig4a_top-1.png)<!-- -->

We could also investigate strain-level divergence as a function of species relative abundance. 


```r
ggplot(data = strain_dynamics, aes(x = log10(species + 0.0001), y = log10(phylo_dist), 
                                        group = strain_name, fill = strain_name, color = strain_name )) + 
  geom_point() + geom_smooth(method = "loess") + 
  scale_fill_manual(values = species_scale) +
  scale_color_manual(values = species_scale) +
  theme_classic() +
  scale_y_continuous(limits = c(-7,1), breaks = seq(-6,1,by=2), name = "Phylogenetic distance (log10)") + 
  scale_x_continuous(name = "Species relative abundance (log10)") +
  facet_wrap(~strain_name, nrow = 1) +
  theme(legend.position = "none")
```

```
## `geom_smooth()` using formula 'y ~ x'
```

![](https://github.com/vdblab/PARADIGM/blob/540fc16ef7c5ff808d52f6a25d32d03ca37746f2/RMD/example_Fig4/figures/Fig4a_bottom-1.png)<!-- -->

*Enterococcus faecium* is a species of particular interest for further investigation of the associations drug exposures and strain evolution and divergence, due to the association between *E. faecium* and important allo-HCT outcomes such as infection and graft-versus-host disease. We compared the phylogenetic distance between subsequently collected samples, for sample pairs under antibiotics exposures and sample pairs not under antibiotics exposures. 

For antibiotics exposures, we considered all seven antibiotics investigated in this study.
For non-antibiotics exposures, we considered the top 7 non-antibiotics most strongly associated with genus *Enterococcus* relative abundance based on the response scores predicted by PARADIGM (in Fig. 3b). We excluded ursodiol from this list due to the fact that most data points are exposed to ursodiol. 


```r
antibiotics = c("ciprofloxacin", "aztreonam", "piperacillin_tazobactam", 
                "meropenem", "metronidazole", "sulfamethoxazole_trimethoprim", 
                "cefepime")
non_antibiotics = tblresponse_score$exposure_name[order(abs(tblresponse_score$Enterococcus), decreasing = T )]
non_antibiotics = non_antibiotics[!(non_antibiotics %in% antibiotics )]
non_antibiotics = non_antibiotics[-which(non_antibiotics == "ursodiol")]
non_antibiotics = non_antibiotics[1:7]

antibiotics
```

```
## [1] "ciprofloxacin"                 "aztreonam"                    
## [3] "piperacillin_tazobactam"       "meropenem"                    
## [5] "metronidazole"                 "sulfamethoxazole_trimethoprim"
## [7] "cefepime"
```

```r
non_antibiotics
```

```
## [1] "diphenoxylate_atropine"        "polyethylene_glycol"          
## [3] "levothyroxine_sodium"          "fentanyl"                     
## [5] "methotrexate"                  "antithymocyte_globulin_rabbit"
## [7] "cyclosporine"
```

```r
strain_dynamics_efaecium = strain_dynamics %>% 
  filter(strain_name == "Enterococcus_faecium") %>% 
  mutate(abx_exposed = rowSums(across(all_of(antibiotics))) ) %>% 
  mutate(abx_exposed_yn = case_when(abx_exposed > 0 ~ "Exposed",
                                    TRUE ~ "Not exposed")) %>% 
  mutate(nonabx_exposed = rowSums(across(all_of(non_antibiotics))) ) %>% 
  mutate(nonabx_exposed_yn = case_when(nonabx_exposed > 0 ~ "Exposed",
                                    TRUE ~ "Not exposed"))

ggplot(data = strain_dynamics_efaecium, aes(x = factor(abx_exposed_yn, level = c("Not exposed", "Exposed")), 
                                                     y = log10(phylo_dist), 
                                                     fill = factor(abx_exposed_yn, level = c("Not exposed", "Exposed"))) ) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  scale_fill_manual(values = c("light grey", "red")) +
  scale_y_continuous(limits = c(-7,1), breaks = seq(-6,1,by=2), name = "Phylogenetic distance (log10)") +
  theme_classic() +
  stat_compare_means() +
  xlab("Antibiotics") +
  theme(legend.position = "none")
```

![](https://github.com/vdblab/PARADIGM/blob/540fc16ef7c5ff808d52f6a25d32d03ca37746f2/RMD/example_Fig4/figures/Fig4b-1.png)<!-- -->

```r
ggplot(data = strain_dynamics_efaecium, aes(x = factor(nonabx_exposed_yn, level = c("Not exposed", "Exposed")), 
                                                     y = log10(phylo_dist), 
                                                     fill = factor(nonabx_exposed_yn, level = c("Not exposed", "Exposed"))) ) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  scale_fill_manual(values = c("light grey", "blue")) +
  scale_y_continuous(limits = c(-7,1), breaks = seq(-6,1,by=2), name = "Phylogenetic distance (log10)") +
  theme_classic() +
  stat_compare_means() +
  xlab("Non-antibiotics") +
  theme(legend.position = "none")
```

![](https://github.com/vdblab/PARADIGM/blob/540fc16ef7c5ff808d52f6a25d32d03ca37746f2/RMD/example_Fig4/figures/Fig4b-2.png)<!-- -->
