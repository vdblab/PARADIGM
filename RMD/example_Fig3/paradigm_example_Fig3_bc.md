---
title: 'PARADIGM example: Figure 3b-c'
author: "Chi L. Nguyen"
date: "10/17/2022"
output:
  html_document:
    keep_md: yes
---



Source script runs elastic net regularized regression to learn the coefficients of each drug exposure to cluster self and attractor transitions, then convert cluster-drug associations into drug-taxon response scores.   

Table details:  
'tbldaily' Training data of microbiome cluster dynamics.  
'tblpatient' Patient characteristics.  
'tbleuclidean_distance' A 10x10 matrix of pair-wise cluster euclidean distance.  
'tblcounts' Counts and taxonomic classifications of ASVs.  
'tblsample' Sample characteristics.  


```r
tbldaily = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tbldaily_sampling_PARADIGM_input_post_filter_2039samples_102422.csv")
head(tbldaily)
```

```
##   PatientID sampleid.x sampleid.y n10.x n10.y day.x day.y dday allopurinol
## 1      1648      1648A      1648B     3     2    -3    -2    1       FALSE
## 2      1648      1648B      1648C     2     3    -2    -1    1       FALSE
## 3      1648      1648C      1648D     3     6    -1     0    1       FALSE
## 4      1648      1648D      1648E     6     6     0     1    1       FALSE
## 5  FMT.0251  FMT.0251A  FMT.0251B     2     1    -1     0    1       FALSE
## 6  FMT.0251  FMT.0251C  FMT.0251D     1     1     3     4    1       FALSE
##   amlodipine_besylate antithymocyte_globulin_rabbit aprepitant aztreonam
## 1               FALSE                         FALSE      FALSE     FALSE
## 2               FALSE                         FALSE      FALSE     FALSE
## 3               FALSE                         FALSE      FALSE     FALSE
## 4               FALSE                         FALSE      FALSE     FALSE
## 5               FALSE                         FALSE      FALSE     FALSE
## 6               FALSE                         FALSE      FALSE     FALSE
##   baclofen busulfan cefepime ciprofloxacin cyclophosphamide cyclosporine
## 1    FALSE    FALSE    FALSE         FALSE            FALSE        FALSE
## 2    FALSE    FALSE    FALSE          TRUE            FALSE        FALSE
## 3    FALSE    FALSE    FALSE          TRUE            FALSE        FALSE
## 4    FALSE    FALSE    FALSE          TRUE            FALSE        FALSE
## 5    FALSE    FALSE    FALSE          TRUE            FALSE        FALSE
## 6    FALSE    FALSE    FALSE         FALSE             TRUE        FALSE
##   diphenoxylate_atropine docusate_sodium enoxaparin entecavir estradiol
## 1                  FALSE           FALSE      FALSE     FALSE     FALSE
## 2                  FALSE           FALSE      FALSE     FALSE     FALSE
## 3                  FALSE           FALSE      FALSE     FALSE     FALSE
## 4                  FALSE           FALSE      FALSE     FALSE     FALSE
## 5                  FALSE           FALSE      FALSE     FALSE     FALSE
## 6                  FALSE           FALSE      FALSE     FALSE     FALSE
##   famotidine fentanyl fludarabine furosemide gabapentin hydralazine
## 1      FALSE    FALSE        TRUE      FALSE      FALSE       FALSE
## 2      FALSE    FALSE        TRUE      FALSE      FALSE       FALSE
## 3      FALSE    FALSE       FALSE      FALSE      FALSE       FALSE
## 4      FALSE    FALSE       FALSE      FALSE      FALSE        TRUE
## 5      FALSE    FALSE       FALSE      FALSE      FALSE       FALSE
## 6      FALSE    FALSE       FALSE       TRUE      FALSE       FALSE
##   hydrocortisone hydromorphone hydroxyzine insulin isavuconazonium_sulfate
## 1          FALSE         FALSE       FALSE   FALSE                   FALSE
## 2          FALSE         FALSE       FALSE   FALSE                   FALSE
## 3          FALSE         FALSE       FALSE   FALSE                   FALSE
## 4          FALSE         FALSE       FALSE   FALSE                   FALSE
## 5          FALSE         FALSE       FALSE   FALSE                   FALSE
## 6          FALSE         FALSE       FALSE   FALSE                   FALSE
##   labetalol letermovir levetiracetam levothyroxine_sodium loperamide loratadine
## 1     FALSE      FALSE          TRUE                FALSE      FALSE      FALSE
## 2     FALSE      FALSE         FALSE                FALSE      FALSE      FALSE
## 3     FALSE      FALSE         FALSE                FALSE      FALSE      FALSE
## 4     FALSE      FALSE         FALSE                FALSE      FALSE      FALSE
## 5     FALSE      FALSE         FALSE                 TRUE      FALSE      FALSE
## 6     FALSE      FALSE         FALSE                 TRUE       TRUE      FALSE
##   melphalan meropenem mesna methotrexate methylprednisolone metoclopramide
## 1     FALSE     FALSE FALSE        FALSE              FALSE          FALSE
## 2     FALSE     FALSE FALSE        FALSE              FALSE          FALSE
## 3     FALSE     FALSE FALSE        FALSE              FALSE          FALSE
## 4     FALSE     FALSE FALSE        FALSE              FALSE          FALSE
## 5     FALSE     FALSE FALSE        FALSE              FALSE          FALSE
## 6     FALSE     FALSE  TRUE        FALSE              FALSE          FALSE
##   metoprolol metronidazole morphine_sulfate mycophenolate_mofetil olanzapine
## 1      FALSE         FALSE            FALSE                 FALSE      FALSE
## 2      FALSE         FALSE            FALSE                 FALSE      FALSE
## 3      FALSE         FALSE            FALSE                 FALSE      FALSE
## 4      FALSE         FALSE            FALSE                 FALSE      FALSE
## 5       TRUE         FALSE            FALSE                 FALSE      FALSE
## 6       TRUE         FALSE            FALSE                 FALSE      FALSE
##   oxycodone palifermin palonosetron piperacillin_tazobactam polyethylene_glycol
## 1     FALSE      FALSE        FALSE                   FALSE               FALSE
## 2     FALSE      FALSE         TRUE                   FALSE               FALSE
## 3     FALSE      FALSE        FALSE                   FALSE               FALSE
## 4     FALSE      FALSE         TRUE                   FALSE               FALSE
## 5     FALSE      FALSE        FALSE                   FALSE               FALSE
## 6     FALSE      FALSE         TRUE                    TRUE               FALSE
##   posaconazole prochlorperazine senna simethicone sirolimus sucralfate
## 1        FALSE             TRUE FALSE       FALSE     FALSE      FALSE
## 2        FALSE             TRUE FALSE       FALSE     FALSE      FALSE
## 3        FALSE             TRUE FALSE       FALSE     FALSE      FALSE
## 4        FALSE             TRUE FALSE       FALSE     FALSE      FALSE
## 5        FALSE            FALSE FALSE       FALSE     FALSE      FALSE
## 6        FALSE            FALSE FALSE       FALSE     FALSE      FALSE
##   sulfamethoxazole_trimethoprim tacrolimus tamsulosin thiotepa ursodiol
## 1                         FALSE       TRUE      FALSE    FALSE     TRUE
## 2                         FALSE       TRUE      FALSE    FALSE     TRUE
## 3                         FALSE       TRUE      FALSE    FALSE     TRUE
## 4                         FALSE       TRUE      FALSE    FALSE     TRUE
## 5                         FALSE      FALSE      FALSE    FALSE     TRUE
## 6                         FALSE      FALSE      FALSE    FALSE     TRUE
##   voriconazole zolpidem
## 1        FALSE    FALSE
## 2        FALSE    FALSE
## 3        FALSE    FALSE
## 4        FALSE    FALSE
## 5        FALSE    FALSE
## 6        FALSE    FALSE
```

```r
dim(tbldaily)
```

```
## [1] 2039   70
```

```r
length(unique(tbldaily$PatientID))
```

```
## [1] 454
```

```r
tblpatient = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tblpatient_cohort_characteristics_master_table_deid_MSKCC_778_discovery_423_validation_142_Duke_102422.csv")
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
table(tblpatient$set)
```

```
## 
##       discovery      validation validation_Duke 
##             778             423             142
```

```r
tbleuclidean_distance = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tbleuclidean_distance_10clusters_post_filter_7454samples_102422.csv")
head(tbleuclidean_distance)
```

```
##          V1        V2        V3       V4       V5       V6       V7       V8
## 1  0.000000  7.289069 10.609251 13.51186 15.58425 11.28030 19.38569 12.21990
## 2  7.289069  0.000000  7.660139 16.40821 19.61082 13.76219 21.46263 17.77651
## 3 10.609251  7.660139  0.000000 18.48578 21.60996 10.02171 22.80795 20.09254
## 4 13.511856 16.408213 18.485782  0.00000 12.17649 14.34784 17.13110 11.25883
## 5 15.584251 19.610819 21.609956 12.17649  0.00000 16.55367 18.38770 10.67781
## 6 11.280301 13.762195 10.021712 14.34784 16.55367  0.00000 19.33602 14.31118
##         V9      V10
## 1 21.75277 31.17886
## 2 24.00877 33.17864
## 3 26.06393 34.45165
## 4 10.86221 29.31250
## 5 19.61815 29.32778
## 6 22.19882 31.02081
```

```r
tblcounts = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tblcounts_master_table_deid_MSKCC_9167_Duke_473_post_filter_102422.csv")
head(tblcounts)
```

```
##   sampleid      oligos_id asv_key count count_total  kingdom     phylum   class
## 1    1671K 1671K..pool971   asv_1   217       70762 Bacteria Firmicutes Bacilli
## 2    1671L 1671L..pool971   asv_1    54       74588 Bacteria Firmicutes Bacilli
## 3    1671M 1671M..pool971   asv_1   128       78377 Bacteria Firmicutes Bacilli
## 4    1671N 1671N..pool971   asv_1   116       58586 Bacteria Firmicutes Bacilli
## 5    1671O 1671O..pool971   asv_1    41       81884 Bacteria Firmicutes Bacilli
## 6    1671Q 1671Q..pool971   asv_1    66       80089 Bacteria Firmicutes Bacilli
##              ordr          family        genus              species   color
## 1 Lactobacillales Enterococcaceae Enterococcus Enterococcus_faecium #129246
## 2 Lactobacillales Enterococcaceae Enterococcus Enterococcus_faecium #129246
## 3 Lactobacillales Enterococcaceae Enterococcus Enterococcus_faecium #129246
## 4 Lactobacillales Enterococcaceae Enterococcus Enterococcus_faecium #129246
## 5 Lactobacillales Enterococcaceae Enterococcus Enterococcus_faecium #129246
## 6 Lactobacillales Enterococcaceae Enterococcus Enterococcus_faecium #129246
##   color_shotgun
## 1       #0A572A
## 2       #0A572A
## 3       #0A572A
## 4       #0A572A
## 5       #0A572A
## 6       #0A572A
```

```r
tblsample = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tblsample_cohort_master_table_deid_MSKCC_9167_Duke_473_post_filter_102422.csv")
head(tblsample)
```

```
##   sampleid      oligos_id PatientID cluster_assignment day_relative_to_hct
## 1    1000A 1000A..pool483      1000                  2                  -9
## 2    1000B 1000B..pool483      1000                  4                  -4
## 3    1000C 1000C..pool483      1000                  5                   6
## 4    1000D 1000D..pool483      1000                  5                   9
## 5    1000E 1000E..pool483      1000                  5                  13
## 6    1001A 1001A..pool483      1001                 NA                  -3
##          set simpson_reciprocal     Blautia Enterococcus Erysipelatoclostridium
## 1  discovery           13.38430 0.406969925  0.000000000            0.005193935
## 2  discovery            6.24706 0.023031435  0.000000000            0.161220044
## 3  discovery            2.08790 0.000000000  0.001073784            0.000000000
## 4  discovery            2.39139 0.001466583  0.001780851            0.000000000
## 5  discovery            1.48186 0.000000000  0.003923877            0.000000000
## 6 validation           14.71820 0.371913367  0.001793351            0.065250379
##   institution      tsne1      tsne2 shotgun_available tsne1_shotgun
## 1       MSKCC -0.7462970 -0.0408929                no            NA
## 2       MSKCC  0.0951663  0.1467790                no            NA
## 3       MSKCC  0.2074600  0.4915400                no            NA
## 4       MSKCC  0.2446640  0.2800610               yes    -0.7709889
## 5       MSKCC  0.2555170  0.4738010                no            NA
## 6       MSKCC -0.6486000  0.0798572                no            NA
##   tsne2_shotgun
## 1            NA
## 2            NA
## 3            NA
## 4     -0.351067
## 5            NA
## 6            NA
```

We created two data frames, one containing the response variable (cluster self/attractor transtions) and one containing the predictors (time of sample collection and drug exposures). 


```r
cluster_dynamic = tbldaily[,1:5]
head(cluster_dynamic) 
```

```
##   PatientID sampleid.x sampleid.y n10.x n10.y
## 1      1648      1648A      1648B     3     2
## 2      1648      1648B      1648C     2     3
## 3      1648      1648C      1648D     3     6
## 4      1648      1648D      1648E     6     6
## 5  FMT.0251  FMT.0251A  FMT.0251B     2     1
## 6  FMT.0251  FMT.0251C  FMT.0251D     1     1
```

```r
environmental_exposures = tbldaily[,c(6,9:ncol(tbldaily))]
head(environmental_exposures)
```

```
##   day.x allopurinol amlodipine_besylate antithymocyte_globulin_rabbit
## 1    -3       FALSE               FALSE                         FALSE
## 2    -2       FALSE               FALSE                         FALSE
## 3    -1       FALSE               FALSE                         FALSE
## 4     0       FALSE               FALSE                         FALSE
## 5    -1       FALSE               FALSE                         FALSE
## 6     3       FALSE               FALSE                         FALSE
##   aprepitant aztreonam baclofen busulfan cefepime ciprofloxacin
## 1      FALSE     FALSE    FALSE    FALSE    FALSE         FALSE
## 2      FALSE     FALSE    FALSE    FALSE    FALSE          TRUE
## 3      FALSE     FALSE    FALSE    FALSE    FALSE          TRUE
## 4      FALSE     FALSE    FALSE    FALSE    FALSE          TRUE
## 5      FALSE     FALSE    FALSE    FALSE    FALSE          TRUE
## 6      FALSE     FALSE    FALSE    FALSE    FALSE         FALSE
##   cyclophosphamide cyclosporine diphenoxylate_atropine docusate_sodium
## 1            FALSE        FALSE                  FALSE           FALSE
## 2            FALSE        FALSE                  FALSE           FALSE
## 3            FALSE        FALSE                  FALSE           FALSE
## 4            FALSE        FALSE                  FALSE           FALSE
## 5            FALSE        FALSE                  FALSE           FALSE
## 6             TRUE        FALSE                  FALSE           FALSE
##   enoxaparin entecavir estradiol famotidine fentanyl fludarabine furosemide
## 1      FALSE     FALSE     FALSE      FALSE    FALSE        TRUE      FALSE
## 2      FALSE     FALSE     FALSE      FALSE    FALSE        TRUE      FALSE
## 3      FALSE     FALSE     FALSE      FALSE    FALSE       FALSE      FALSE
## 4      FALSE     FALSE     FALSE      FALSE    FALSE       FALSE      FALSE
## 5      FALSE     FALSE     FALSE      FALSE    FALSE       FALSE      FALSE
## 6      FALSE     FALSE     FALSE      FALSE    FALSE       FALSE       TRUE
##   gabapentin hydralazine hydrocortisone hydromorphone hydroxyzine insulin
## 1      FALSE       FALSE          FALSE         FALSE       FALSE   FALSE
## 2      FALSE       FALSE          FALSE         FALSE       FALSE   FALSE
## 3      FALSE       FALSE          FALSE         FALSE       FALSE   FALSE
## 4      FALSE        TRUE          FALSE         FALSE       FALSE   FALSE
## 5      FALSE       FALSE          FALSE         FALSE       FALSE   FALSE
## 6      FALSE       FALSE          FALSE         FALSE       FALSE   FALSE
##   isavuconazonium_sulfate labetalol letermovir levetiracetam
## 1                   FALSE     FALSE      FALSE          TRUE
## 2                   FALSE     FALSE      FALSE         FALSE
## 3                   FALSE     FALSE      FALSE         FALSE
## 4                   FALSE     FALSE      FALSE         FALSE
## 5                   FALSE     FALSE      FALSE         FALSE
## 6                   FALSE     FALSE      FALSE         FALSE
##   levothyroxine_sodium loperamide loratadine melphalan meropenem mesna
## 1                FALSE      FALSE      FALSE     FALSE     FALSE FALSE
## 2                FALSE      FALSE      FALSE     FALSE     FALSE FALSE
## 3                FALSE      FALSE      FALSE     FALSE     FALSE FALSE
## 4                FALSE      FALSE      FALSE     FALSE     FALSE FALSE
## 5                 TRUE      FALSE      FALSE     FALSE     FALSE FALSE
## 6                 TRUE       TRUE      FALSE     FALSE     FALSE  TRUE
##   methotrexate methylprednisolone metoclopramide metoprolol metronidazole
## 1        FALSE              FALSE          FALSE      FALSE         FALSE
## 2        FALSE              FALSE          FALSE      FALSE         FALSE
## 3        FALSE              FALSE          FALSE      FALSE         FALSE
## 4        FALSE              FALSE          FALSE      FALSE         FALSE
## 5        FALSE              FALSE          FALSE       TRUE         FALSE
## 6        FALSE              FALSE          FALSE       TRUE         FALSE
##   morphine_sulfate mycophenolate_mofetil olanzapine oxycodone palifermin
## 1            FALSE                 FALSE      FALSE     FALSE      FALSE
## 2            FALSE                 FALSE      FALSE     FALSE      FALSE
## 3            FALSE                 FALSE      FALSE     FALSE      FALSE
## 4            FALSE                 FALSE      FALSE     FALSE      FALSE
## 5            FALSE                 FALSE      FALSE     FALSE      FALSE
## 6            FALSE                 FALSE      FALSE     FALSE      FALSE
##   palonosetron piperacillin_tazobactam polyethylene_glycol posaconazole
## 1        FALSE                   FALSE               FALSE        FALSE
## 2         TRUE                   FALSE               FALSE        FALSE
## 3        FALSE                   FALSE               FALSE        FALSE
## 4         TRUE                   FALSE               FALSE        FALSE
## 5        FALSE                   FALSE               FALSE        FALSE
## 6         TRUE                    TRUE               FALSE        FALSE
##   prochlorperazine senna simethicone sirolimus sucralfate
## 1             TRUE FALSE       FALSE     FALSE      FALSE
## 2             TRUE FALSE       FALSE     FALSE      FALSE
## 3             TRUE FALSE       FALSE     FALSE      FALSE
## 4             TRUE FALSE       FALSE     FALSE      FALSE
## 5            FALSE FALSE       FALSE     FALSE      FALSE
## 6            FALSE FALSE       FALSE     FALSE      FALSE
##   sulfamethoxazole_trimethoprim tacrolimus tamsulosin thiotepa ursodiol
## 1                         FALSE       TRUE      FALSE    FALSE     TRUE
## 2                         FALSE       TRUE      FALSE    FALSE     TRUE
## 3                         FALSE       TRUE      FALSE    FALSE     TRUE
## 4                         FALSE       TRUE      FALSE    FALSE     TRUE
## 5                         FALSE      FALSE      FALSE    FALSE     TRUE
## 6                         FALSE      FALSE      FALSE    FALSE     TRUE
##   voriconazole zolpidem
## 1        FALSE    FALSE
## 2        FALSE    FALSE
## 3        FALSE    FALSE
## 4        FALSE    FALSE
## 5        FALSE    FALSE
## 6        FALSE    FALSE
```

```r
dim(cluster_dynamic)
```

```
## [1] 2039    5
```

```r
dim(environmental_exposures)
```

```
## [1] 2039   63
```

We calculated coefficients of predictors to cluster self and attractor transitions.  
Expected output is a matrix, with predictors column-wise and coefficient values row-wise.  
The parameters are explained in the source script. 


```r
self_coefficient_matrix = self_transition(cluster_dynamic = cluster_dynamic, 
                                          environmental_exposures = environmental_exposures, 
                                          seed = 2208, 
                                          numFold = 10, 
                                          numCluster = 10)
```

```
## [1] "Processing cluster 1"
## [1] "Processing cluster 2"
## [1] "Processing cluster 3"
## [1] "Processing cluster 4"
## [1] "Processing cluster 5"
## [1] "Processing cluster 6"
## [1] "Processing cluster 7"
```

```
## Warning in nominalTrainWorkflow(x = x, y = y, wts = weights, info = trainInfo, :
## There were missing values in resampled performance measures.
```

```
## [1] "Processing cluster 8"
## [1] "Processing cluster 9"
## [1] "Processing cluster 10"
```

```
## Warning: from glmnet C++ code (error code -91); Convergence for 91th lambda
## value not reached after maxit=100000 iterations; solutions for larger lambdas
## returned
```

```r
head(self_coefficient_matrix)
```

```
##    Intercept        day.x allopurinol amlodipine_besylate
## 1 -0.1400418 -0.023793141   0.0000000                   0
## 2  0.3113556  0.000000000   0.0000000                   0
## 3 -0.4937865  0.051151730  -0.3063351                   0
## 4 -0.2590320  0.000000000   0.0000000                   0
## 5  0.6088448  0.000000000   0.0000000                   0
## 6  0.6137571  0.009113115   0.0000000                   0
##   antithymocyte_globulin_rabbit aprepitant   aztreonam  baclofen   busulfan
## 1                   -0.03576882  -2.217768 -0.00619185 0.3216409  0.0000000
## 2                    0.00000000   0.000000  0.00000000 0.0000000  0.0000000
## 3                    0.00000000   1.010212  0.00000000 0.6965204  0.3797357
## 4                    0.00000000   0.000000  0.00000000 0.0000000  0.0000000
## 5                    0.00000000   0.000000  0.00000000 0.0000000  0.0000000
## 6                    0.00000000   0.000000  0.00000000 0.0000000 -0.1228584
##     cefepime ciprofloxacin cyclophosphamide cyclosporine diphenoxylate_atropine
## 1  0.0000000    0.00000000        0.0000000            0              0.0000000
## 2  0.0000000   -0.03612195        0.0000000            0              0.0000000
## 3 -1.4367447    0.00000000       -0.5907536            0              0.9579282
## 4  0.0000000    0.00000000        0.0000000            0              0.0000000
## 5  0.0000000    0.00000000        0.0000000            0              0.0000000
## 6 -0.5416531    0.12066703        0.0000000            0              0.0000000
##   docusate_sodium  enoxaparin entecavir  estradiol famotidine   fentanyl
## 1       0.0000000  0.00000000 -0.214744 -0.8031322  0.0000000  0.4419992
## 2       0.0000000  0.00000000  0.000000  0.0000000  0.0000000  0.0000000
## 3       0.0000000 -0.17406742  0.106967  0.1984475 -0.1903401 -0.3176473
## 4       0.0000000  0.00000000  0.000000  0.0000000  0.0000000  0.0000000
## 5       0.0000000  0.00000000  0.000000  0.0000000  0.0000000  0.0000000
## 6      -0.7604281 -0.03074963  0.000000  0.0000000  0.0000000  0.0000000
##   fludarabine furosemide gabapentin hydralazine hydrocortisone hydromorphone
## 1  0.00000000 -0.4749400  0.0000000   0.0000000      0.7942646   0.000000000
## 2  0.09854763  0.0000000  0.0000000   0.0000000      0.0000000   0.000000000
## 3 -0.29416150 -0.3008423 -0.3208483  -2.1574534      0.0000000  -0.004920562
## 4  0.00000000  0.0000000  0.0000000   0.0000000      0.0000000   0.000000000
## 5  0.00000000  0.0000000  0.0000000   0.0000000      0.0000000   0.000000000
## 6  0.00000000  0.0000000  0.0000000   0.2531269      0.0000000   0.000000000
##   hydroxyzine   insulin isavuconazonium_sulfate labetalol letermovir
## 1  0.00000000 0.5139214              -0.6488319 0.0000000  0.0000000
## 2  0.00000000 0.0000000               0.0000000 0.0000000  0.0000000
## 3  0.90682755 0.2468979               0.5018957 0.0000000 -0.5058638
## 4  0.00000000 0.0000000               0.0000000 0.0000000  0.0000000
## 5  0.00000000 0.0000000               0.0000000 0.0000000  0.0000000
## 6 -0.09713641 0.0000000               0.0000000 0.5655705  0.0000000
##   levetiracetam levothyroxine_sodium loperamide loratadine melphalan  meropenem
## 1     0.3471355          0.003991994  0.0000000  0.0000000  0.000000 -0.3320321
## 2     0.0000000          0.000000000 -0.5083219  0.0000000  0.000000  0.0000000
## 3     0.0000000          0.399355833  0.0000000  0.1111664 -1.735586  0.0000000
## 4     0.0000000          0.000000000  0.0000000  0.0000000  0.000000  0.0000000
## 5     0.0000000          0.000000000  0.0000000  0.0000000  0.000000  0.0000000
## 6     0.0000000          0.000000000  0.0000000  0.0000000  0.000000  0.0000000
##         mesna methotrexate methylprednisolone metoclopramide metoprolol
## 1  1.51207404    0.3221244         -0.5149653        0.00000          0
## 2  0.00000000    0.0000000          0.0000000        0.00000          0
## 3 -0.02901782    0.1086481          0.0000000       -1.61142          0
## 4  0.00000000    0.0000000          0.0000000        0.00000          0
## 5  0.00000000    0.0000000          0.0000000        0.00000          0
## 6  0.00000000    0.0000000          0.0000000        0.00000          0
##   metronidazole morphine_sulfate mycophenolate_mofetil olanzapine  oxycodone
## 1      0.000000      -0.11803720             0.0000000   0.000000 -0.2958560
## 2      0.000000       0.00000000             0.0000000   0.000000  0.0000000
## 3     -1.340693      -0.28097359            -0.5414176  -2.004263  0.1606446
## 4      0.000000       0.06344341             0.0000000   0.000000  0.0000000
## 5      0.000000       0.00000000             0.0000000   0.000000  0.0000000
## 6     -1.259357       0.41126479             0.0000000   0.000000  0.0000000
##   palifermin palonosetron piperacillin_tazobactam polyethylene_glycol
## 1 -0.6647088    0.0000000              -0.5958084           0.7458183
## 2  0.0000000    0.0000000               0.0000000           0.0000000
## 3 -0.1641323    0.0300603              -0.7259889           0.0000000
## 4  0.0000000    0.0000000               0.0000000           0.0000000
## 5  0.0000000   -0.1548487               0.0000000           0.0000000
## 6  0.0000000    0.0000000              -0.9585926          -0.2680720
##   posaconazole prochlorperazine     senna simethicone  sirolimus  sucralfate
## 1   -0.1803234        0.0000000 0.0000000   0.3690371 0.30464983  0.00000000
## 2    0.0000000        0.0000000 0.0000000   0.0000000 0.00000000  0.00000000
## 3    0.2780898       -0.2131087 0.2429833  -0.7476522 0.06712467 -2.00523679
## 4    0.0000000        0.0000000 0.0000000  -0.1040165 0.00000000  0.00000000
## 5    0.0000000        0.0000000 0.0000000   0.0000000 0.00000000  0.00000000
## 6    0.0000000        0.0000000 0.0000000   0.0000000 0.00000000  0.02631096
##   sulfamethoxazole_trimethoprim tacrolimus tamsulosin   thiotepa   ursodiol
## 1                     0.0501414  0.0000000  0.0000000  0.0000000 -0.1011974
## 2                     0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
## 3                    -0.2987858 -0.1255332  0.2123955 -0.2931381  1.0068113
## 4                     0.3225567  0.0000000  0.0000000  0.0000000  0.0000000
## 5                     0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
## 6                     0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
##   voriconazole    zolpidem
## 1    0.0000000  0.00000000
## 2    0.0000000  0.00000000
## 3    0.0000000  0.40294878
## 4    0.0000000 -0.05527504
## 5    0.0000000  0.00000000
## 6   -0.3569107  0.00000000
```

```r
attractor_coefficient_matrix = attractor_transition(cluster_dynamic = cluster_dynamic, 
                                                    environmental_exposures = environmental_exposures,
                                                    euclidean_distance = tbleuclidean_distance, 
                                                    seed = 2208, 
                                                    numFold = 10, 
                                                    numCluster = 10)
```

```
## [1] "Processing cluster 1"
## [1] "Processing cluster 2"
## [1] "Processing cluster 3"
## [1] "Processing cluster 4"
## [1] "Processing cluster 5"
## [1] "Processing cluster 6"
## [1] "Processing cluster 7"
## [1] "Processing cluster 8"
## [1] "Processing cluster 9"
## [1] "Processing cluster 10"
```

```r
head(attractor_coefficient_matrix)
```

```
##   Intercept         day.x allopurinol amlodipine_besylate
## 1 -1.973970  0.0024143469  0.37127167         -0.31156259
## 2 -1.922221 -0.0459529895  0.09175195          0.00000000
## 3 -2.531545 -0.0053531321  0.27563035          0.00000000
## 4 -1.968206  0.0001008266  0.00000000          0.16078410
## 5 -2.725170  0.0000000000 -1.03915959          0.13253785
## 6 -2.379389  0.0001503296  0.08280710          0.09530755
##   antithymocyte_globulin_rabbit aprepitant   aztreonam     baclofen   busulfan
## 1                   1.383496963  1.1660303 -2.19434709 -2.840607262 -0.7449130
## 2                   0.133711930  0.5965643 -0.13393394  0.000000000  0.2228301
## 3                   0.003428886  0.0000000  0.01393778 -0.007913443  0.0000000
## 4                   0.000000000  0.4356610  0.33955588  0.223335594  0.0000000
## 5                  -0.990098414 -1.1074744  0.16548909  1.100825670 -0.2168365
## 6                   0.000000000 -0.4451633  0.83636647 -0.442272182 -0.2374264
##      cefepime ciprofloxacin cyclophosphamide cyclosporine
## 1  0.86310676    0.09603773       0.80297978  -0.81528197
## 2 -0.06782254    0.00000000       0.00000000   0.00000000
## 3  0.00000000    0.11454457       0.35293367   0.00000000
## 4 -0.34248662    0.01331914       0.00000000  -0.03579379
## 5  0.46996976    0.17709941       0.07032465   0.00000000
## 6  0.00000000    0.66225382       0.59822211  -0.02181930
##   diphenoxylate_atropine docusate_sodium  enoxaparin   entecavir  estradiol
## 1              0.9900328      0.83625990 -0.15179951  0.00000000 -0.8473275
## 2              0.0000000      0.14694179  0.33638017 -0.01405314  0.0000000
## 3             -0.1741712      0.07570974  0.01192229  0.00000000  0.7465923
## 4              0.0000000      0.00000000  0.00000000  0.39660629  0.0000000
## 5             -0.6317765     -0.70991778 -0.15483785 -0.12677465 -0.2458593
## 6              0.0000000     -0.06879213  0.00000000 -0.15837606  0.0000000
##   famotidine    fentanyl fludarabine   furosemide  gabapentin  hydralazine
## 1 -1.9107019 -2.53492779   0.2294302  0.112196797  1.00972941 -0.508936940
## 2  0.1518973  0.00000000   0.7136978  0.000000000  0.00000000  0.002964125
## 3  0.0000000  0.00000000  -0.2408313  0.099104968  0.11362060 -0.045768092
## 4 -0.2465371 -0.03541637   0.0000000  0.000000000 -0.09736823 -0.319290057
## 5  0.4484885  0.26068299   0.0946610  0.526643625  0.98911985  0.000000000
## 6 -0.2524451  0.27853502  -0.2387156 -0.007825215 -0.31858258 -0.049622492
##   hydrocortisone hydromorphone hydroxyzine     insulin isavuconazonium_sulfate
## 1      0.8536001   -0.03209942  -0.3077884 -0.05831239              -0.5953267
## 2      0.0000000    0.00000000  -0.0294153  0.00000000               0.0000000
## 3     -0.3677780   -0.07324446   0.1544093 -0.07392604               0.0000000
## 4     -0.4396051    0.00000000   0.1756145 -0.30789816              -0.1093762
## 5      1.6365809    0.00000000   0.5791671 -0.23838351              -1.0289357
## 6      0.0000000    0.00000000   0.1734844 -0.18440588              -0.5009345
##    labetalol letermovir levetiracetam levothyroxine_sodium  loperamide
## 1 -2.5392615 -0.1168265     0.7128895           0.26891560 -0.26473369
## 2  0.0000000  0.0000000     0.0000000          -0.04435581 -0.26129108
## 3  0.0000000  0.0000000     0.2049792           0.05722022 -0.04926700
## 4  0.4197026 -0.1545197    -0.2677232           0.00000000  0.00000000
## 5  0.0000000  0.3680820    -0.5340559          -0.99409279  0.76204881
## 6  0.0000000  0.1928353     0.0000000           0.00000000 -0.07753176
##    loratadine  melphalan   meropenem      mesna methotrexate methylprednisolone
## 1  0.01467884 -0.4050364 -0.47217005 -0.8307043  -1.61800707          0.0000000
## 2 -0.06323171  0.7084729 -0.01999656  0.0000000  -0.20776692          0.0000000
## 3  0.28709705  0.2022733 -0.13025924  0.0000000   0.00000000          0.0000000
## 4 -0.36193144 -0.1265850  0.22647719  0.3511030  -0.08842402          0.0000000
## 5  0.21683229  0.5122316 -0.38509758  0.6144305   0.00000000          0.4414713
## 6 -0.29394599 -0.1007561  0.27806810  0.0000000   0.00000000          0.0000000
##   metoclopramide  metoprolol metronidazole morphine_sulfate
## 1      0.9523486 -0.02076376    0.00000000      -0.20751954
## 2      0.0000000  0.20546983   -0.07028993       0.00000000
## 3     -0.7107244  0.06890582   -0.20229337       0.00000000
## 4      0.4431214 -0.05607834   -0.01864135       0.00000000
## 5      0.0000000  0.06146263    0.00000000      -0.48817721
## 6     -0.4045248 -0.31286695    0.12464604       0.04635525
##   mycophenolate_mofetil  olanzapine   oxycodone   palifermin palonosetron
## 1           -0.12044963  0.32115244  0.06881598 -0.007790792   -0.3278813
## 2            0.00000000  0.00000000  0.03756716  0.000000000    0.0000000
## 3           -0.03236915 -0.19301818  0.02730154  0.000000000    0.0000000
## 4            0.20260088 -0.03378235  0.00000000  0.346656196    0.0000000
## 5           -0.23618380 -0.17241310 -0.41130415 -0.877093317    0.4046232
## 6           -0.12219981  0.00000000  0.00000000  0.000000000   -0.1223130
##   piperacillin_tazobactam polyethylene_glycol posaconazole prochlorperazine
## 1           -0.4329072626           1.1175089    0.0000000      -0.42320657
## 2           -0.0173749318          -0.2873404    0.0000000       0.00000000
## 3           -0.2728136529           1.1077059    0.0000000       0.00000000
## 4            0.4121799755           0.1785268    0.0000000      -0.14478136
## 5           -0.0008810412          -0.6553278    0.9780014      -0.20026052
## 6           -0.2086284144          -0.3319308    0.0000000       0.05160934
##        senna simethicone   sirolimus  sucralfate sulfamethoxazole_trimethoprim
## 1 -0.2024203  -0.3157736 -0.34970592 -0.58587161                   -0.09166086
## 2  0.6319945   0.0000000  0.03447566 -0.02365071                    0.09485971
## 3  0.1706772   0.0000000  0.07682494  0.00000000                    0.33674994
## 4  0.0000000   0.0000000 -0.21648803 -0.39092338                   -0.16524739
## 5 -0.3443936   0.4351109 -0.62000720  0.97753934                   -0.99368935
## 6 -0.5297137   0.1175808 -0.23268086  0.00000000                   -0.03725320
##    tacrolimus  tamsulosin   thiotepa     ursodiol voriconazole    zolpidem
## 1 -0.42423565  0.08315223 -1.2754316 0.9449275156  -0.26396256 -0.10700922
## 2 -0.04038497  0.00000000  0.0000000 0.0000000000   0.00000000 -0.12232914
## 3  0.11470440  0.00000000  0.8813242 0.0006732827   0.00000000  0.33505767
## 4  0.01411641  0.00000000 -0.1587500 0.4119686131   0.06060163  0.00000000
## 5 -0.14957615 -0.58450964  0.6369294 0.0000000000   0.09328590  0.03051026
## 6  0.07736819  0.20889095  0.0000000 0.0000000000  -0.11206904  0.05214616
##   euclidean_distance_vector
## 1               -0.10731697
## 2               -0.04700235
## 3               -0.04800007
## 4               -0.07323456
## 5               -0.09852747
## 6               -0.05859611
```

To plot Fig. 3b and Fig.S3, we arranged drug exposures by drug class. 


```r
self_coefficient_to_plot = t(self_coefficient_matrix)
attractor_coefficient_to_plot = t(attractor_coefficient_matrix)

tblgraph_drug_exposure_classification = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tblgraph_drug_exposure_classification.csv") 

drug_group_annotation = tblgraph_drug_exposure_classification %>% select(group, ind)
rownames(drug_group_annotation) = tblgraph_drug_exposure_classification$exposure_name
```

We could now plot a heatmap for the coefficients of the elastic net regularized regression model, which indicates the magnitude and direction of association between a given exposures and cluster self transitions...


```r
set = c("#FFFFFF", RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(12, "Set3")[3:7])
set = list(group = set)
names(set$group) = unique(tblgraph_drug_exposure_classification$group)

ind_df = match(rownames(drug_group_annotation), rownames(self_coefficient_to_plot)) 

out_self = pheatmap(self_coefficient_to_plot[c(1,2,ind_df),], col = colorRampPalette(c("navy", "white", "firebrick3"), bias = 1.3)(150), 
                    annotation_row = drug_group_annotation, annotation_colors = set, 
                    cluster_cols = F, cluster_rows = F, breaks = seq(-3,4,by=0.05), main = "Self Transition")
```

![](https://github.com/ChiLNguyen/PARADIGM/blob/3f72929276630d26495616e4ad520ec43fd9f69f/RMD/example_Fig3/figures/FigS3_self_manual_export.png)<!-- -->

And cluster attractor transitions. 


```r
out_attractor = pheatmap(attractor_coefficient_to_plot[c(1,2,ind_df),], col = colorRampPalette(c("navy", "white", "firebrick3"), bias = 1.3)(150), 
                         annotation_row = drug_group_annotation, annotation_colors = set,
                         cluster_cols = F, cluster_rows = F, breaks = seq(-3,4,by=0.05), main = "Attractor Transition")    
```

![](https://github.com/ChiLNguyen/PARADIGM/blob/3f72929276630d26495616e4ad520ec43fd9f69f/RMD/example_Fig3/figures/FigS3_attractor_manual_export.png)<!-- -->

We trained the model and learned the associations between a drug exposure and a given cluster self and attractor transitions. A negative coefficient value indicates that a drug exposure is associated with decreased cluster self/attractor transition probability, and a positive coefficient value indicates that a drug exposure is associated with increased cluster self/attracotr transition probability.  
  
To convert drug-cluster associations into drug-taxon associations, we calculated bacteria response scores, which measure the magnitude and direction of drug-taxon associations. Deriving the bacteria response score depends on cluster transition probability matrix and cluster-specific mean measurement values of the taxonomic feature of interest. The source script contains a function to calculate cluster transition probability matrix based on the self and attractor coefficient matrix.  
  
First, we determined the mean taxonomic relative abundance (or alpha-diversity) per cluster. Users could explore other microbiome features of interest (if available in either the 'tblsample' or the 'tblcounts' tables).  
  
In our study, we focused on microbiome features that have been associated with allo-HCT outcomes, including relative abundances of genus Blautia, Enterococcus, Erysipelatoclostridium, and alpha-diversity (Simpson reciprocal).  


```r
drugList = colnames(environmental_exposures)[-1]
microbiomeFeatureList = c("Blautia", "Enterococcus", "Erysipelatoclostridium", "simpson_reciprocal")
feature_profile_by_cluster = list()
for (m in microbiomeFeatureList){
  
  feature_profile = tblsample %>% 
    filter(set == "discovery") %>% 
    group_by(cluster_assignment) %>% 
    summarise(across(all_of(m), mean)) 
  
  feature_profile_by_cluster = append(feature_profile_by_cluster, list(feature_profile[[m]]) )
} 
```

We then calculated the bacteria response scores for 62 investigated medications across these four microbiome features. 


```r
response_score_feature_aggregate = data.frame(exposure_name = drugList)

for (i in 1:length(microbiomeFeatureList)){
  
  feature_profile_by_cluster_i = feature_profile_by_cluster[[i]]
  feature_profile_name = microbiomeFeatureList[i]
  
  response_score_aggregate = c() 
  
  for (drug in drugList){
    
    response_score = bacteria_response_score(self_coef = self_coefficient_matrix, 
                                             attractor_coef = attractor_coefficient_matrix, 
                                             exposure_name = drug,
                                             feature_profile_by_cluster = feature_profile_by_cluster_i, 
                                             euclidean_distance = tbleuclidean_distance, 
                                             d = 0, 
                                             numCluster = 10)
                                                  
    
    response_score_aggregate = c(response_score_aggregate, response_score)
                          
  }
  response_score_feature_aggregate[feature_profile_name] = response_score_aggregate
}

head(response_score_feature_aggregate)
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

The 'response_score_feature_aggregate' table contains bacteria response score values, with each row being a drug exposure and each column being a microbiome feature. A negative response score indicates that a drug exposure is associated with decreased taxon relative abundance (or alpha-diversity), and a positive response score indicates that a drug exposure is associated with increased taxon relative abundance (or alpha-diversity). 


```r
drug_level = rownames(drug_group_annotation)
response_score_feature_aggregate_to_plot = melt(response_score_feature_aggregate, id = "exposure_name")
response_score_feature_aggregate_to_plot = response_score_feature_aggregate_to_plot %>% 
  mutate(sign_coef = case_when(value < 0 ~ "-", 
                               TRUE ~ "+"))

ggplot(data = response_score_feature_aggregate_to_plot, 
       aes(x=variable, y=factor(exposure_name, levels = rev(drug_level)), fill=value )) +
  geom_tile(stat = "identity") + theme_classic() + 
  theme(plot.title=element_text(size=16,face="bold"),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank()) +
  geom_text(aes(label = sign_coef), color = "black", size = 4) +
  coord_flip() + coord_fixed() +
  scale_fill_gradientn(colours = colorRampPalette(c("navy", "white", "firebrick3"), bias = 1.5)(30), 
                       limits = c(-3, 5), breaks = seq(-2,4,by=2)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.y = element_blank(), axis.title.x = element_blank())
```

```
## Coordinate system already present. Adding new coordinate system, which will replace the existing one.
```

![](https://github.com/ChiLNguyen/PARADIGM/blob/3f72929276630d26495616e4ad520ec43fd9f69f/RMD/example_Fig3/figures/FigS4_response_score_manual_export.png)<!-- -->

As a sanity check, we correlated the bacteria response scores for alpha-diversity and for Enterococcus relative abundance. We expected and observed a significant negative correlation, indicating that drug exposures which are associated with increased Enterococcus relative abundnace, are associated with decreased alpha-diversity, and vice-versa. In allo-HCT, Enterococcus expansion typically leads to a low-diversity state of Enterococcus domination and dysbiosis. 


```r
ggplot(data = response_score_feature_aggregate, 
       aes(x=Enterococcus, y=simpson_reciprocal)) +
  geom_point() +
  geom_smooth(method = "lm", linetype = "dashed") +
  stat_cor() +
  theme_classic() +
  xlab("Enterococcus \n response scores") +
  ylab("Alpha-diversity \n response scores")
```

```
## `geom_smooth()` using formula 'y ~ x'
```

![](https://github.com/ChiLNguyen/PARADIGM/blob/1068daee45b20d6bdf489d179e88fd4fe9cb6bb4/RMD/example_Fig3/figures/Fig3d-1.png)<!-- -->

