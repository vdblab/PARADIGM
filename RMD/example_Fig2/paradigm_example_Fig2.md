PARADIGM example: Figure 2
================
Chi L. Nguyen
10/17/2022

Cohort characteristics are outlined in ‘tblsample’ table; ASV counts and
taxonomic profiles are outlined in ‘tblcounts’ table. ‘tblsample’ table
contains samples from the entire study, which is divided into the MSKCC
discovery cohort, the MSKCC validation cohort and the Duke validation
cohort.

``` r
tblcounts = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tblcounts_master_table_deid_MSKCC_9674_Duke_500_101222.csv")

tblsample = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tblsample_cohort_master_table_deid_MSKCC_9674_Duke_500_101222.csv")

head(tblcounts)
```

    ##   sampleid      oligos_id asv_key count count_total  kingdom     phylum   class
    ## 1    1533D 1533D..pool938   asv_1    16       21089 Bacteria Firmicutes Bacilli
    ## 2    1533E 1533E..pool938   asv_1     7       11685 Bacteria Firmicutes Bacilli
    ## 3    1533F 1533F..pool938   asv_1     7       12841 Bacteria Firmicutes Bacilli
    ## 4    1533G 1533G..pool938   asv_1    12       14669 Bacteria Firmicutes Bacilli
    ## 5    1533H 1533H..pool938   asv_1   799        5671 Bacteria Firmicutes Bacilli
    ## 6    1581A 1581A..pool938   asv_1    39       17800 Bacteria Firmicutes Bacilli
    ##             order          family        genus              species   color
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

``` r
head(tblsample)
```

    ##   sampleid      oligos_id
    ## 1    1000A 1000A..pool483
    ## 2    1000B 1000B..pool483
    ## 3    1000C 1000C..pool483
    ## 4    1000D 1000D..pool483
    ## 5    1000E 1000E..pool483
    ## 6     1001  1001..pool535
    ##                                                               PatientID
    ## 1                                                                  1000
    ## 2                                                                  1000
    ## 3                                                                  1000
    ## 4                                                                  1000
    ## 5                                                                  1000
    ## 6 pt_with_samples_1001_1002_1003_1004_1005_1006_1007_1008_1048_1121_152
    ##   cluster_assignment day_relative_to_hct        set simpson_reciprocal
    ## 1                  2                  -9  discovery           13.38430
    ## 2                  4                  -4  discovery            6.24706
    ## 3                  5                   6  discovery            2.08790
    ## 4                  5                   9  discovery            2.39139
    ## 5                  5                  13  discovery            1.48186
    ## 6                 NA                  -6 validation            8.48545
    ##       Blautia Enterococcus Erysipelatoclostridium institution      tsne1
    ## 1 0.406969925  0.000000000            0.005193935       MSKCC -0.7462970
    ## 2 0.023031435  0.000000000            0.161220044       MSKCC  0.0951663
    ## 3 0.000000000  0.001073784            0.000000000       MSKCC  0.2074600
    ## 4 0.001466583  0.001780851            0.000000000       MSKCC  0.2446640
    ## 5 0.000000000  0.003923877            0.000000000       MSKCC  0.2555170
    ## 6 0.322351421  0.000000000            0.136950904       MSKCC -0.8626910
    ##        tsne2 shotgun_available tsne1_shotgun tsne2_shotgun
    ## 1 -0.0408929                No            NA            NA
    ## 2  0.1467790                No            NA            NA
    ## 3  0.4915400                No            NA            NA
    ## 4  0.2800610               Yes    -0.7709889     -0.351067
    ## 5  0.4738010                No            NA            NA
    ## 6  0.0937450                No            NA            NA

Taxonomic and cluster characteristics as plotted in Figure 2 are from
the MSKCC discovery cohort.

``` r
tblsample_discovery = tblsample %>% 
  filter(set == "discovery")
```

Since the entire cohort is large, we provided the tSNE coordinates we
pre-calculated using the R package ‘Rtsne’. tSNE analysis was performed
on the pairwise Bray-Curtis beta-diversity matrix (which was
pre-calculated using QIIME). We first identified the most abundant genus
per sample.

``` r
genus_abundance_tbl = tblcounts %>% 
  filter(oligos_id %in% tblsample_discovery$oligos_id) %>% 
  group_by(oligos_id, genus) %>% 
  summarise(.,total_frequency = sum(count/count_total)) %>% 
  spread(genus, total_frequency, fill = 0) %>% 
  group_by(oligos_id)
```

    ## `summarise()` has grouped output by 'oligos_id'. You can override using the
    ## `.groups` argument.

``` r
dominant_genus = colnames(genus_abundance_tbl[,-1])[apply(genus_abundance_tbl[,-1], 1, which.max)]
tblsample_discovery$dominant_genus = dominant_genus
```

We then assigned a color for each genus using the pre-generated color
scheme in ‘tblcounts’ table.

``` r
color_scheme_genus = tblcounts %>% 
  select(color, genus) %>% 
  filter(!duplicated(.))
color_scheme_genus_vector = color_scheme_genus$color
names(color_scheme_genus_vector) = color_scheme_genus$genus

ggplot(tblsample_discovery, aes(x=tsne1,y=tsne2,col=factor(dominant_genus))) +
  geom_point() +
  scale_color_manual(name = "Legend", values=color_scheme_genus_vector) +
  theme_classic() +
  coord_fixed() +
  theme(legend.position = "none")
```

![](paradigm_example_Fig2_files/figure-gfm/Fig2a-1.png)<!-- -->

Taxonomic profiles (relative abundance) of shotgun sequencing data from
MetaPhlAn 3.0are outlined in a separate table. Each row is a taxonmic
classification, each column is a sample. Some samples with numeric names
will have the letter “X” preceeding the numeric names, and the character
“X” needs to be removed.

``` r
tblshotgun = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tblshotgun_MetaPhlAn_classification.csv")
colnames(tblshotgun) = gsub("X(\\d)", "\\1", colnames(tblshotgun))

tblshotgun[1:6,1:10]
```

    ##                                                                                    clade_name
    ## 1                                                                                     UNKNOWN
    ## 2                                                                                  k__Archaea
    ## 3                                                                 k__Archaea|p__Euryarchaeota
    ## 4                                              k__Archaea|p__Euryarchaeota|c__Methanobacteria
    ## 5                        k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales
    ## 6 k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae
    ##                   NCBI_tax_id 894A    801A 788B 786A 784B 759A 737B 697B
    ## 1                          -1    0 0.00000    0    0    0    0    0    0
    ## 2                        2157    0 0.03016    0    0    0    0    0    0
    ## 3                  2157|28890    0 0.03016    0    0    0    0    0    0
    ## 4           2157|28890|183925    0 0.03016    0    0    0    0    0    0
    ## 5      2157|28890|183925|2158    0 0.03016    0    0    0    0    0    0
    ## 6 2157|28890|183925|2158|2159    0 0.03016    0    0    0    0    0    0

The pre-calculated tSNE coordinates for shotgun data are outlined in
‘tblsample’ table. We selected shotgun samples in the discovery cohort,
and selected species classification.

``` r
tblshotgun_species = tblshotgun[grepl("s__", tblshotgun$clade_name) , ]
tblshotgun_species$clade_name = gsub(".*s__", "", tblshotgun_species$clade_name)
tblshotgun_species_t = t(tblshotgun_species[,-c(1,2)])
colnames(tblshotgun_species_t) = tblshotgun_species$clade_name
tblshotgun_species_relative_abundance = apply(tblshotgun_species_t, c(1,2), function(x){as.numeric(x)/100})
shotgun_discovery_ind = which(rownames(tblshotgun_species_relative_abundance) %in% tblsample_discovery$sampleid)
tblshotgun_species_relative_abundance = tblshotgun_species_relative_abundance[shotgun_discovery_ind,]

tblshotgun_species_relative_abundance[1:6, 1:6]
```

    ##           Methanobrevibacter_oralis Methanobrevibacter_smithii
    ## 682D                              0                          0
    ## 276D                              0                          0
    ## FMT.0099O                         0                          0
    ## FMT.0099A                         0                          0
    ## FMT.0094A                         0                          0
    ## FMT.0093Z                         0                          0
    ##           Methanosphaera_stadtmanae Actinobaculum_sp_oral_taxon_183
    ## 682D                              0                               0
    ## 276D                              0                               0
    ## FMT.0099O                         0                               0
    ## FMT.0099A                         0                               0
    ## FMT.0094A                         0                               0
    ## FMT.0093Z                         0                               0
    ##           Actinomyces_cardiffensis Actinomyces_denticolens
    ## 682D                             0                       0
    ## 276D                             0                       0
    ## FMT.0099O                        0                       0
    ## FMT.0099A                        0                       0
    ## FMT.0094A                        0                       0
    ## FMT.0093Z                        0                       0

We then identified the most abundant species per shotgun sample.

``` r
shotgun_dominant_species = colnames(tblshotgun_species_relative_abundance)[apply(tblshotgun_species_relative_abundance, 1, which.max)]
shotgun_dominant_species_df = cbind.data.frame(sampleid = rownames(tblshotgun_species_relative_abundance),
                                               dominant_species = shotgun_dominant_species)
```

We extracted the tSNE coordinates for shotgun sample, assigned a color
for each species using the pre-generated color scheme in ‘tblcounts’
table.

``` r
tblshotgun_plot_tsne = tblsample %>% 
  filter(shotgun_available == "Yes") %>% 
  filter(set == "discovery") %>% 
  left_join(shotgun_dominant_species_df, by = "sampleid")

color_scheme_species = tblcounts %>% 
  select(color_shotgun, species) %>% 
  filter(!duplicated(.))
color_scheme_species_vector = color_scheme_species$color_shotgun
names(color_scheme_species_vector) = color_scheme_species$species

ggplot(tblshotgun_plot_tsne , aes(x=tsne1_shotgun,y=tsne2_shotgun,col=factor(dominant_species))) +
  geom_point() +
  scale_color_manual(name = "Legend", values=color_scheme_species_vector) +
  theme_classic() +
  coord_fixed() +
  theme(legend.position = "none")
```

![](paradigm_example_Fig2_files/figure-gfm/Fig2b-1.png)<!-- -->

We pre-defined kmeans cluster by applying kmeans() functions to the
pairwise Bray-Curtis beta-diversity matrix for samples in the MSKCC
discovery cohort. Kmeans cluster assignments are outlined in the
‘tblsample’ table.

``` r
colors_asv_beta <- c("#A50026","#D73027","#F46D43","#C2A5CF","#3690C0",
                     "#FDAE61","#5AAE61","#B3E2CD","#762A83","#1B7837")

ggplot(tblsample_discovery, aes(x=tsne1,y=tsne2,col=factor(cluster_assignment))) +
  geom_point() +
  scale_color_manual(name = "Legend", values=colors_asv_beta) +
  theme_classic() +
  coord_fixed() +
  theme(legend.position = "none")
```

![](paradigm_example_Fig2_files/figure-gfm/Fig2c-1.png)<!-- -->

We then visualized the taxonomic characteristics of each cluster by
plotting a heatmap of the top 20 most abundant genera (based on 16S
sequencing data) in the MSKCC discovery cohort. First, we identified the
top 20 most abundant genera and their relative abundance per sample.

``` r
col_choice = names(sort(colMeans(genus_abundance_tbl[,-1]),decreasing=T)[1:20])

top20_genus_tbl = tblsample_discovery %>% 
  arrange(cluster_assignment) %>% 
  select(oligos_id, cluster_assignment) %>%  
  left_join(genus_abundance_tbl %>% select(oligos_id, all_of(col_choice) ), by = "oligos_id")
```

We then visualized this information with a heatmap.

``` r
mat_colors <- list(group = colors_asv_beta)
names(mat_colors$group) <- c(1:10)

hm = data.frame(top20_genus_tbl[,-c(1:2)])
hm = apply(hm, c(1,2), as.numeric)
rownames(hm) = top20_genus_tbl$oligos_id

group = data.frame(group = top20_genus_tbl$cluster_assignment)
rownames(group) = top20_genus_tbl$oligos_id

pheatmap.type(log10(hm+0.00001),annRow=group,show_rownames=FALSE, annotation_colors = mat_colors)
```

![](paradigm_example_Fig2_files/figure-gfm/Fig2d-1.png)<!-- -->

We could also visualize the taxonomic characteristics of each cluster by
plotting a heatmap of the top 20 most abundant species (based on shotgun
metageonomic sequencing) in the MSKCC discovery cohort. First, we
identified the top 20 most abundant species in the shotgun data.

``` r
col_choice_shotgun = names(sort(colMeans(tblshotgun_species_relative_abundance),decreasing=T)[1:20])

tblshotgun_species_relative_abundance = data.frame(tblshotgun_species_relative_abundance)
tblshotgun_species_relative_abundance$sampleid = rownames(tblshotgun_species_relative_abundance)

top20_species_shotgun_tbl = tblsample_discovery %>% 
  filter(sampleid %in% tblshotgun_species_relative_abundance$sampleid) %>% 
  arrange(cluster_assignment) %>% 
  select(sampleid, cluster_assignment) %>%  
  left_join(tblshotgun_species_relative_abundance %>% 
              select(sampleid, all_of(col_choice_shotgun) ), by = "sampleid")
```

We then visualized this information with a heatmap.

``` r
mat_colors <- list(group = colors_asv_beta)
names(mat_colors$group) <- c(1:10)

hm = data.frame(top20_species_shotgun_tbl)
hm = hm[,-c(1,2)]
hm = apply(hm, c(1,2), as.numeric)
colnames(hm) = colnames(top20_species_shotgun_tbl)[-c(1,2)]
rownames(hm) = top20_species_shotgun_tbl$sampleid

group = data.frame(group = top20_species_shotgun_tbl$cluster_assignment)
rownames(group) = top20_species_shotgun_tbl$sampleid

pheatmap.type(log10(hm+0.00001),annRow=group,show_rownames=FALSE, annotation_colors = mat_colors)
```

![](paradigm_example_Fig2_files/figure-gfm/Fig2e-1.png)<!-- -->

Alpha-diversity of each sample (based on 16S sequencing data) are
outlined in the ‘tblsample’ table. We could plot alpha-diversity of each
sample per cluster in a boxplot.

``` r
ggplot(tblsample_discovery, aes(y = simpson_reciprocal, x = factor(cluster_assignment), fill = factor(cluster_assignment))) + 
  geom_boxplot() +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(name = "Simpson reciprocal", breaks=seq(0,60,by=20), limits = c(0,60)) +
  scale_x_discrete(name = "Cluster index", breaks=seq(1,10,by=1)) +
  scale_fill_manual(name = "Cluster", values = colors_asv_beta) +
  theme_classic() +
  geom_hline(yintercept = median(tblsample_discovery$simpson_reciprocal), linetype = 2) 
```

![](paradigm_example_Fig2_files/figure-gfm/Fig2f-1.png)<!-- -->

To examine cluster dynamics over time, we calculated cluster frequency
in a weekly interval between day -28 and 27 relative to HCT.

``` r
tblsample_discovery_time = tblsample_discovery %>% 
  filter(day_relative_to_hct %in% -28:27) %>% 
  mutate(day_relative_to_hct_f = factor(day_relative_to_hct, levels = c(-28:27) ))

count_tab <- as.matrix(table(tblsample_discovery_time$day_relative_to_hct_f, tblsample_discovery_time$cluster_assignment) )

relative_frequency_time_matrix = data.frame()

for (i in seq(1,50,by=7)){
  tab = count_tab[i:(i+6),]
  sum = colSums(tab)
  rel_freq = sum/sum(tab)
  relative_frequency_time_matrix = rbind.data.frame(relative_frequency_time_matrix, rel_freq)
}

relative_frequency_time_matrix = as.data.frame(t(relative_frequency_time_matrix))
colnames(relative_frequency_time_matrix) <- c("-28 to -21", "-20 to -14", "-13 to -7", "-6 to 0", "1 to 7", "8 to 14", "15 to 21", "22 to 28")
rownames(relative_frequency_time_matrix) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
m2 <- melt(cbind(relative_frequency_time_matrix, ind = rownames(relative_frequency_time_matrix)), id.vars = c('ind'))
```

We visualized cluster dynamics over time by a stacked barplot.

``` r
ggplot(m2, aes(x = variable, y = value, fill = factor(ind))) + 
  geom_bar(stat = "identity") +
  ggtitle("Relative frequency k = 10") +
  scale_y_continuous("Percent") +
  scale_x_discrete("") +
  theme_classic() + 
  scale_fill_manual(values = colors_asv_beta) +
  theme(legend.position = "none", axis.text.x = element_text(angle=45,hjust=1))
```

![](paradigm_example_Fig2_files/figure-gfm/Fig2g-1.png)<!-- -->
