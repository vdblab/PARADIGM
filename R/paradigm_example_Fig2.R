library(tidyverse)
library(pheatmap)
library(Pigengene)

tblsample = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tblsample_cohort_master_table_deid_MSKCC_9674_Duke_500_101222.csv")
tblcounts = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tblcounts_master_table_deid_MSKCC_9674_Duke_500_101222.csv")

color_scheme_genus = tblcounts %>% 
  select(color, genus) %>% 
  filter(!duplicated(.))
color_scheme_genus_vector = color_scheme_genus$color
names(color_scheme_genus_vector) = color_scheme_genus$genus

tblsample_discovery = tblsample %>% 
  filter(set == "discovery")

## find the most abundant genus in each sample 

genus_abundance_tbl = tblcounts %>% 
  filter(oligos_id %in% tblsample_discovery$oligos_id) %>% 
  group_by(oligos_id, genus) %>% 
  summarise(.,total_frequency = sum(count/count_total)) %>% 
  spread(genus, total_frequency, fill = 0) %>% 
  group_by(oligos_id)

dominant_genus = colnames(genus_abundance_tbl[,-1])[apply(genus_abundance_tbl[,-1], 1, which.max)]
tblsample_discovery$dominant_genus = dominant_genus

## plot Fig. 2a

ggplot(tblsample_discovery, aes(x=tsne1,y=tsne2,col=factor(dominant_genus))) +
  geom_point() +
  scale_color_manual(name = "Legend", values=color_scheme_genus_vector) +
  theme_classic() +
  coord_fixed() +
  theme(legend.position = "none")

## plot Fig. 2b

tblshotgun = read.csv("tblshotgun_MetaPhlAn_classification.csv")
colnames(tblshotgun) = gsub("X(\\d)", "\\1", colnames(tblshotgun))

tblshotgun_species = tblshotgun[grepl("s__", tblshotgun$clade_name) , ]
tblshotgun_species$clade_name = gsub(".*s__", "", tblshotgun_species$clade_name)
tblshotgun_species_t = t(tblshotgun_species[,-c(1,2)])
colnames(tblshotgun_species_t) = tblshotgun_species$clade_name
tblshotgun_species_relative_abundance = apply(tblshotgun_species_t, c(1,2), function(x){as.numeric(x)/100})
shotgun_discovery_ind = which(rownames(tblshotgun_species_relative_abundance) %in% tblsample_discovery$sampleid)
tblshotgun_species_relative_abundance = tblshotgun_species_relative_abundance[shotgun_discovery_ind,]

shotgun_dominant_species = colnames(tblshotgun_species_relative_abundance)[apply(tblshotgun_species_relative_abundance, 1, which.max)]
shotgun_dominant_species_df = cbind.data.frame(sampleid = rownames(tblshotgun_species_relative_abundance),
                                               dominant_species = shotgun_dominant_species)
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

## plot Fig. 2c

colors_asv_beta <- c("#A50026","#D73027","#F46D43","#C2A5CF","#3690C0",
                     "#FDAE61","#5AAE61","#B3E2CD","#762A83","#1B7837")

ggplot(tblsample_discovery, aes(x=tsne1,y=tsne2,col=factor(cluster_assignment))) +
  geom_point() +
  scale_color_manual(name = "Legend", values=colors_asv_beta) +
  theme_classic() +
  coord_fixed() +
  theme(legend.position = "none")

## plot Fig. 2d

col_choice = names(sort(colMeans(genus_abundance_tbl[,-1]),decreasing=T)[1:20])

top20_genus_tbl = tblsample_discovery %>% 
  arrange(cluster_assignment) %>% 
  select(oligos_id, cluster_assignment) %>%  
  left_join(genus_abundance_tbl %>% select(oligos_id, all_of(col_choice) ), by = "oligos_id")

mat_colors <- list(group = colors_asv_beta)
names(mat_colors$group) <- c(1:10)

hm = data.frame(top20_genus_tbl[,-c(1:2)])
hm = apply(hm, c(1,2), as.numeric)
rownames(hm) = top20_genus_tbl$oligos_id

group = data.frame(group = top20_genus_tbl$cluster_assignment)
rownames(group) = top20_genus_tbl$oligos_id

pheatmap.type(log10(hm+0.00001),annRow=group,show_rownames=FALSE, annotation_colors = mat_colors)

## plot Fig. 2e

col_choice_shotgun = names(sort(colMeans(tblshotgun_species_relative_abundance),decreasing=T)[1:20])

tblshotgun_species_relative_abundance = data.frame(tblshotgun_species_relative_abundance)
tblshotgun_species_relative_abundance$sampleid = rownames(tblshotgun_species_relative_abundance)

top20_species_shotgun_tbl = tblsample_discovery %>% 
  filter(sampleid %in% tblshotgun_species_relative_abundance$sampleid) %>% 
  arrange(cluster_assignment) %>% 
  select(sampleid, cluster_assignment) %>%  
  left_join(tblshotgun_species_relative_abundance %>% 
              select(sampleid, all_of(col_choice_shotgun) ), by = "sampleid")

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

## plot Fig. 2f

ggplot(tblsample_discovery, aes(y = simpson_reciprocal, x = factor(cluster_assignment), fill = factor(cluster_assignment))) + 
  geom_boxplot() +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(name = "Simpson reciprocal", breaks=seq(0,60,by=20), limits = c(0,60)) +
  scale_x_discrete(name = "Cluster index", breaks=seq(1,10,by=1)) +
  scale_fill_manual(name = "Cluster", values = colors_asv_beta) +
  theme_classic() +
  geom_hline(yintercept = median(tblsample_discovery$simpson_reciprocal), linetype = 2) 

## plot Fig. 2g

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

ggplot(m2, aes(x = variable, y = value, fill = factor(ind))) + 
  geom_bar(stat = "identity") +
  ggtitle("Relative frequency k = 10") +
  scale_y_continuous("Percent") +
  scale_x_discrete("") +
  theme_classic() + 
  scale_fill_manual(values = colors_asv_beta) +
  theme(legend.position = "none", axis.text.x = element_text(angle=45,hjust=1))






