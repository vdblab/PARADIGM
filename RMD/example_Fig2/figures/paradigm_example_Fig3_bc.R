library(pheatmap)

#' Source script runs elastic net regularized regression to learn the coefficients of each drug exposure to cluster self
#' and attractor transitions, then convert cluster-drug associations into drug-taxon response scores. 

source("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Script/github_deposit/paradigm_drug_cluster_associations_v4_preSpecified_PatientFold.R", echo=F)

#' Table details: 
#' tbldaily.csv Training data of microbiome cluster dynamics.
#' tblpatient.csv Patient characteristics.
#' tbleuclidean_distance.csv A 10x10 matrix of pair-wise cluster euclidean distance.
#' tblcounts.csv Counts and taxonomic classifications of ASVs.
#' tblsample.csv Sample characteristics. 

tbldaily = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tbldaily_sampling_PARADIGM_input_101122.csv")
tblpatient = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tblpatient_cohort_characteristics_master_table_deid_MSKCC_1227_Duke_142.csv")
tbleuclidean_distance = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tbleuclidean_distance_10clusters.csv")
tblcounts = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tblcounts_master_table_deid_MSKCC_9674_Duke_500_101222.csv")
tblsample = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tblsample_cohort_master_table_deid_MSKCC_9674_Duke_500_101222.csv")

## create two data frames, one containing the response variable (cluster self/attractor transtions)
## and one containing the predictors (drug exposures, time of sample collection). 

cluster_dynamic = tbldaily[,1:5]
environmental_exposures = tbldaily[,c(6,9:ncol(tbldaily))]

## calculate coefficients of predictors to cluster self and attractor transitions.
## output is a matrix, with predictors column-wise and coefficient values row-wise. 

self_coefficient_matrix = self_transition(cluster_dynamic = cluster_dynamic, 
                                          environmental_exposures = environmental_exposures, 
                                          seed = 2208, 
                                          numFold = 10, 
                                          numCluster = 10)

attractor_coefficient_matrix = attractor_transition(cluster_dynamic = cluster_dynamic, 
                                                    environmental_exposures = environmental_exposures,
                                                    euclidean_distance = tbleuclidean_distance, 
                                                    seed = 2208, 
                                                    numFold = 10, 
                                                    numCluster = 10)

## plot Fig. 3b

self_coefficient_to_plot = t(self_coefficient_matrix)
attractor_coefficient_to_plot = t(attractor_coefficient_matrix)
tblgraph_drug_exposure_classification = read.csv("tblgraph_drug_exposure_classification.csv") 
drug_group_annotation = tblgraph_drug_exposure_classification %>% select(group, ind)
rownames(drug_group_annotation) = tblgraph_drug_exposure_classification$exposure_name

set = c("#FFFFFF", RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(12, "Set3")[3:7])
set = list(group = set)
names(set$group) = unique(tblgraph_drug_exposure_classification$group)

ind_df = match(rownames(drug_group_annotation), rownames(self_coefficient_to_plot)) 

out_self = pheatmap(self_coefficient_to_plot[c(1,2,ind_df),], col = colorRampPalette(c("navy", "white", "firebrick3"), bias = 1.3)(150), 
                    annotation_row = drug_group_annotation, annotation_colors = set, 
                    cluster_cols = F, cluster_rows = F, breaks = seq(-3,4,by=0.05), main = "Self Transition")

out_attractor = pheatmap(attractor_coefficient_to_plot[c(1,2,ind_df),], col = colorRampPalette(c("navy", "white", "firebrick3"), bias = 1.3)(150), 
                         annotation_row = drug_group_annotation, annotation_colors = set,
                         cluster_cols = F, cluster_rows = F, breaks = seq(-3,4,by=0.05), main = "Attractor Transition")

## calculate drug-taxon response score using the coefficient matrices of cluster self and attractor transitions. 

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

## reiterate the script across drug exposures of interest and microbiome features of interest. 

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

## plot Fig.3c

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
  scale_fill_gradientn(colours = colorRampPalette(c("navy", "white", "firebrick3"), bias = 1)(30)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.y = element_blank(), axis.title.x = element_blank())

## plot Fig.3d

ggplot(data = response_score_feature_aggregate, 
       aes(x=Enterococcus, y=simpson_reciprocal)) +
  geom_point() +
  geom_smooth(method = "lm", linetype = "dashed") +
  stat_cor() +
  theme_classic() +
  xlab("Enterococcus \n response scores") +
  ylab("Alpha-diversity \n response scores")
  

  
