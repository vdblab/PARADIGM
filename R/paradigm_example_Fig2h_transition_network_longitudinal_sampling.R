library(tidyverse)

tblsample = read.csv("tblsample_cohort_master_table_deid_MSKCC_9674_Duke_500.csv")

##create a table with sample pairs to calculate transitions among subsequently collected samples 
cluster_on_t.end <- function(taxonomic_table = taxonomic_table){
  return(cbind.data.frame(oligos_id = taxonomic_table$oligos_id, subsequent_cluster = c(taxonomic_table$cluster_assignment[-1],NA)) )
}

day_on_t.end <- function(taxonomic_table){
  return(cbind.data.frame(oligos_id = taxonomic_table$oligos_id, day.y = c(taxonomic_table$day_relative_to_hct[-1],NA)) )
} 

subsequent_cluster = tblsample %>%
  arrange(PatientID) %>% 
  arrange(day_relative_to_hct) %>% 
  group_split(PatientID) %>%
  map_dfr(~cluster_on_t.end(.))

day.y = taxonomic_table %>%
  arrange(PatientID) %>% 
  arrange(day_relative_to_hct) %>% 
  group_split(PatientID) %>%
  map_dfr(~day_on_t.end(.))

pairing_samples = tblsample %>% 
  left_join(day.y, by = "oligos_id") %>% 
  left_join(subsequent_cluster, by = "oligos_id") %>% 
  filter(!is.na(subsequent_cluster)) %>% 
  mutate(dday = day.y - day_relative_to_hct) %>% 
  filter(dday <= 7) %>% 
  filter(set == "discovery")

transition_matrix = table(pairing_samples$cluster_assignment, pairing_samples$subsequent_cluster)
transition_matrix = transition_matrix/rowSums(transition_matrix)

##network map
library(igraph)

colors_asv_beta <- c("#A50026","#D73027","#F46D43","#C2A5CF","#3690C0",
                     "#FDAE61","#5AAE61","#B3E2CD","#762A83","#1B7837")

g <- graph.adjacency(transition_matrix, mode="directed", weighted=TRUE)
val2rgb <- colorRamp(rev(RColorBrewer::brewer.pal(11, "RdYlBu")[c(1:4,8:10)]), bias = 1)
# val2rgb <- colorRamp(rainbow(5)[1:6], bias = 2)
size = table(taxonomic_table$cluster_assignment[taxonomic_table$set == "discovery_set"])
edge.curve = c(0.15)
layout <- matrix(c(1,5, -0.5,2, 1.5,1.5, 0.55,0.3, 2,0.65, 3,1, 3,2, 3,3, 6,6, 7,7), nrow=10, ncol=2, byrow=T)

edge.colored = rep(NA, 100)
edge.colored[which(E(g)$weight < 0.1)] <- 'grey'
edge.colored[which(E(g)$weight >= 0.1)] <- rgb(val2rgb(E(g)$weight[which(E(g)$weight >= 0.1)]) / 255)

edge.widthed = rep(NA,100)
edge.widthed[which(edge.colored == "grey")] <- 0.2
edge.widthed[which(is.na(edge.widthed))] <- E(g)$weight[which(is.na(edge.widthed))] * 8

plot.igraph(
  x               = g,
  xlim            = c(-1.5, 1.5),
  ylim            = c(-1.5, 1.5),
  layout          = layout_in_circle,
  vertex.size     = sqrt(size), vertex.color = colors_asv_beta, edge.arrow.size = 0.25,
  edge.width      = edge.widthed,
  edge.color      = edge.colored,
  edge.curved     = edge.curve, 
  edge.loop.angle = scales::rescale(seq_along(E(g)), to=c(0, -2 * 3.14)) )

legend(
  x      = 1.5, 
  y      = 1, 
  fill   = c("grey", "#679DC9", "#FDAB60", "#DC3A2C"), 
  legend = c("<= 10%", "25%", "50%", "80%"), 
  title  = "Transition\nFrequency", 
  bty    = "n" )

##create a heatmap 

pheatmap((transition_matrix*100),cluster_rows= FALSE, cluster_cols = FALSE,
         display_numbers = round(transition_matrix,2), breaks = seq(0,100,by=1))

##calculate self transition frequency by cluster 
self_transition_vector = c()

for(i in 1:10){
  self_transition_vector = c(self_transition_vector, transition_matrix[i,i])
}

mean(self_transition_vector[8:10])
mean(self_transition_vector)
