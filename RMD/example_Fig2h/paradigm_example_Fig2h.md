PARADIGM example: Figure 2h- cluster transition network
================
Chi L. Nguyen
10/17/2022

Cohort characteristics are outlined in ‘tblsample’ table.

``` r
tblsample = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tblsample_cohort_master_table_deid_MSKCC_9674_Duke_500_101222.csv")
```

We created sample-pair-wise table for patients with longitudinal
sampling in the MSKCC discovery cohort.

``` r
cluster_on_t_end <- function(taxonomic_table){
  return(cbind.data.frame(oligos_id = taxonomic_table$oligos_id, subsequent_cluster = c(taxonomic_table$cluster_assignment[-1],NA), 
                          day.y = c(taxonomic_table$day_relative_to_hct[-1],NA)) )
}

subsequent_cluster = tblsample %>%
  filter(set == "discovery") %>% 
  arrange(PatientID) %>% 
  arrange(day_relative_to_hct) %>% 
  group_split(PatientID) %>%
  map_dfr(~cluster_on_t_end(.))
head(subsequent_cluster)
```

    ##        oligos_id subsequent_cluster day.y
    ## 1 1000A..pool483                  4    -4
    ## 2 1000B..pool483                  5     6
    ## 3 1000C..pool483                  5     9
    ## 4 1000D..pool483                  5    13
    ## 5 1000E..pool483                 NA    NA
    ## 6 1002B..pool483                  5     3

``` r
pairing_samples = tblsample %>% 
  left_join(subsequent_cluster, by = "oligos_id") %>% 
  filter(!is.na(subsequent_cluster)) %>% 
  mutate(dday = day.y - day_relative_to_hct) %>% 
  filter(dday <= 7) %>% 
  filter(set == "discovery")
head(pairing_samples)
```

    ##   sampleid      oligos_id PatientID cluster_assignment day_relative_to_hct
    ## 1    1000A 1000A..pool483      1000                  2                  -9
    ## 2    1000C 1000C..pool483      1000                  5                   6
    ## 3    1000D 1000D..pool483      1000                  5                   9
    ## 4    1002B 1002B..pool483      1002                  1                  -1
    ## 5    1002C 1002C..pool483      1002                  5                   3
    ## 6    1004A 1004A..pool483  FMT.0002                  1                 -10
    ##         set simpson_reciprocal     Blautia Enterococcus Erysipelatoclostridium
    ## 1 discovery           13.38430 0.406969925  0.000000000            0.005193935
    ## 2 discovery            2.08790 0.000000000  0.001073784            0.000000000
    ## 3 discovery            2.39139 0.001466583  0.001780851            0.000000000
    ## 4 discovery            5.38780 0.112959720  0.000000000            0.071706558
    ## 5 discovery            5.15961 0.005445227  0.156630365            0.000000000
    ## 6 discovery            9.71336 0.114697162  0.000000000            0.061245235
    ##   institution     tsne1      tsne2 shotgun_available tsne1_shotgun
    ## 1       MSKCC -0.746297 -0.0408929                No            NA
    ## 2       MSKCC  0.207460  0.4915400                No            NA
    ## 3       MSKCC  0.244664  0.2800610               Yes    -0.7709889
    ## 4       MSKCC -0.475266 -0.1420040                No            NA
    ## 5       MSKCC  0.359112  0.3055980                No            NA
    ## 6       MSKCC -0.516421 -0.3004220                No            NA
    ##   tsne2_shotgun subsequent_cluster day.y dday
    ## 1            NA                  4    -4    5
    ## 2            NA                  5     9    3
    ## 3     -0.351067                  5    13    4
    ## 4            NA                  5     3    4
    ## 5            NA                  7     7    4
    ## 6            NA                  3    -3    7

We then calculated the cluster transition frequency among subsequently
collected samples.

``` r
transition_matrix = table(pairing_samples$cluster_assignment, pairing_samples$subsequent_cluster)
transition_matrix = transition_matrix/rowSums(transition_matrix)
transition_matrix
```

    ##     
    ##                1           2           3           4           5           6
    ##   1  0.408256881 0.203363914 0.071865443 0.062691131 0.021406728 0.068807339
    ##   2  0.120393120 0.474201474 0.126535627 0.058968059 0.011056511 0.052825553
    ##   3  0.060263653 0.120527307 0.506591337 0.043314501 0.011299435 0.124293785
    ##   4  0.052953157 0.026476578 0.018329939 0.384928717 0.073319756 0.054989817
    ##   5  0.015479876 0.021671827 0.012383901 0.095975232 0.603715170 0.024767802
    ##   6  0.046709130 0.014861996 0.116772824 0.065817410 0.031847134 0.554140127
    ##   7  0.021686747 0.024096386 0.019277108 0.043373494 0.033734940 0.036144578
    ##   8  0.039001560 0.009360374 0.007800312 0.101404056 0.076443058 0.039001560
    ##   9  0.012552301 0.014644351 0.000000000 0.156903766 0.037656904 0.002092050
    ##   10 0.005945303 0.009512485 0.001189061 0.007134364 0.010701546 0.004756243
    ##     
    ##                7           8           9          10
    ##   1  0.022935780 0.081039755 0.032110092 0.027522936
    ##   2  0.033169533 0.047911548 0.035626536 0.039312039
    ##   3  0.039548023 0.048964218 0.018832392 0.026365348
    ##   4  0.030549898 0.158859470 0.177189409 0.022403259
    ##   5  0.043343653 0.102167183 0.043343653 0.037151703
    ##   6  0.042462845 0.065817410 0.036093418 0.025477707
    ##   7  0.474698795 0.086746988 0.024096386 0.236144578
    ##   8  0.063962559 0.536661466 0.079563183 0.046801872
    ##   9  0.023012552 0.098326360 0.606694561 0.048117155
    ##   10 0.120095125 0.023781213 0.010701546 0.806183115

We could visualize cluster transitions in a network map.

``` r
colors_asv_beta <- c("#A50026","#D73027","#F46D43","#C2A5CF","#3690C0",
                     "#FDAE61","#5AAE61","#B3E2CD","#762A83","#1B7837")

g <- graph.adjacency(transition_matrix, mode="directed", weighted=TRUE)
val2rgb <- colorRamp(rev(RColorBrewer::brewer.pal(11, "RdYlBu")[c(1:4,8:10)]), bias = 1)
size = table(tblsample$cluster_assignment[tblsample$set == "discovery"])
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
```

![](paradigm_example_Fig2h_files/figure-gfm/Fig2h_network_map-1.png)<!-- -->

We could also visualize cluster transitions in a 10x10 heatmap.

``` r
pheatmap((transition_matrix*100),cluster_rows= FALSE, cluster_cols = FALSE,
         display_numbers = round(transition_matrix,2), breaks = seq(0,100,by=1))
```

![](paradigm_example_Fig2h_files/figure-gfm/Fig2h_heat_map-1.png)<!-- -->
