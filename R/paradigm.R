library(tidyverse)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(R.utils)
library(caret)
library(e1071)
library(pheatmap)
library(glmnet)
library(randomcoloR)
library(fastDummies)

tbldaily = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tbldaily_sampling_PARADIGM_input_101122.csv")
tblpatient = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tblpatient_cohort_characteristics_master_table_deid_MSKCC_1227_Duke_142.csv")
tbleuclidean_distance = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tbleuclidean_distance_10clusters.csv")
tblcounts = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tblcounts_master_table_deid_MSKCC_9674_Duke_500_101222.csv")
tblsample = read.csv("~/Desktop/Backup from server /backup_Oct2020/First rotation - VDB/Data - ASV/Sept7_rebuttal/deposited dataset/tblsample_cohort_master_table_deid_MSKCC_9674_Duke_500_101222.csv")

## create two data frames, one containing the response variable (cluster self/attractor transtions)
## and one containing the predictors (drug exposures, time of sample collection). 

tbldaily_clinical = tbldaily %>% left_join(tblpatient %>% select(PatientID = PatientID, 
                                                        intensity, 
                                                        simplesource, 
                                                        disease_simple), by = "PatientID")
cluster_dynamic = tbldaily_clinical[,1:5]
environmental_exposures = tbldaily_clinical[,c(6,9:ncol(tbldaily_clinical))]

self_transition_categorical_variable <- function(cluster_dynamic,  
                                                 environmental_exposures, 
                                                 seed = 2208, 
                                                 numFold = 10, 
                                                 numCluster = 10){
                           
  
  self_coef <- data.frame()
  
  for (n in 1:numCluster){
    
    print(paste("Processing cluster", n, sep = " ") )
    
    ##select rows to include in the training data 
    
    row_ind = which(cluster_dynamic$n10.x == n)
    
    cluster_data <- cluster_dynamic[row_ind, ]
    self <- cluster_data$n10.y == n
    y = as.matrix(self)
    x = environmental_exposures[row_ind, ]
    
    categorical_variable_id = colnames(x)[which(unlist(lapply(x, class)) %in% c("character", "factor"))]
    
    x = dummy_cols(x, select_columns = categorical_variable_id,
                   remove_selected_columns = TRUE)
    
    ##create pre-specified folds based on unique PatientID 
    
    set.seed(seed)
    uniquePatientID = unique(cluster_data$PatientID)
    foldlist = groupKFold(cluster_data$PatientID, k = numFold)
    
    fit <- train(x, as.factor(as.numeric(y)), method = "glmnet",
                 trControl = trainControl("cv", number = numFold, index = foldlist), 
                 penalty.factor = c(1, rep(1, (ncol(x)-1))), 
                 tuneLength = 10)
    
    self_coef<- rbind.data.frame(self_coef, 
                                 as.vector(coef(fit$finalModel, fit$bestTune$lambda)) )
    
  }
  
  colnames(self_coef) = c("Intercept", colnames(x))
  return(self_coef)
}

attractor_transition_categorical_variable <- function(cluster_dynamic, 
                                                      environmental_exposures,
                                                      euclidean_distance, 
                                                      seed = 2208, 
                                                      numFold = 10, 
                                                      numCluster = 10){
                                                    
                                 
  
  attractor_coef <- data.frame()
  
  for (n in 1:numCluster){
    
    print(paste("Processing cluster", n, sep = " ") )
    
    ##select rows to include in the training data 
    
    row_ind = which(cluster_dynamic$n10.x != n)
    
    cluster_data <- cluster_dynamic[row_ind, ]
    nonself <- cluster_data$n10.y == n
    euclidean_distance_vector <- euclidean_distance[cbind(cluster_data$n10.x, n)]
    y = as.matrix(nonself)
    x = environmental_exposures[row_ind, ]
    
    categorical_variable_id = colnames(x)[which(unlist(lapply(x, class)) %in% c("character", "factor"))]
    
    x = dummy_cols(x, select_columns = categorical_variable_id,
                   remove_selected_columns = TRUE)
    
    x = cbind(x, euclidean_distance_vector = euclidean_distance_vector)
    
    ##create pre-specified folds based on unique PatientID 
    
    set.seed(seed)
    uniquePatientID = unique(cluster_data$PatientID)
    foldlist = groupKFold(cluster_data$PatientID, k = numFold)
    
    fit <- train(x, as.factor(as.numeric(y)), method = "glmnet",
                 trControl = trainControl("cv", number = numFold, index = foldlist), 
                 penalty.factor = c(1, rep(1, (ncol(x)-1))),
                 tuneLength = 10)
    
    attractor_coef <- rbind.data.frame(attractor_coef, 
                                       as.vector(coef(fit$finalModel, fit$bestTune$lambda)) )
  }
  
  colnames(attractor_coef) = c("Intercept", colnames(x))
  return(attractor_coef)
}

self_coefficient_matrix_with_categorical = self_transition_categorical_variable(cluster_dynamic = cluster_dynamic,
                                                                                environmental_exposures = environmental_exposures, 
                                                                                seed = 2208, 
                                                                                numFold = 10, 
                                                                                numCluster = 10)
                                     

attractor_coefficient_matrix_with_categorical = attractor_transition_categorical_variable(cluster_dynamic = cluster_dynamic, environmental_exposures = environmental_exposures,
                                                                                          euclidean_distance = tbleuclidean_distance, 
                                                                                          seed = 2208, 
                                                                                          numFold = 10, 
                                                                                          numCluster = 10)
                                                    
