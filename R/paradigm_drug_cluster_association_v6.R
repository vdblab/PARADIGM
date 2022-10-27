#' PARameters Associated with DynamIcs of Gut Microbiota 
#' 
#' Predict of environmental exposures on microbiome cluster dynamics and changes in taxonomic relative abundance. 
#' @param cluster_dynamic Training data of microbiome cluster dynamics. 
#' Must contain sequentially sequenced pairs of sample. 
#' Must have the exact same rows (subjects/samples) as environmental_exposure. 
#' @param environmental_exposure Environmental exposure data for model training. 
#' Must have the exact same rows (subjects/samples) as cluster_dynamic. 
#' Must contain numeric values. 
#' @param seed Set specific number as seed for reproducibility. 
#' @param numFold Number of cross-validation folds for elastic net logistic regression training 
#' @param numCluster Number of unique clusters in the cluster_dynamic data. 
#' @param euclidean_distance A matrix contains the euclidean distance between each pair of clusters. 
#' Must be of numCluster x numCluster dimension.
#' This parameter is only applicable for attractor transition (self-transition distance is 0).  
#' @param exposure_name Name of the environmental exposure to calculate response score. 
#' @param feature_profile_by_cluster A vector of clusters' mean value of microbiome feature of interest. 
#' Must have the same length as numCluster. 

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

#self transition for each cluster separately 
self_transition <- function(cluster_dynamic, 
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
    
    if (length(categorical_variable_id) > 0){
      
      x = dummy_cols(x, select_columns = categorical_variable_id,
                     remove_selected_columns = TRUE)
    } else{
      
      x = x
    }
    
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

attractor_transition <- function(cluster_dynamic, 
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
    
    if (length(categorical_variable_id) > 0){
      
      x = dummy_cols(x, select_columns = categorical_variable_id,
                     remove_selected_columns = TRUE)
    } else{
      
      x = x
    }
    
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

make_tm <- function(self_coef, 
                    attractor_coef, 
                    exposure_name, 
                    on_drug, 
                    euclidean_distance,
                    d = 0, 
                    numCluster = 10){
  
  drug_col = which(colnames(self_coef) == exposure_name)
  temp_matrix = matrix(nrow = numCluster,
                       ncol= numCluster)
  
  for (i in 1:10){
    w_self <- on_drug*self_coef[i,drug_col] + self_coef[i,1] + d*self_coef[i,2]
    
    p_self = exp(w_self)/(1 + exp(w_self))
    
    temp_matrix[i,i] <- p_self;
    
    for (j in setdiff(1:10,i)){
      
      distance <- euclidean_distance[i,j]
      
      w_attractor <- on_drug*attractor_coef[j,drug_col] + attractor_coef[j,1] +
        distance*attractor_coef[j, ncol(attractor_coef)] + d*attractor_coef[j,2] 
      
      p_attractor <- exp(w_attractor)/(1+exp(w_attractor))
      
      temp_matrix[i,j] <- p_attractor;
    }
  }
  temp_matrix = temp_matrix/rowSums(temp_matrix)
  
  return(temp_matrix)
}


bacteria_response_score <- function(self_coef, 
                                    attractor_coef, 
                                    exposure_name,
                                    feature_profile_by_cluster, 
                                    euclidean_distance, 
                                    d = 0, 
                                    numCluster = 10){
  
  tm_on = make_tm(self_coef, attractor_coef, exposure_name, on_drug = T, euclidean_distance, d)
  tm_off = make_tm(self_coef, attractor_coef, exposure_name, on_drug = F, euclidean_distance, d)
  
  score = 0
  for (i in 1:numCluster){
    for (j in 1:numCluster){
      c = (tm_on[i,j] - tm_off[i,j]) * log(feature_profile_by_cluster[j]/feature_profile_by_cluster[i])
      score = score + c
    }
  }
  
  return(score)
}

