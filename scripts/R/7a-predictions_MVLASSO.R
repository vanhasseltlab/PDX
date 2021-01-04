#Prediction of KG, KD, KR from omics data through monte carlo cross-validation

##Prediction with multple lasso (mgaussian)
###Prediction of growth curves

#libraries
library(abind)
library(RxODE)
library(gridExtra)
library(glmnet)
library(tidyverse)

# functions
source("scripts/R/functions.R")

#import data
load("data/clean/lasso_data_filtered.Rdata")
load("data/clean/outcomes_allData.Rdata")

#create result folder
if (!dir.exists("results/MVLasso_predictions/")) {
  dir.create("results/MVLasso_predictions/")
}

#Choose model
include_KR <- FALSE # TRUE or FALSE

####Using LOO CV to predict TGcurves.
treatments <- unique(outcomes$Treat_name)
treatments <- treatments[treatments != "untreated"]
X_ <- X_cn
prefix_output <- paste(c("KGKD", "KGKDKR")[include_KR + 1], "- ")

####Perform 10-fold CV and extract predictions per lambda. Retrieve minimizing errors.####
#size of Monte Carlo
M <- 20

#save elements of monte carlo CV
smooth_curves <- list()
predictions <- list()
lambda_mat <- NULL

for (i in 1:length(treatments)) {
  dat_Y <- outcomes[outcomes$Treat_name == treatments[i], ]
  rownames(dat_Y) <- dat_Y$ID_name
  ids <- intersect(rownames(X_), dat_Y$ID_name)
  n <- length(ids)
  outcome <- c("kg", "kd", "kr")
  
  #Extract estimated KG, KD and KR
  
  #Model with or without KR
  if (var(dat_Y[ids, "kr"]) > 0 & length(unique(unlist(dat_Y[ids, "kr"]))) > 2 & include_KR) {
    Y <- as.matrix(cbind(log(dat_Y[ids, outcome[1:2]]), kr = dat_Y[ids, outcome[3]]))
  } else {
    outcome <- outcome[1:2]
    Y <- as.matrix(log(dat_Y[ids, outcome]))
  }
  colnames(Y)[1:2] <- paste0("log", outcome[1:2])
  rownames(Y) <- ids
  
  #Extract same ID's from X
  X <- as.matrix(X_[ids, ])
  
  #check for predictiveness - select lambda
  
  #gather different lambda's
  lambda_i <- as.data.frame(matrix(NA, ncol = 6, nrow = M))
  cv_pred <- array(0, dim = c(nrow(Y), ncol(Y), M))
  rownames(cv_pred) <- ids
  missing_calcs <- 0
  smooth_curve <- NULL
  
  cat(paste("\nTreatment", treatments[i], "\n"))
  for (j in 1:M) {
    #print progress
    cat('\r', "Progression Monte Carlo at", round(j/M*100), "%")
    
    set.seed(123 + j)
    
    #Multivariate outcome option "MV"
    #10-fold CV mvnorm with error handling
    glm_cv <- tryCatch(cv.glmnet(X, Y, family = "mgaussian", alpha = 1, nfolds = 10, standardize.response = T),
                       error = function(e) FALSE)
    
    #save results output
    if (class(glm_cv) != "cv.glmnet") {
      print(paste("Run", j, "for treatment", treatments[i], "did not run."))
      cv_pred[, , j] <- NA
      lambda_i[j, ] <- NA
      missing_calcs <- missing_calcs + 1
      next
    }
    cv_pred[, , j] <- predict(glm_cv, X, s = "lambda.min")[, , 1]
    lambda_i[j, ] <- data.frame(V1 = glm_cv$lambda.min, V2 = which(glm_cv$lambda == glm_cv$lambda.min),
                                V3 = glm_cv$lambda.1se, V4 = which(glm_cv$lambda == glm_cv$lambda.1se),
                                V5 = j, V6 = paste0(outcome, collapse = ""))
    if (is.null(smooth_curve)) {
      smooth_curve <- cbind(glm_cv$lambda, glm_cv$cvm, glm_cv$cvsd, glm_cv$nzero)
    } else {
      smooth_curve[, 2:3] <- smooth_curve[, 2:3] + cbind(glm_cv$cvm, glm_cv$cvsd)
    }
  }
  
  colnames(smooth_curve) <- c("lambda", "cverror", "cverror_sd", "nonzero")
  smooth_curve[, 2:4] <- smooth_curve[, 2:4]/(M - missing_calcs)
  lambda_mat <- rbind(lambda_mat, lambda_i)
  
  #bring predictions back to actual size
  cv_pred[, 1:2, ] <-  exp(cv_pred[, 1:2, ])
  colnames(cv_pred) <- outcome
  
  
  smooth_curves[[i]] <- smooth_curve
  predictions[[i]] <- cv_pred
}

ind_treat <- which(lambda_mat$replication == M)[which(lambda_mat$replication == M)%%M == 0]
rep_treat <- c(ind_treat[1], abs(diff(ind_treat)))
lambda_mat <- data.frame(lambda_mat, treatment = rep(treatments, times = rep_treat))
colnames(lambda_mat)[1:6] <- c("lambda_min", "ind_lambda_min", "lambda_1se", "ind_lambda_1se", "replication", "parameter")
names(predictions) <- treatments
names(smooth_curves) <- treatments

save(lambda_mat, file = paste0("results/MVLasso_predictions/", prefix_output, "lambda_repeatedCV.Rdata"))
save(predictions, smooth_curves, file = paste0("results/MVLasso_predictions/", prefix_output, "predictions_repeatedCV.Rdata"))

#Extract sABCs
load(file = paste0("results/MVLasso_predictions/", prefix_output, "lambda_repeatedCV.Rdata"))
load(file = paste0("results/MVLasso_predictions/", prefix_output, "predictions_repeatedCV.Rdata"))

#Retrieve intercept model either from file or predict intercept in LOOCV ####
if (file.exists("results/MVLasso_predictions/intercept_model.Rdata")) {
  load("results/MVLasso_predictions/intercept_model.Rdata")
} else {
  #Calculate parameters from intercept/null model
  dat_intercept <- outcomes[, c("Treat_name", "ID_TREAT", "kg", "kd", "kr", "base", "ID_name")] %>% 
    mutate(kg_hat = 0, kd_hat = 0, kr_hat = 0, area_scaled100 = 0, area_scaled56 = 0)
  
  for (treat in unique(dat_intercept$Treat_name)) {
    dat <- dat_intercept[dat_intercept$Treat_name == treat, ]
    
    for (i in 1:nrow(dat)) {
      dat[i, ]$kg_hat <- exp(mean(log(dat[-i, ]$kg)))
      dat[i, ]$kd_hat <- exp(mean(log(dat[-i, ]$kd)))
      dat[i, ]$kr_hat <- median(dat[-i, ]$kr)
      
    }
    dat_intercept[dat_intercept$Treat_name == treat, ] <- dat
  }
  save(dat_intercept, file = "results/MVLasso_predictions/intercept_model.Rdata")
}

mod1 <- RxODE({
  d/dt(V) <- kg*V - kd*exp(-kr*t)*V
})
days <- seq(0, 100, by = 0.5)
ev <- eventTable(amount.units = 'mm^3', time.units = 'days')%>%
  add.sampling(days)

sABC_treat_total <- data.frame()

#loop through treatments
for (treat in treatments) {
  cv_pred <- predictions[[treat]]
  if (ncol(cv_pred) == 2) {
    cv_pred <- abind(cv_pred, array(0, replace(dim(cv_pred), 2, 1)), along = 2)
    colnames(cv_pred)[3] <- "kr"
  }
  cv_pred[is.na(cv_pred)] <- 0
  
  lambdas <- lambda_mat[lambda_mat[, "treatment"] == treat, ]
  M_set <- which(!duplicated(lambdas[, "lambda_min"]))

  sABC_treat <- data.frame()
  #loop through individuals
  for (id in rownames(cv_pred)) {
    
    pred_id <- cv_pred[id, ,]
    intercept_id <- dat_intercept[dat_intercept$ID_name == id & dat_intercept$Treat_name == treat, 
                                  c("kg_hat", "kd_hat", "kr_hat")]
    names(intercept_id) <- c("kg", "kd", "kr")
    
    outcome_id <- outcomes[outcomes$ID_name == id & outcomes$Treat_name == treat, ]
    y_base <- outcome_id$base
    
    out_intercept <- mod1$solve(intercept_id, ev, y_base)
    colnames(out_intercept)[2] <- "intercept"
    out_truth <- mod1$solve(outcome_id[, c("kg", "kd", "kr")], ev, y_base)
    colnames(out_truth)[2] <- "truth"
    
    curves <- data.frame()
    average_curve <- rep(0, length(days))
    weights <- table(pred_id[1, ])/M
    
    for (j in 1:length(M_set)){
      m <- M_set[j]
      input_pars_pred <- c(pred_id[1, m], pred_id[2, m], pred_id[3, m])
      out_pred <- mod1$solve(input_pars_pred, ev, y_base)
      colnames(out_pred) <- c("time", "predicted")
      curves <- rbind.data.frame(curves, data.frame(out_pred, m = m))
      
      average_curve <- average_curve + out_pred[, "predicted"]*weights[j]
    }
    
    average_curve <- data.frame(time = out_pred[, "time"], volume = average_curve)
    
    times_interest <- c(7, 14, 28, 56, 75, 100)
    
    sABC_treat <- rbind.data.frame(sABC_treat, data.frame(CalculateSABC(average_curve, out_intercept, out_truth, times_interest), 
                                                          ID_name = id, Treat_name = treat))
    
  }
  sABC_treat_total <- rbind.data.frame(sABC_treat_total, sABC_treat)
}

save(sABC_treat_total, file = paste0("results/MVLasso_predictions/", prefix_output, "sABC_over_time.Rdata"))