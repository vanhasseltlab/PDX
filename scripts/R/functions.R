#Required packages
library(ggplot2)
library(readtext)
library(readr)

#Data preparation

#Script 5 read NONMEM files and visualize models
ExtractDataNONMEM <- function(model, dict, folder = getwd(), index = NULL) {
  allData <- NULL
  colnamesrun <- c("ID", "TREAT", "TIME", "DV", "DOSE", "CMT", "TID", "AMT", "EVID", "OCC", "IPRED",
                   "PRED", "kg", "kr" ,"kd" ,"RES" ,"WRES" ,"CWRES")
  for(rec1 in unique(dict$TREAT)){
    file_t <- paste0(folder, "/treat_", rec1, "/treat_", rec1,  "_", model, ".tab")
    if (!file.exists(file_t)) {
      print(paste("File", file_t, "does not exist"))
      next
    }
    output <- read.table(file_t, header = FALSE)
    names(output) <- colnamesrun
    output <- output[output$EVID == 0, ]
    allData <- rbind(allData, output)
  }
  allData$ID_TREAT <- paste(allData$ID, allData$TREAT, sep = "_")
  allData <- merge(allData, dict[, c("ID_TREAT", "ID_name", "Treat_name")], by = "ID_TREAT", 
                   all.x = TRUE)
  if (!is.null(index)) {
    allData <- allData[allData$ID_name %in% index, ]
  }
  return(allData)
}

GetRunInfo <- function(allData, model, folder) {
  run_info <- as.data.frame(matrix(NA, nrow = length(unique(allData$TREAT)), ncol = 15))
  colnames(run_info) <- c("model", "OFV", "Term_text", "Term_detail", "Cov_step_succes", "TH1", 
                          "TH2", "TH3", "OM1", "OM2", "OM3", "SI1", "SI2", "shrinkage_kd", "shrinkage_kr")
  run_info$TREAT <- unique(allData$TREAT)
  run_info$Treat_name <- unique(allData$Treat_name)
  
  #start
  for (i in 1:length(unique(allData$TREAT))) {
    treat <- unique(allData$TREAT)[i]
    file_t <- paste0(folder, "/treat_", treat, "/treat_", treat, "_", model, ".lst")
    recKD <- readtext(file_t, verbosity = 0)
    lines_text <- unlist(strsplit(recKD$text, "\n", fixed = TRUE))
    lines_text <- trimws(lines_text, which = "both")
    lines_text <- lines_text[lines_text != ""]
    run_info$OFV[i] <- parse_number(lines_text[grepl("#OBJV", lines_text)])
    run_info[i, c("TH1", "TH2", "TH3")] <- sapply(strsplit(lines_text[grep("TH 1", lines_text) + 1][1], split = "\\s+"), parse_number)
    omegas_ind <- grep("OMEGA - COV", lines_text)[1]
    run_info$OM1[i] <- parse_number(lines_text[omegas_ind + 3])
    run_info$OM2[i] <- parse_number(strsplit(lines_text[omegas_ind + 5], split = "\\s+")[[1]][3])
    run_info$OM3[i] <- parse_number(strsplit(lines_text[omegas_ind + 7], split = "\\s+")[[1]][4])
    sigmas_ind <- grep("SIGMA - COV", lines_text)[1]
    run_info$SI1[i] <- parse_number(lines_text[sigmas_ind + 3])
    run_info$SI2[i] <- parse_number(strsplit(lines_text[sigmas_ind + 5], split = "\\s+")[[1]][3])
    
    term_ind <- grep("0MINIMIZATION ", lines_text)
    run_info$Term_text[i] <- substring(lines_text[term_ind], 2)
    run_info$Term_detail[i] <- paste(lines_text[term_ind + (1:2)], collapse = " ")
    run_info$Cov_step_succes[i] <- ifelse(
      any(sapply(c("0S MATRIX", "0R MATRIX", "PSEUDO INVERSE OF", "EIGENVALUES NO.", "0MINIMIZATION TERMINATED"), 
                 function(x){grepl(x, lines_text)})), "NO", "YES")
    
    #find shrinkage
    run_info[i, c("shrinkage_kd", "shrinkage_kr")] <- parse_number(str_split(lines_text[grepl("ETASHRINKSD(%)", lines_text, fixed = TRUE)], "  ")[[1]][2:3])
  }
  run_info$model <- model
  run_info$converged <- as.numeric(run_info$Term_text == "MINIMIZATION SUCCESSFUL")
  run_info[, paste("converged", model, sep = "_")] <- run_info$converged
  return(run_info)
}

AddIndividualOFV <- function(allData = NULL, model, folder) {
  philes <- NULL
  for(rec_t in unique(dict$TREAT)){
    file_t <- paste0(folder, "/treat_", rec_t, "/treat_", rec_t,  "_", model, ".phi")
    if (!file.exists(file_t)) {
      print(paste("File", file_t, "does not exist"))
      next
    }
    phile <- read_table2(file_t, skip = 1, col_types = cols())
    phile <- cbind.data.frame(phile[, c("ID", "OBJ")], rec_t)
    names(phile) <- c("ID", "OFV", "TREAT")
    phile$ID_TREAT <- paste(phile$ID, phile$TREAT, sep = "_")
    if (any(is.na(phile$OFV))) {
      print("NA's in file")
    }
    philes <- rbind(philes, phile)
  }
  
  if (!is.null(allData)) {
    allData <- merge(allData, philes[, c("OFV", "ID_TREAT")], by = "ID_TREAT", all.x = TRUE)
    return(allData)
  }
  names(philes)[2] <- paste("OFV", model, sep = "_")
  return(philes[, c(paste("OFV", model, sep = "_"), "ID_TREAT")])
}

CalculateMSEIPRED <- function(allData, model = "") {
  mse <- allData %>% group_by(ID_TREAT) %>% summarize(MSE_IPRED = mean((IPRED - DV)^2), MSE_PRED = mean((PRED - DV)^2),
                                                      kr = kr[1])
  names(mse)[-1] <- paste0(names(mse)[-1], model)
  return(mse)
  
  #alternative retrieve KRs
  
}

PlotGrowthCurves <- function(allData, untreated, sorted = TRUE, KR = TRUE, model = FALSE, treatments = NULL) {
  if (!KR & !model) {
    allData$ID_kd <- paste0(allData$ID_name, "\nKD = ", allData$kd)
  }
  if (KR & !model) {
    allData$ID_kd <- paste0(allData$ID_name, "\nKD = ", allData$kd, "\nKR = ", allData$kr)
  }
  if (model) {
    allData$ID_kd <- paste0(allData$model, ", ", allData$ID_name, "\nKD = ", allData$kd, "\nKR = ", allData$kr)
  }
  allData$ID_kd <- as.factor(allData$ID_kd)
  
  if (sorted) {
    allData$ID_kd <- reorder(allData$ID_kd, -allData$kd)
  }
  predict_col <- "#F46E32"
  untreated_col <- "#5CB1EB"
  observation_col <- "#001158"
  
  if (is.null(treatments)) {
    treatments <- sort(unique(allData$Treat_name))
  }
  
  p_list <- list()
  for (treat in treatments) {
    dat1 <- allData[allData$Treat_name == treat, ]
    dat2 <- untreated[untreated$ID %in% dat1$ID, ]
    if (treat == "untreated") {
      dat1$ID_kd <- as.factor(dat1$ID_name)
      dat2 <- merge(dat2, unique(dat1[c("ID", "ID_kd")]), by = "ID")
      plot_col <- untreated_col
    } else {
      dat2 <- merge(dat2, unique(dat1[, c("ID", "ID_kd")]), by = "ID")
      plot_col <- predict_col
    }
    if (is.null(allData$converged)){
      title_plot <- paste0("Treatment ", treat)
    } else {
      title_plot <- paste0("Treatment ", treat, c(", NOT CONVERGED", ", CONVERGED")[dat1$converged[1] + 1])
    }
    p_list[[treat]] <- ggplot(dat1, aes(x = TIME, group = ID_kd)) +
      geom_point(aes(y = DV), colour = observation_col) +
      geom_line(data = dat2, aes(y = IPRED_NAT), colour = untreated_col, linetype = "dashed") +
      geom_line(aes(y = IPRED), colour = plot_col) +
      facet_wrap(~ ID_kd) +
      ggtitle(title_plot) +
      ylim(0, max(1200, max(dat1$DV))) +
      theme_bw() +
      theme(plot.title = element_text(size = 16, face = "bold"))
  }
  return(p_list)
}

#Script 6 data preparation for the lasso
CreateXandY <- function(pdx_, set_, filtercolumns = TRUE) {
  FilterColumn <- function(X) {
    X <- X[!duplicated(rownames(X)), ]
    
    uni <- apply(X, 2, function(x) length(unique(x)))
    ind_uni <- names(uni)[uni < 2]
    
    abz <- apply(X, 2, function(x) sum(x != 0))
    ind_nval <- names(abz)[abz < 2]
    
    # cor_Y <- apply(XY_$X[, -(1:54)], 2, function(x) cor(x, XY_$Y$kd))
    # cor_Y[is.na(cor_Y)] <- 0
    # ind_cor <- names(cor_Y)[abs(cor_Y) < 0.0025]
    
    ind_remove <- unique(c(ind_nval, ind_uni))
    return(ind_remove)
  }
  # extract matching outcomes
  outcomes_ <- outcomes[outcomes$ID_name %in% intersect(set_outcomes, set_), ]
  # create X an
  X_ <- t(as.matrix(pdx_[, -1], rownames.force = T))
  colnames(X_) <- pdx_$Sample
  X_ <- X_[as.character(outcomes_$ID_name), ]
  
  if (filtercolumns) {
    X_ <- X_[, !colnames(X_) %in% FilterColumn(X_)]
  }
  outcomes_$Treat_name <- as.factor(outcomes_$Treat_name)
  outcomes_$Treat_name <- relevel(outcomes_$Treat_name, ref = "untreated")
  X_meds <- model.matrix( ~ Treat_name, data = outcomes_)[, -1] #same order as outcome
  colnames(X_meds) <- gsub("Treat_name", "" , colnames(X_meds))
  # filter treatments in data
  X_meds <- X_meds[, colSums(X_meds) > 0]
  # combine medication and omics data
  X <- cbind(X_meds, X_)
  sets <- list(medicine = 1:ncol(X_meds), omics = (ncol(X_meds) + 1):ncol(X))
  return(list(X = X, Y = outcomes_, sets = sets))
}


# Data analysis

#Script 7a prediction
CalculateSABC <- function(model, intercept, truth, times) {
  sABC <- function(mod, truth) {
    mean(abs(mod - truth))/mean(truth)
  }
  
  df_sABC <- data.frame()
  for (ti in times) {
    ind_times <- model[, "time"] <= ti
    df_sABC <- rbind.data.frame(df_sABC, data.frame(sABC_intercept = sABC(intercept[ind_times, 2], truth[ind_times, 2]), 
                                                    sABC_model = sABC(model[ind_times, 2], truth[ind_times, 2]),
                                                    time = ti))
  }
  return(df_sABC)
}


#Script 7b variable selection
#Group lasso data preparation
SaveGroupLists <- function(filename) {
  name <- tail(unlist(strsplit(filename, "/|\\.")), 2)[1]
  txt_groups <- readLines(filename)
  group_list <- strsplit(txt_groups, "\t")
  names_groups <- lapply(group_list, function(x) x[1])
  #Structure grouplist accoording to package
  groups_cn <- groups_rna <- list()
  
  for (i in 1:length(group_list)) {
    group_list[[i]] <- group_list[[i]][group_list[[i]] != ""]
    
    group_list[[i]] <- group_list[[i]][-1]
    
    #Remove genes not in cn or rna
    groups_cn[[i]] <- intersect(group_list[[i]], cn_vars)
    groups_rna[[i]] <- intersect(group_list[[i]], rna_vars)
  }
  names(groups_cn) <- names(groups_rna) <- names_groups
  
  save(groups_cn, groups_rna, file = paste0("data/clean/groups/grlasso_", name, ".Rdata"))
}

#Group lasso
GroupLasso <- function(groups_omics, X_omic, outcomes, name_file = "") {
  #Which omics in pathways?
  omic_path <- unique(unlist(groups_omics))
  
  X <- as.matrix(X_omic[, omic_path])
  
  treatments <- unique(outcomes$Treat_name)
  treatments <- treatments[treatments != "untreated"]
  variables <- colnames(X)
  
  group_sizes <- lapply(groups_omics, length)
  
  beta <- data.frame()
  
  #Go through all treatments
  for (i in 1:length(treatments)) {
    
    #print progress
    print(paste0("Treatment ", i, " out of " , length(treatments), ": ",treatments[i]))
    
    dat_Y <- outcomes[outcomes$Treat_name == treatments[i], ]
    rownames(dat_Y) <- dat_Y$ID_name
    ids <- intersect(rownames(X), dat_Y$ID_name)
    n <- length(ids)
    for (outcome in c("kd", "kr")) {
      if (outcome == "kr") {
        if (var(dat_Y[ids, outcome]) == 0) next
        Y <- as.matrix(dat_Y[ids, outcome])
      } else if(outcome == "kd") {
        Y <- as.matrix(log(dat_Y[ids, outcome]))
      }
      rownames(Y) <- ids
      
      
      if (length(unique(Y)) > 3) {
        #Group lasso cross-validation
        print("Start CV")
        glasso_k_cv <- cv.grpregOverlap(X = X[ids, ], y = Y, group = groups_omics, 
                                        penalty = "grLasso", nfolds = 10)
        print("End CV")
        ind_min_lambda <- which(glasso_k_cv$lambda == glasso_k_cv$lambda.min)
        min_error <- glasso_k_cv$cve[ind_min_lambda]
        
        #Group lasso on min_lambda
        
        glasso_k_cv$fit$beta[, ind_min_lambda]
        
        ind_beta <- which(glasso_k_cv$fit$beta[-1, ind_min_lambda] != 0)
        
        #Save results in data frame 'beta'
        if(length(ind_beta) > 0) {
          print("extract selected groups")
          
          groups_select <- rowSums(glasso_k_cv$fit$incidence.mat[, ind_beta])
          
          #isolate beta's
          betas <- glasso_k_cv$fit$beta[ind_beta + 1, ind_min_lambda]
          betas <- data.frame(omic = names(betas), beta = betas)
          
          ind_selected_groups <- which(groups_select == group_sizes)
          selected_groups <- groups_omics[ind_selected_groups]
          
          #Save beta's
          for (j in names(selected_groups)) {
            beta <- rbind(beta, merge(data.frame(omic = selected_groups[[j]], group = j, 
                                                 outcome = outcome, treatment = treatments[i],
                                                 min_cve = min_error), 
                                      betas, all.x = T))
          }
        } #end any beta if
      } #end var(Y) if
    } #end outcome for loop
  } #end treatment for loop
  
  return(beta)
}