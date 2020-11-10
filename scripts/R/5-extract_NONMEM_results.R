#Data preparation: outcomes
#Input: results/NONMEM/*, dictionary.csv, index_PDXs_with_omics_data.Rdata
#Output: selected values for KG, KR, KD (outcomes_allData.Rdata)

#libraries
library(tidyr)
library(ggplot2)
library(readtext)
library(readr)

####Import data####
#Dictionary
dict <- read.csv("data/dictionary.csv", header = TRUE, sep = ",", quote = "", 
                 stringsAsFactors = FALSE)
load("data/index_PDX_omics_data.Rdata") #from script 4

#Functions
source("scripts/functions.R")


####Extract data from NONMEM####
dirs <- sort(list.dirs(path = "results/NONMEM/", full.names = FALSE, recursive = FALSE), decreasing = TRUE)
runName <- grep("KDKR", dirs, value = TRUE)[1]
folder_nonmem_result <- paste0("results/NONMEM/", runName)

allDataKD1 <- ExtractDataNONMEM(1, dict, folder = folder_nonmem_result, index = set_omics)
allDataKDKR1 <- ExtractDataNONMEM(2, dict, folder_nonmem_result, index = set_omics)
allDataKD2 <- ExtractDataNONMEM(3, dict, folder_nonmem_result, index = set_omics) #no additive error
allDataKDKR2 <- ExtractDataNONMEM(4, dict, folder_nonmem_result, index = set_omics) #no additive error

####Get metaresults per run####
runinfo_KD1 <- GetRunInfo(allDataKD1, 1, folder = folder_nonmem_result)
runinfo_KDKR1 <- GetRunInfo(allDataKDKR1, 2, folder = folder_nonmem_result)
runinfo_KD2 <- GetRunInfo(allDataKD2, 3, folder = folder_nonmem_result)
runinfo_KDKR2 <- GetRunInfo(allDataKDKR2, 4, folder = folder_nonmem_result)

####Combine data and convergence information####
allDataKD1 <- merge(allDataKD1, runinfo_KD1[, c("Treat_name", "converged", "model")])
allDataKDKR1 <- merge(allDataKDKR1, runinfo_KDKR1[, c("Treat_name", "converged", "model")])
allDataKD2 <- merge(allDataKD2, runinfo_KD2[, c("Treat_name", "converged", "model")])
allDataKDKR2 <- merge(allDataKDKR2, runinfo_KDKR2[, c("Treat_name", "converged", "model")])

####Extract OFV per individual####
KD1_OFV <- AddIndividualOFV(model = 1, folder = folder_nonmem_result)
KDKR1_OFV <- AddIndividualOFV(model = 2, folder = folder_nonmem_result)
KD2_OFV <- AddIndividualOFV(model = 3, folder = folder_nonmem_result)
KDKR2_OFV <- AddIndividualOFV(model = 4, folder = folder_nonmem_result)

all_OFV <- list(KD1_OFV = KD1_OFV, KDKR1_OFV = KDKR1_OFV, KD2_OFV = KD2_OFV, KDKR2_OFV = KDKR2_OFV)
all_info <- list(runinfo_KD1 = runinfo_KD1, runinfo_KDKR1 = runinfo_KDKR1, runinfo_KD2 = runinfo_KD2,
                 runinfo_KDKR2 = runinfo_KDKR2)
OFVs <- dict[, c("ID_name", "Treat_name", "ID_TREAT")]
for (i in 1:length(all_OFV)) {
  OFVs <- merge(OFVs, all_OFV[[i]], all = TRUE)
  OFVs <- merge(OFVs, all_info[[i]][, c("Treat_name", paste("converged", i, sep = "_"))], all.x = TRUE)
}
OFVs <- OFVs[apply(OFVs[, c(4,6,8,10)], 1, function(x) !all(is.na(x))), ]
OFVs <- OFVs[OFVs$ID_name %in% set_omics, ]
nr_obs <- as.data.frame(table(allDataKD1$ID_TREAT))


#Remove PDX with 4 or less observations
ind_small <- nr_obs$Var1[nr_obs$Freq <= 3]
OFVs <- OFVs[!OFVs$ID_TREAT %in% ind_small, ]

####Extract MSE(IPRED) per individual)####
library(tidyverse)
mseKD1 <- CalculateMSEIPRED(allDataKD1, model = 1)
mseKDKR1 <- CalculateMSEIPRED(allDataKDKR1, model = 2)
mseKD2 <- CalculateMSEIPRED(allDataKD2, model = 3)
mseKDKR2 <- CalculateMSEIPRED(allDataKDKR2, model = 4)

mseAll <- mseKD1 %>% left_join(mseKDKR1) %>% left_join(mseKD2) %>% left_join(mseKDKR2)
mseAll <- mseAll[!mseAll$ID_TREAT %in% ind_small, ]

#order mseAll the same as OFVs
mseAll <- mseAll[order(match(mseAll$ID_TREAT,OFVs$ID_TREAT)), ]
all(OFVs$ID_TREAT == mseAll$ID_TREAT)

####Calculate LRT, simpel model - complex model####
OFVs$LRT12 <- OFVs$OFV_1 - OFVs$OFV_2
OFVs$LRT34 <- OFVs$OFV_3 - OFVs$OFV_4
OFVs$LRT32 <- OFVs$OFV_3 - OFVs$OFV_2
OFVs$LRT14 <- OFVs$OFV_1 - OFVs$OFV_4


####Pick the usefull treatments####
convergence <- unique(OFVs[, c(1, 5, 7, 9, 11)])
convergence$pick <- NA
convergence$pick[rowSums(convergence[2:5], na.rm = T) == 0] <- "non"

ind_only <- rowSums(convergence[2:5], na.rm = T) == 1
convergence$pick[ind_only] <- apply(convergence[ind_only, 2:5], 1, function(x) which(x == 1))
convergence$pick[convergence$Treat_name == "untreated"] <- "0"
convergence$pick[rowSums(convergence[, 2:3], na.rm = T) == 2 & 
                   is.na(convergence$pick)] <- "LRT12"
convergence$pick[rowSums(convergence[, c(2, 5)], na.rm = T) == 2 & 
                   is.na(convergence$pick)] <- "LRT14"
convergence$pick[rowSums(convergence[, c(3, 4)], na.rm = T) == 2 & 
                   is.na(convergence$pick)] <- "LRT32"
convergence$pick[rowSums(convergence[, 4:5], na.rm = T) == 2 & 
                   is.na(convergence$pick)] <- "LRT34"
convergence$pick[convergence[, 2] == 1 & is.na(convergence$pick)] <- "1"


used_models <- data.frame()
for (i in 1:nrow(convergence)) {
  if (convergence$pick[i] == "non" | convergence$pick[i] == 0) {
    next
  }
  if (nchar(convergence$pick[i]) > 1) {
    used_models <- rbind(used_models, 
                         data.frame(Treat_name = convergence$Treat_name[i], 
                                    model = as.numeric(unlist(strsplit(substring(convergence$pick[i], 4, 5), "")))))
    
  } else {
    used_models <- rbind(used_models, 
                         data.frame(Treat_name = convergence$Treat_name[i], 
                                    model = convergence$pick[i]))
  }
}

shrinkage_nonmem <- merge(used_models, do.call(rbind, lapply(all_info, function(x) {
  x[, -ncol(x)]}))[, c("Treat_name", "model", "shrinkage_kd", "shrinkage_kr")], 
  all.x = TRUE) %>% 
  mutate(shrinkage_kr = ifelse(model %in% c(1, 3), NA, shrinkage_kr))
save(shrinkage_nonmem, file = "data/clean/shrinkage.Rdata")



####Collect chosen observations in allData####
OFVs$model <- character(nrow(OFVs))
rownames(OFVs) <- OFVs$ID_TREAT
for (treat in convergence$Treat_name) {
  ind_t <- which(OFVs$Treat_name == treat)
  column <- convergence$pick[convergence$Treat_name == treat]
  if (column %in% c("non", "2", "4")) {
    OFVs$model[ind_t] <- "non"
    print(paste("Treatment", treat, "has no converged models. Observations left out of final analysis"))
    next
  }
  if (column == "0") {
    OFVs$model[ind_t] <- "0"
    
    next
  }
  if (column %in% as.character(1:4)) {
    OFVs$model[ind_t] <- column
    next
  }
  LRs <- substr(column, start = 4, stop = 4)
  LRc <- substr(column, start = 5, stop = 5)
  LRT_t <- OFVs[ind_t, column]
  
  #Check the MSE's
  col_mse_s <- grep(paste0("IPRED", LRs), colnames(mseAll))
  col_mse_c <- grep(paste0("IPRED", LRc), colnames(mseAll))
  s_better <- mseAll[ind_t, col_mse_s] <= mseAll[ind_t, col_mse_c]
  
  #Check KR's
  col_kr_c <- grep(paste0("kr", LRc), colnames(mseAll))
  s_better <- (mseAll[ind_t, col_kr_c] > 1 | s_better)
  
  ind_crit <- LRT_t < qchisq(0.95, 1, lower.tail = T)
  OFVs$model[ind_t][(ind_crit | s_better)] <- LRs
  OFVs$model[ind_t][(!ind_crit & !s_better)] <- LRc
}

#Include untreated data
allDataKD1[allDataKD1$Treat_name == "untreated", c("kd", "kr", "model")] <- 0
#allDataKD1[allDataKD1$Treat_name == "untreated", c("model")] <- 0

allData <- rbind.data.frame(allDataKD1[allDataKD1$ID_TREAT %in% OFVs$ID_TREAT[OFVs$model %in% c("0", "1")], ], 
                            allDataKDKR1[allDataKDKR1$ID_TREAT %in% OFVs$ID_TREAT[OFVs$model == "2"], ],
                            allDataKD2[allDataKD2$ID_TREAT %in% OFVs$ID_TREAT[OFVs$model == "3"], ],
                            allDataKDKR2[allDataKDKR2$ID_TREAT %in% OFVs$ID_TREAT[OFVs$model == "4"], ])
allData$model <- NULL
allData$converged <- NULL
allData <- merge(OFVs[, c("ID_TREAT", "model")], allData, by = "ID_TREAT", all.y = T)

#Remove bad fitting treatment, later comment: WHY?
# allData <- allData[allData$Treat_name != "TAS266", ]

outcomes <- allData %>% 
  group_by(ID_TREAT) %>% 
  filter(TIME == min(TIME)) %>% 
  select(c("ID_TREAT", "kd", "kr", "ID_name", "Treat_name", "model", "kg", "DV")) %>% 
  mutate(base = DV, DV = NULL)

####Save outcome values####
save(outcomes, allData, file = "data/clean/outcomes_allData.Rdata")
save(outcomes, file = "data/clean/outcomes.Rdata")


####Visualize growth curves for visual check####
dirs <- sort(list.dirs(path = "results/NONMEM/", full.names = FALSE, recursive = FALSE), decreasing = TRUE)
runName <- grep("KG", dirs, value = TRUE)[1]

untreated <- read.table(paste0("results/NONMEM/", runName, "/m0_natural_growth.tab"), header = 1, skip = 1)
names(untreated)[11] <- "IPRED_NAT"

plotsChosen <- PlotGrowthCurves(allData, untreated, sorted = TRUE, KR = TRUE, model = TRUE)

####Save to files####
if (!dir.exists("results/figures")) {
  dir.create("results/figures")
}

pdf(file = "results/figures/modelled_growth_curves.pdf", width = 20, height = 16)
for (treat in names(plotsChosen)) {
  print(plotsChosen[[treat]])
}
dev.off()