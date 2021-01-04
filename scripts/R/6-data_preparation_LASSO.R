#Data preparation: X and Y, filter variables
#Input: outcome data (data/clean/outcomes_allData.Rdata from script 5), raw omic data (data/raw/nm.3954-S2.xlsx)
#Output: data ready for LASSO, in X and y format (data/clean/lasso_data_filtered.Rdata)

# libraries
library(readxl)

# functions
source("scripts/R/functions.R")

#Import data
load("data/clean/outcomes_allData.Rdata")
# wide format RNA
pdx_rna <- read_excel("data/raw/nm.3954-S2.xlsx", sheet = "RNAseq_fpkm")
# wide format copy numbers
pdx_cn <- read_excel("data/raw/nm.3954-S2.xlsx", sheet = "copy number")

#Sets with indices
set_outcomes <- as.character(unique(outcomes$ID_name))
set_rna <- names(pdx_rna)[-1]
set_cn <- names(pdx_cn)[-1]

# keep only omics with outcomes
rownames(pdx_rna) <- pdx_rna$Sample
pdx_rna <- pdx_rna[, c("Sample", intersect(set_outcomes, set_rna))]

#ATTENTION:cn contains multiple measurements of same variable: 
# now put all in model with different names
while (length(pdx_cn$Sample[duplicated(pdx_cn$Sample)]) > 0) {
  pdx_cn$Sample[duplicated(pdx_cn$Sample)] <- paste(pdx_cn$Sample[duplicated(pdx_cn$Sample)],
                                                    "i", sep = "")
  print("i")
}
#Remove FocalCNScore and ArmLevelCNScore from cn
pdx_cn <- pdx_cn[!pdx_cn$Sample %in% c("FocalCNScore", "ArmLevelCNScore"), ]

rownames(pdx_cn) <- pdx_cn$Sample
pdx_cn <- pdx_cn[, c("Sample", intersect(set_outcomes, set_cn))]

#Part III: Create X and y for the omics sets and medication (57 medicines and untreated as baseline) #######
XY_rna <- CreateXandY(pdx_rna, set_rna)
XY_cn <- CreateXandY(pdx_cn, set_cn)
save(XY_rna, XY_cn, file = "data/clean/XY_data_filtered.Rdata")

#For LASSO per treatment
X_cn <- t(pdx_cn[, -1])
colnames(X_cn) <- pdx_cn$Sample
X_rna <- t(pdx_rna[, -1])
colnames(X_rna) <- pdx_rna$Sample

save(X_cn, X_rna, outcomes, file = "data/clean/lasso_data_filtered.Rdata")





