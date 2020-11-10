##Grouping for group lasso

#package grpregOverlap
library(grpregOverlap)

# functions
source("scripts/functions.R")

#####list with indices or colnames of variables in group####
load("data/clean/lasso_data_filtered.Rdata")

#Get groups from txt files
SaveGroupLists("data/raw/WikiPathways_2019_Human.txt")
SaveGroupLists("data/raw/GO_Molecular_Function_2018.txt")
SaveGroupLists("data/raw/KEGG_2019_Human.txt")

cn_vars <- colnames(X_cn)
rna_vars <- colnames(X_rna)

####
####Apply group lasso####
start_time <- Sys.time()
set.seed(12345)
#Retrieve seed of CN Wiki trhough .Random.seed
load("data/clean/groups/grlasso_GO_Molecular_Function_2018.Rdata")
beta_cn_go <- GroupLasso(groups_cn, X_cn, outcomes, name_file = "CN_GO")
beta_rna_go <- GroupLasso(groups_rna, X_rna, outcomes, name_file = "RNA_GO")

load("data/clean/groups/grlasso_KEGG_2019_Human.Rdata")
beta_cn_kegg <- GroupLasso(groups_cn, X_cn, outcomes, name_file = "CN_KEGG")
beta_rna_kegg <- GroupLasso(groups_rna, X_rna, outcomes, name_file = "RNA_KEGG")


load("data/clean/groups/grlasso_WikiPathways_2019_Human.Rdata")
beta_cn_wiki <- GroupLasso(groups_cn, X_cn, outcomes, name_file = "CN_Wiki")
beta_rna_wiki <- GroupLasso(groups_rna, X_rna, outcomes, name_file = "RNA_Wiki")


save(beta_cn_go, beta_rna_go, beta_cn_kegg, beta_rna_kegg, beta_cn_wiki, beta_rna_wiki, file = "results/beta_grplasso.Rdata")
print(start_time)
print(Sys.time())
print(Sys.time() - start_time)

