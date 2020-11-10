#Data preperation: omics and outcomes
#Input: raw data (nm.3954-S2.xlsx)
#Output: index_PDXs_with_omics_data.Rdata

# libraries
library(readxl)

#Part I: Import omics data ######
# wide format RNA
pdx_rna <- read_excel("data/raw/nm.3954-S2.xlsx", sheet = "RNAseq_fpkm", n_max = 0)
# wide format copy numbers
pdx_cn <- read_excel("data/raw/nm.3954-S2.xlsx", sheet = "copy number", n_max = 0)

#Make data for filtering step:
set_rna <- names(pdx_rna)[-1]
set_cn <- names(pdx_cn)[-1]

# combine both sets (in this case set_cn is fully in set_rna
set_omics <- union(set_rna, set_cn)

save(set_omics, file = "data/index_PDX_omics_data.Rdata")
