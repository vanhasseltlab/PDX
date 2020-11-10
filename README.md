# PDX data analysis

Scripts and pipeline for data analsysis in corresponding manuscript: "Identification of high-dimensional omics-derived predictors for tumor growth dynamics using machine learning and pharmacometric modeling"


### Scripts
Data preparation and data analsysis should be conducted in order of script names. Scripts 7a and 7b are independent of each other and can be run in any order. Scripts 2 and 3 require NONMEM installed on the machine.

1. data_preparation_NONMEM
2. TGI_natural_growth_NONMEM
3. TGI_treatments_NONMEM
4. data_preparation_omics
5. extract_NONMEM_results
6. data_preparation_LASSO
7. 
    1. predictions_MVLASSO
    2. pathway_selection_groupLASSO
  

### Data

* PDX study data 'nm.3954-S2.xlsx' from https://www.nature.com/articles/nm.3954
* Wikipathways'WikiPathways_2019_Human.txt' from https://maayanlab.cloud/Enrichr/#stats

Data can be stored in data/raw for replication of the study.
