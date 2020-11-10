#Run nonmem models for treatment response for every treatment
# Input: outcome data (datanonmem.csv), natural growth curves from NONMEM (m0_natural_growth.tab), 
#        dictionary (dictionary.csv), nonmem models (m1.mod, m2.mod, m3.mod, m4.mod)
# Output: folder with all NONMEM models for all growth curves (folder: "run_dateymd_KDKR") 


# libraries
library(doParallel)

# Read data
dictio <- read.csv("data/dictionary.csv", header = TRUE, sep = ",", quote = "")
data <- read.csv("data/clean/datanonmem.csv")

#Extract most recent run for KG measurements
dirs <- sort(list.dirs(path = "results/NONMEM/", full.names = FALSE, recursive = FALSE), decreasing = TRUE)
runName <- grep("KG", dirs, value = TRUE)[1]

kg_import <- read.table(paste0("results/NONMEM/", runName, "/m0_natural_growth.tab"), header = TRUE, 
                        stringsAsFactors = FALSE, colClasses = "numeric", skip = 1)

# merge KG from NONMEM results with all treatment data
id_kg <- kg_import[, c("ID", "KG", "BASE")]
names(id_kg)[1] <- "Subject"
id_kg <- id_kg[!duplicated(id_kg$Subject), ]
merged <- merge(data, id_kg, by = "Subject")
alldata <- merged[order(merged$Subject, merged$Treatment, merged$Time, merged$Volume), ]


# Create folder for entire set of runs accross treatments
runName <- paste("run", format(Sys.time(), "%Y%m%d"), "KDKR", sep = "_")
dir.create(paste0("results/NONMEM/", runName))



setwd(paste0("data/NONMEM/", runName))


# Function for analysis in nonmem
runAnalysis <- function(d) {
  runID <- paste("treat_", d, sep="")
  treatment_folder <- paste0("results/NONMEM/", runName, "/", runID, "/")
  dir.create(treatment_folder)
  
  # copy model files to treatment folder
  file.copy(paste0("scripts/NONMEM/m", 1:4, ".mod"), treatment_folder)

  
  setwd(treatment_folder)
  
  # rename model files
  sapply(1:4, function(x) {
    file.rename(paste0("m", x,".mod"), paste0("treat_", d, "_", x, ".mod"))
  })
  
  # make a variable that saves the treatment name
  dname <- unique(as.character(dictio$Treat_name[dictio$TREAT == d]))
  
  # modify the data file thats used in the nonmem model, so that it contains the data from the right treatment
  data_treatment <- alldata[alldata$Treatment == d, ]
  write.table(x = data_treatment, "data_treatment.csv",
              quote = FALSE, row.names = F, sep = ",", na = ".")
  
  # add a print statement to keep track of for loops progress
  print(paste0("### Starting treatment ", d, ": ", dname, " ###"))
  # run nonmem model (only when nonmem is installed on PC!)
  system(paste("execute", paste("treat_", d, "_", 1:4, ".mod", sep = "", collapse = " "), 
               "-threads=4 -clean=2"), intern = FALSE, wait = TRUE)
  
  # rename output files
  sapply(1:4, function(x) {
    file.rename(paste0("m", x,".tab"), paste0("treat_", d, "_", x, ".tab"))
  })
 
  setwd("../../../..")
}

# Run function in parallel
no_cores <- floor((detectCores() - 1)/4) # already 4 threads per run
registerDoParallel(cores = no_cores)
cl <- makeCluster(no_cores, type="FORK") 
parLapply(cl, unique(alldata$Treatment), runAnalysis)

stopCluster(cl)