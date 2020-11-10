#Run nonmem model for natural growth curve (untreated)
# Input: outcome data (data_u.csv), dictionary (dictionary.csv)
# Output: folder with NONMEM model for natural growth curves (results/NONMEM/KG_untreated)

runName <- paste("run", format(Sys.time(), "%Y%m%d"), "KG_untreated", sep = "_")

#Create result directories
if (!dir.exists("results")) {
  dir.create("results")
} 
if (!dir.exists("results/NONMEM")) {
  dir.create("results/NONMEM")
}

dir.create(paste0("results/NONMEM/", runName))


file.copy("scripts/NONMEM/m0_natural_growth.mod", paste0("results/NONMEM/", runName))
file.copy("data/clean/data_u.csv", paste0("results/NONMEM/", runName))

setwd(paste0("results/NONMEM/", runName))

#command for running nonmem
command <- "execute m0_natural_growth -threads=4 -clean=2"

print(paste("starting KG estimation at", Sys.time()))
system(command, show.output.on.console = FALSE, intern = FALSE, wait = TRUE)
print(paste("ended at", Sys.time()))

print("removing extra data_u.csv copy")
file.remove("data_u.csv")


print("done!")