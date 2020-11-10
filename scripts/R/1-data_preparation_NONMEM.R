# Create input for nonmem and dictionary file
# Input: nm.3954-S2.xlsx (raw data)
# output: dictionary (data/dictionary.csv), nonmem formatted outcome PDX data (data/clean/datanonmem.csv and data/clean/data_u.csv)
library(readxl)
library(tidyverse)

raw_data <- read_excel("data/raw/nm.3954-S2.xlsx", sheet = "PCT raw data")

# rename columns for correct format
nonmem_data <- raw_data[, c("Model", "Treatment", "Days Post T0", "Volume (mm3)")]
names(nonmem_data)<- c("Subject", "Treatment", "Time", "Volume")

nonmem_data <- nonmem_data %>% 
  mutate(Dose = ifelse(Treatment == "untreated", 0, 1),
         Cmt = 2,
         TID = as.numeric(as.factor(paste0(Subject, as.factor(Treatment)))),
         Treat_name = as.character(Treatment),
         Treatment = as.numeric(as.factor(Treatment)),
         ID_name = Subject,
         Subject = as.numeric(as.factor(Subject)),
         AMT = 0, # dose column
         EVID = 0) # event ID column

dictionary_nonmem <- nonmem_data %>%
  select(c(Subject, Treatment, Treat_name, ID_name, TID)) %>%
  unique() %>%
  rename(ID = Subject, TREAT = Treatment, TID_original = TID) %>% 
  mutate(ID_TREAT = paste(ID, TREAT, sep = "_"))

#Save disctionary to be able to convert numeric IDs to names
write.table(dictionary_nonmem, file = "data/dictionary.csv", quote = FALSE, row.names = F, sep = ",", na = ".")

nonmem_data <- nonmem_data[, c('Subject', 'Treatment', 'Time', 'Volume', 'Dose', 'Cmt', 'TID', "AMT", "EVID")]

# add dosing + reset records to dataset for nonmem
nonmem_rows <- nonmem_data %>% select(-c(Time, Volume)) %>% unique() %>% 
  mutate(Time = 0, Cmt = 1, Volume = NA, EVID = 4, AMT = as.numeric(Dose == 1)) %>% 
  select(names(nonmem_data))

nonmem_data_final <- rbind(nonmem_data, nonmem_rows) %>% as.data.frame() %>% 
  arrange(TID, Time, -EVID) %>% 
  mutate(EVID = ifelse(AMT == 0 & EVID == 4,  3, EVID))

# add occasions column for nonmem (occasion is a treatment within subject)
nonmem_data_final$occ <- NA
for(s in unique(nonmem_data_final$Subject)) {
  c = 1
  for(t in unique(nonmem_data_final$Treatment[nonmem_data_final$Subject == s])) {
    nonmem_data_final$occ[nonmem_data_final$Subject == s & nonmem_data_final$Treatment == t] <- c
    c = c + 1
  }
}


# nonmem dataset all treatments
write.table(x = nonmem_data_final,file = "data/clean/datanonmem.csv", quote = FALSE, row.names = F, sep = ",", na = ".")


# create untreated dataset for calculating KG
untreated <- dictionary_nonmem$TREAT[dictionary_nonmem$Treat_name == "untreated"][1]
nonmem_data_untreated <- nonmem_data_final[nonmem_data_final$Treatment == untreated, ]

write.table(nonmem_data_untreated, file = "data/clean/data_u.csv", quote = FALSE, row.names = F, sep = ",", na = ".")
