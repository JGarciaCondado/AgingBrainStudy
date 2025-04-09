# Preprocess longitudinal data for ADNI

# Load necessary libraries
library(dplyr)

# Load baseline, longitudinal and deltas
data_baseline <- read.csv("~/Documents/Projects/AgingBrainStudy/data/final/adni/processed/baseline.csv", 
                          na.strings = c("", "NA")) 
data_long <- read.csv("~/Documents/Projects/AgingBrainStudy/data/final/adni/adni_long_040725.csv")
data_deltas <- read.csv("~/Documents/Projects/AgingBrainStudy/analysis/final/adni/ageml/model_age/predicted_age.csv") %>%
  rename(delta = delta_all)

# Join with baseline data that has already removed those with less than 2 timepoints
data_long <- merge(data_long, data_baseline, by = "ID")
data_long <- merge(data_long, data_deltas[, c('ID', 'delta')], by = "ID")

# Group data and calculate time from MRI for each test and PACC_change
data_long <- data_long %>%
  arrange(ID, Date) %>%
  group_by(ID) %>%
  mutate(Date = as.Date(Date), mri_date = as.Date(mri_date)) %>%
  mutate(time_mri = as.numeric(Date - mri_date) / 365.25) %>%
  mutate(PACC_change = zPACC - zPACC[1]) %>%
  mutate(delta_group = ifelse(delta < 0, 'younger', 'older'))

# Rename zPACC to PACC
data_long <- data_long %>%
  rename(PACC = zPACC) 

# Remove unnecessary columns
data_long <- data_long %>%
  select(-c(X, VISCODE, Cohort, PHASE, Date, Age, MonthsFromBaseline_raw,MonthsFromBaseline_rounded, lPACC, CDRGLOB, Diagnosis, Baseline_Dx, mri_date))

# Save the processed data
write.csv(data_long, "~/Documents/Projects/AgingBrainStudy/data/final/adni/processed/longitudinal.csv")
