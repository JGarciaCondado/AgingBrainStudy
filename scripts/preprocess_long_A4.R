# Preprocess longitudinal data for A4

library(stringr)
library(tidyr)
library(dplyr)
library(janitor)
library(plyr)
library(readxl)

# Load baseline, longitudinal and deltas
data_baseline <- read.csv("~/Documents/Projects/AgingBrainStudy/data/final/a4/processed/baseline.csv", 
                          na.strings = c("", "NA")) 
data_long <- read.csv("~/Documents/Projects/AgingBrainStudy/data/final/a4/a4_long_041125.csv")
data_deltas <- read.csv("~/Documents/Projects/AgingBrainStudy/analysis/final/a4/ageml/model_age/predicted_age.csv") %>%
  dplyr::rename(delta = delta_all)

# Join with baseline data that has already removed those with less than 2 timepoints
data_long <- merge(data_long, data_baseline, by = "ID")
data_long <- merge(data_long, data_deltas[, c('ID', 'delta')], by = "ID")

# Group data and calculate time from MRI for each test and PACC_change
data_long <- data_long %>%
  arrange(ID, MonthsFromBaseline_raw) %>%
  group_by(ID) %>%
  filter(!is.na(zPACC)) %>%
  mutate(time_mri = (MonthsFromBaseline_raw - time_baseline_to_mri) / 12) %>%
  mutate(PACC_change = zPACC - zPACC[1]) %>%
  mutate(delta_group = ifelse(delta < 0, 'younger', 'older'))

# Rename zPACC to PACC
data_long <- data_long %>%
  dplyr::rename(PACC = zPACC) 

# Remove unnecessary columns
data_long <- data_long %>%
  select(-c(X, MonthsFromBaseline_raw, lPACC, CDRSUM, CDRGLOB, time_baseline_to_mri))

# DOSE --------------------------------------------------------------------

dose <- read_excel("~/Documents/Projects/AgingBrainStudy/data/final/a4/raw_data/A4_DOSE_blind.xlsx") %>%
  dplyr::rename(ID = BID) %>% 
  dplyr::mutate(ID = paste0(ID, "_a4"))
dose_visit <- read.csv("~/Documents/Projects/AgingBrainStudy/data/final/a4/raw_data/visit_dose.csv", na.strings = c("", " ", "NaN", "NA", NA))

dose_visit <- dose_visit %>%
  dplyr::rename(ID = USUBJID, VISCODE = VISITCD) %>%
  dplyr::select(-c(VISIT, AVISIT)) %>% 
  dplyr::mutate(ID = paste0(ID, "_a4"))

dose <- merge(dose, dose_visit, by = c("ID", "VISCODE"), all.x = TRUE)
dose <- merge(dose, data_baseline %>% dplyr::select(ID, TX), by = "ID", all.x = TRUE)

# Fill in epoch type in between
dose <- dose %>%
  group_by(ID) %>%
  arrange(ID, VISCODE) %>%
  tidyr::fill(EPOCH, .direction = "down")

dose_trt <- dose %>%
  mutate(Dose = case_when(
    TX == "Placebo" & EPOCH != "OPEN LABEL TREATMENT" ~ 0,
    EPOCH == "SCREENING" ~ 0,
    is.na(BLINDDOSE) ~ 0,
    TRUE ~ BLINDDOSE
  )) %>%
  group_by(ID) %>%
  arrange(ID, VISCODE) %>%
  dplyr::mutate(Cumulative_Dose = cumsum(Dose)) %>%
  mutate(Cumulative_Dose_Scaled = Cumulative_Dose / 10000) %>%
  dplyr::select(ID, VISCODE, Cumulative_Dose_Scaled, TX)
  

# Merge data ----------------------------------------------------------------

data_long <- data_long %>%
  dplyr::left_join(dose_trt %>% dplyr::select(ID, VISCODE, Cumulative_Dose_Scaled), by = c("ID" = "ID", "VISCODE" = "VISCODE")) %>%
  arrange(ID, time_mri)

# Fill in for other studies
data_long <- data_long %>%
  dplyr::mutate(Cumulative_Dose_Scaled = case_when(
    cohort == "LEARN/SF" ~ 0,
    cohort == "A4 Treated" | cohort == "A4 Placebo" & VISCODE == 1 ~ 0,
    TRUE ~ Cumulative_Dose_Scaled
  ))

# Save as CSV
write.csv(data_long, "~/Documents/Projects/AgingBrainStudy/data/final/a4/processed/longitudinal.csv", row.names = FALSE)

