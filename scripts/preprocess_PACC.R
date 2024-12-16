library(stringr)
library(tidyr)
library(dplyr)
library(janitor)
library(plyr)
library(readxl)

adqs <- read.csv("~/Documents/Projects/AgingBrainStudy/data/cognition/ADQS.csv", na.strings = c("", " ", "NaN", "NA", NA))
visits <- read.csv("~/Documents/Projects/AgingBrainStudy/data/Visits/SV.csv", na.strings = c("", " ", "NaN", "NA", NA))
subj <- read.csv("~/Documents/Projects/AgingBrainStudy/data/SUBJINFO.csv", na.strings = c("", " ", "NaN", "NA", NA))
amyloid <- read.csv("~/Documents/Projects/AgingBrainStudy/data/pet_imaging/imaging_SUVR_amyloid.csv", na.strings = c("", " ", "NaN", "NA", NA))
pet <- read.csv("~/Documents/Projects/AgingBrainStudy/data/pet_imaging/imaging_PET_VA.csv", na.strings = c("", " ", "NaN", "NA", NA))

# Visits ------------------------------------------------------------------

visits <- select(visits, BID, VISCODE = VISITCD, SVSTDTC_DAYS_CONSENT)

# Rescreens ---------------------------------------------------------------

rescreens <- subj %>% 
  filter(!is.na(PREVBID)) %>% 
  dplyr::rename(BID1 = PREVBID, BID2 = BID) %>% 
  select(BID1, BID2)

ids <- amyloid %>%
  filter(brain_region != '' & VISCODE == 2) %>%
  pivot_wider(id_cols = 'BID', names_from = brain_region, values_from = suvr_cer) %>% 
  left_join(pet, by='BID') %>%
  left_join(rescreens, by = c('BID' = 'BID1')) %>%
  mutate( # update PET BIDs to second BID if necessary
    BID = case_when( 
      !is.na(BID2) ~ BID2,
      TRUE ~ BID)) %>%
  arrange(BID, BID2) %>% 
  filter(!duplicated(BID, fromLast = TRUE)) %>% 
  select(BID)


# PACC --------------------------------------------------------------------

tests <- c("MMSE", "DIGIT", "LMIIa", "FCSRT96", "PACC") # select tests to keep under QSTESTCD

data <- adqs %>%
  filter(BID %in% ids$BID) %>%
  filter(QSTESTCD %in% tests) %>%
  pivot_wider(names_from = QSTESTCD, values_from = c(QSSTRESN, QSVERSION)) %>%
  group_by(BID, VISITCD) %>%
  fill(all_of(contains(tests)), .direction = "updown") %>%
  select(SUBSTUDY, BID, TX, AGEYR, SEX, EDCCNTU, RACE, ETHNIC,
         APOE = APOEGNPRSNFLG, APOEGN, NP_VISCODE = VISITCD, VISIT, EPOCH, 
         contains(paste0("_", tests)), pacc_consent_time = QSDTC_DAYS_CONSENT) %>%
  distinct() %>%
  remove_empty("cols") %>% # remove empty cols - incase test doesn't have a test version var
  rename_with(.col = starts_with("QSSTRESN"), ~ str_replace(., "QSSTRESN_", "")) %>% # remove prefix from test vars
  group_by(BID) %>%
  dplyr::mutate(
    pacc_consent_time = pacc_consent_time / 365.25,
    pacc_time = pacc_consent_time - pacc_consent_time[1],
    pacc_age = AGEYR + pacc_time,
    bl_pacc_age = min(pacc_age, na.rm = TRUE),
    SEX = case_when(
      SEX == 1 ~ "Female",
      SEX == 2 ~ "Male", 
      TRUE ~ NA),
    ETHNIC = case_when(
      ETHNIC == 50 ~ "Hispanic or Latino",
      ETHNIC == 56 ~ "Not Hispanic or Latino",
      TRUE ~ NA),
    RACE = case_when(
      RACE == 1 ~ "White",
      RACE == 2 ~ "Black or African American",
      RACE == 58 ~ "Asian",
      RACE == 79 ~ "Native Hawaiian or Other Pacific Islander",
      RACE == 84 ~ "American Indian or Alaskan Native",
      RACE == 100 ~ "More than one race",
      TRUE ~ NA),
    APOE = case_when(
      APOE == 1 ~ "e4+",
      APOE == 0 ~ "e4-",
      TRUE ~ NA)
  ) %>%
  filter(!is.na(pacc_consent_time)) %>%
  arrange(BID, pacc_consent_time) %>%
  group_by(BID) %>%
  add_tally(name = "pacc_tp", wt = !is.na(PACC))

# DOSE --------------------------------------------------------------------

dose <- read_excel("~/Documents/Projects/AgingBrainStudy/data/dose/A4_DOSE_blind.xlsx")
dose_visit <- read.csv("~/Documents/Projects/AgingBrainStudy/data/Visits/visit_dose.csv", na.strings = c("", " ", "NaN", "NA", NA))

dose_visit <- dose_visit %>%
  dplyr::rename(BID = USUBJID, VISCODE = VISITCD) %>%
  dplyr::select(-c(VISIT, AVISIT))

dose <- merge(dose, dose_visit, by = c("BID", "VISCODE"), all.x = TRUE)
dose <- merge(dose, subj %>% dplyr::select(BID, TX), by = "BID", all.x = TRUE)

# Fill in epoch type in between
dose <- dose %>%
  group_by(BID) %>%
  arrange(BID, VISCODE) %>%
  tidyr::fill(EPOCH, .direction = "down")

dose_trt <- dose %>%
  mutate(Dose = case_when(
    TX == "Placebo" & EPOCH != "OPEN LABEL TREATMENT" ~ 0,
    EPOCH == "SCREENING" ~ 0,
    is.na(BLINDDOSE) ~ 0,
    TRUE ~ BLINDDOSE
  )) %>%
  group_by(BID) %>%
  arrange(BID, VISCODE) %>%
  dplyr::mutate(Cumulative_Dose = cumsum(Dose)) %>%
  mutate(Cumulative_Dose_Scaled = Cumulative_Dose / 10000) %>%
  dplyr::select(BID, VISCODE, Cumulative_Dose_Scaled)
  

# Merge data ----------------------------------------------------------------

data <- data %>%
  dplyr::left_join(dose_trt, by = c("BID" = "BID", "NP_VISCODE" = "VISCODE")) %>%
  arrange(BID, pacc_consent_time)

# Fill in for other studies
data <- data %>%
  dplyr::mutate(Cumulative_Dose_Scaled = case_when(
    SUBSTUDY != "A4" ~ 0,
    SUBSTUDY == "A4" & EPOCH == "SCREENING" ~ 0,
    TRUE ~ Cumulative_Dose_Scaled
  ))

# Save as CSV
write.csv(data, "~/Documents/Projects/AgingBrainStudy/data/A4/longitudinal_PACC.csv", row.names = FALSE)

