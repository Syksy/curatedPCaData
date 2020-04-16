# Import library

library(tidyverse)

############################################################################## TCGA ############# 
# For now
TCGA_clinic <- read.delim("data-raw/prad_tcga_curated_pdata.txt")

# TCGA recoding
# https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-entity-list&anchor=clinical

colnames(TCGA_clinic)

TCGA <- TCGA_clinic %>% 
  mutate(gleason = case_when(
    gleason_grade == 6                          ~ "<=6",
    gleason_grade == 7                          ~ "7",
    gleason_grade == 8                          ~ "8",
    gleason_grade %in% c(9, 10)                 ~ "9-10"
  )) %>%
  mutate(tstage = case_when(
    T_clinical == 1 ~ "T1",
    T_clinical == 2 ~ "T2",
    T_clinical == 3 ~ "T3",
    T_clinical == 4 ~ "T4"
  )) %>%
  mutate(isup = case_when(
    grade_group == "<=6" ~ 1,
    grade_group == ">=8" ~ 4,
    grade_group == "3+4" ~ 2,
    grade_group == "4+3" ~ 3
  )) %>%
  mutate(capra_psa = case_when( # need to remove 323
    psa <= 6               ~ 0,
    (psa > 6 & psa <= 10)  ~ 1,
    (psa > 10 & psa <= 20) ~ 2,
    (psa > 20 & psa <= 30) ~ 3,
    psa > 30               ~ 4,
    TRUE ~ NA_real_
  )) %>%
  mutate(capra_gleason = case_when(
    gleason__major %in% c(1:3) &
      gleason_minor %in% c(1:3)        ~ 0,
    gleason__major %in% c(1:3) &
      gleason_minor %in% c(4:5)        ~ 1,
    gleason__major %in% c(4:5) &
      gleason_minor %in% c(1:5)        ~ 3,
    TRUE                               ~ NA_real_
  )) %>%
  mutate(capra_tstage = case_when(
    tstage %in% c("T1", "T2")  ~ 0,
    tstage %in% c("T3", "T4")  ~ 1,
    TRUE                       ~ NA_real_
  ))

