# Import library

library(tidyverse)
############################################################################## SUN ############# 

clinical_SUN <- clinical_SUN %>% 
  # PSA is already there as "PSA"
  mutate(gleason = case_when(
    GLEASON_SCORE == "3+3"                                                      ~ "<=6",
    GLEASON_SCORE %in% c("3+4", "4+3")                                          ~ "7",
    GLEASON_SCORE %in% c("3+5", "4+4")                                          ~ "8",
    GLEASON_SCORE %in% c("4+5", "5+4", "5+5")                                   ~ "9-10"
  )) %>%
  mutate(tstage = case_when(
    str_detect(TNMSTAGE, "T2c")                                                 ~ "T2c",
    str_detect(TNMSTAGE, "T3a")                                                 ~ "T3a",
    str_detect(TNMSTAGE, "T3b")                                                 ~ "T3b",
    str_detect(TNMSTAGE, "T4")                                                  ~ "T4"
  )) %>%
  mutate(isup = case_when(
    GLEASON_SCORE == "3+3"                                                      ~ 1,
    GLEASON_SCORE == "3+4"                                                      ~ 2,
    GLEASON_SCORE == "4+3"                                                      ~ 3,
    GLEASON_SCORE %in% c("3+5", "4+4")                                          ~ 4,
    GLEASON_SCORE %in% c("4+5", "5+4", "5+5")                                   ~ 5
  )) %>%
  mutate(capra_psa = case_when( 
    PSA <= 6                                                                    ~ 0,
    (PSA > 6 & PSA <= 10)                                                       ~ 1,
    (PSA > 10 & PSA <= 20)                                                      ~ 2,
    (PSA > 20 & PSA <= 30)                                                      ~ 3,
    PSA > 30                                                                    ~ 4
  )) %>%
  mutate(capra_gleason = case_when(
    GLEASON_SCORE == "3+3"                                                      ~ 0,
    GLEASON_SCORE %in% c("3+4", "3+5")                                          ~ 1,
    GLEASON_SCORE %in% c("4+3", "4+4", "4+5", "5+4", "5+5")                     ~ 3
  )) %>%
  mutate(capra_tstage = case_when(
    str_detect(tstage, "(T1|T2)")                                               ~ 0,
    str_detect(tstage, "(T3|T4)")                                               ~ 1
  )) %>%
  mutate(capra_age = case_when(
    AGE < 50                                                                    ~ 0,
    AGE >= 50                                                                   ~ 1
  )) %>% 
  mutate(damico = case_when( #----------------------------------------Risk calculation start here, need to go in a function
    tstage %in% c("T2c", "T3a", "T3b", "T4") |
      PSA > 20 |
      gleason %in% c("8", "9-10")                                               ~ "High",
    (PSA > 10 & PSA <= 20) |
      gleason == "7"                                                            ~ "Intermediate",
    PSA <= 10 &
      gleason == "<=6"                                                          ~ "Low"
  )) %>%
  mutate(nice = case_when( 
    tstage %in% c("T2c", "T3", "T3a", "T3b", "T3c", "T4", "T4a", "T4b") | 
      PSA > 20 |
      gleason %in% c("8", "9-10")                                               ~ "High",
    tstage == "T2b" |
      gleason == "7" |
      (PSA >= 10 & PSA <= 20)                                                   ~ "Intermediate",
    tstage %in% c("T1", "T1a", "T1b", "T1c", "T2", "T2a") & 
      PSA < 10 & 
      gleason == "<=6"                                                          ~ "Low",
    TRUE                                                                        ~ NA_character_
  )) %>%
  mutate(eau = case_when(
    tstage %in% c("T2c", "T3", "T4", "T3a", "T3b", "T3c", "T4a", "T4b") |
      PSA > 20 |
      gleason %in% c("8", "9-10")                                               ~ "High",
    tstage == "T2b" |
      (PSA >= 10 & PSA <= 20) |
      gleason == "7"                                                            ~ "Intermediate",
    tstage %in% c("T1", "T1a", "T1b", "T1c", "T2a", "T2") &
      PSA < 10 &
      gleason == "<=6"                                                          ~ "Low",
    TRUE                                                                        ~ NA_character_
  )) %>%
  mutate(GUROC = case_when(
    PSA > 20 |
      gleason %in% c("8", "9-10") |
      tstage %in% c("T3", "T3a", "T3b", "T3c", "T4a", "T4b", "T4", "T4c")       ~ "High",
    # wriie low before intermediate to compare with "not otherwise low risk
    PSA <= 10 &
      gleason == "<=6" &
      tstage %in% c("T1", "T1a", "T1b" , "T1c", "T2", "T2a")                    ~ "Low", 
    PSA <= 20 &
      gleason %in% c("<=6", "7") &
      tstage %in% c("T1", "T1a", "T1b" , "T1c", "T2", "T2a", "T2b", "T2c")      ~ "Intermediate", 
    TRUE                                                                        ~ NA_character_
  )) %>% 
  # https://www.ncbi.nlm.nih.gov/pubmed/27483464/ 
  mutate(CPG = case_when( 
    (PSA > 20 & isup == 4) |
      (PSA > 20 & tstage %in% c("T3", "T3a", "T3b", "T3c")) |
      (isup == 4 &tstage %in% c("T3", "T3a", "T3b", "T3c")) |
      isup == 5 |
      tstage %in% c("T4", "T4a", "T4b", "T4c") ~ "Very High",
    PSA > 20 | isup == 4 | tstage %in% c("T3", "T3a", "T3b", "T3c")             ~ "High",
    ((PSA >= 10 & PSA <= 20) & isup  == 2 &
       tstage %in% c("T1", "T1a", "T1b", "T1c", "T2", "T2a", "T2b", "T2c")) |
      (isup == 3 &
         tstage %in% c("T1", "T1a", "T1b", "T1c", "T2", "T2a", "T2b", "T2c"))   ~ "Intermediate Unfavorable",                                             
    isup == 2 | (PSA >= 10 & PSA <= 20) &
      tstage %in% c("T1", "T1a", "T1b", "T1c", "T2", "T2a", "T2b", "T2c")       ~ "Intermediate Favorable",
    PSA < 10 & isup == 1 &
      tstage %in% c("T1", "T1a", "T1b" , "T1c", "T2", "T2a", "T2b", "T2c")      ~ "Low",
    TRUE                                                                        ~ NA_character_
  )) %>%
  mutate(capra_score = rowSums(select(.,capra_psa:capra_age), na.rm = TRUE)) %>%
  mutate(capra = case_when(
    capra_score %in% c(0:2)                                                     ~ "Low",
    capra_score %in% c(3:5)                                                     ~ "Intermediate",
    capra_score %in% c(6:10)                                                    ~ "High",
    TRUE                                                                        ~ NA_character_
  ))


############################################################################## CPC-GENE ############# 
clinical_CPC_GENE <- clinical_CPC_GENE %>%
  mutate(PSA = PSA_MOST_RECENT_RESULTS) %>% # Don't have anything else
  mutate(gleason = case_when(
    GLEASON_SCORE %in% c("3+3", "6")                                             ~ "<=6",
    GLEASON_SCORE %in% c("3+4", "3+4;5", "4+3","4+3;5", "7")                     ~ "7", 
    # Add "3+4;5" here (only 2 patients) and "4+3;5" (only 1 patient)
    GLEASON_SCORE %in% c("3+5", "4+4")                                           ~ "8"
  )) %>% 
  mutate(tstage = case_when(
    CLIN_T_STAGE == "T1B"                                                        ~ "T1b", # Found in dictionary but not in data...
    CLIN_T_STAGE == "T1C"                                                        ~ "T1c",
    CLIN_T_STAGE == "T2"                                                         ~ "T2",
    CLIN_T_STAGE == "T2A"                                                        ~ "T2a",
    CLIN_T_STAGE == "T2B"                                                        ~ "T2b",
    CLIN_T_STAGE == "T2C"                                                        ~ "T2c",
    CLIN_T_STAGE == "T3A"                                                        ~ "T3Aa",
    CLIN_T_STAGE == "T3B"                                                        ~ "T3Ab"
  )) %>%
  mutate(isup = case_when(
    GLEASON_SCORE %in% c("3+3", "6")                                             ~ 1,
    GLEASON_SCORE %in% c("3+4", "3+4;5") |
      (GLEASON_PATTERN_PRIMARY == "3" &
         GLEASON_PATTERN_SECONDARY == "4")                                       ~ 2,
    GLEASON_SCORE %in% c("4+3", "4+3;5") |
      (GLEASON_PATTERN_PRIMARY == "4" &
         GLEASON_PATTERN_SECONDARY == "3")                                       ~ 3
  )) %>%
  mutate(capra_psa = case_when( 
    PSA <= 6                                                                     ~ 0,
    (PSA > 6 & PSA <= 10)                                                        ~ 1,
    (PSA > 10 & PSA <= 20)                                                       ~ 2,
    (PSA > 20 & PSA <= 30)                                                       ~ 3,
    PSA > 30                                                                     ~ 4
  )) %>%
  mutate(capra_gleason = case_when(
    gleason == "<=6"                                                             ~ 0,
    gleason == "7"                                                               ~ 1,
    gleason ==  "8"                                                              ~ 3
  )) %>%
  mutate(capra_tstage = case_when(
    str_detect(tstage, "(T1|T2)")                                                ~ 0,
    str_detect(tstage, "(T3|T4)")                                                ~ 1
  )) %>%
  mutate(capra_age = case_when(
    AGE < 50                                                                     ~ 0,
    AGE >= 50                                                                    ~ 1
  )) %>% 
  mutate(damico = case_when( #----------------------------------------Risk calculation start here, need to go in a function
    tstage %in% c("T2c", "T3a", "T3b", "T4") |
      PSA > 20 |
      gleason %in% c("8", "9-10")                                               ~ "High",
    (PSA > 10 & PSA <= 20) |
      gleason == "7"                                                            ~ "Intermediate",
    PSA <= 10 &
      gleason == "<=6"                                                          ~ "Low"
  )) %>%
  mutate(nice = case_when( 
    tstage %in% c("T2c", "T3", "T3a", "T3b", "T3c", "T4", "T4a", "T4b") | 
      PSA > 20 |
      gleason %in% c("8", "9-10")                                               ~ "High",
    tstage == "T2b" |
      gleason == "7" |
      (PSA >= 10 & PSA <= 20)                                                   ~ "Intermediate",
    tstage %in% c("T1", "T1a", "T1b", "T1c", "T2", "T2a") & 
      PSA < 10 & 
      gleason == "<=6"                                                          ~ "Low",
    TRUE                                                                        ~ NA_character_
  )) %>%
  mutate(eau = case_when(
    tstage %in% c("T2c", "T3", "T4", "T3a", "T3b", "T3c", "T4a", "T4b") |
      PSA > 20 |
      gleason %in% c("8", "9-10")                                               ~ "High",
    tstage == "T2b" |
      (PSA >= 10 & PSA <= 20) |
      gleason == "7"                                                            ~ "Intermediate",
    tstage %in% c("T1", "T1a", "T1b", "T1c", "T2a", "T2") &
      PSA < 10 &
      gleason == "<=6"                                                          ~ "Low",
    TRUE                                                                        ~ NA_character_
  )) %>%
  mutate(GUROC = case_when(
    PSA > 20 |
      gleason %in% c("8", "9-10") |
      tstage %in% c("T3", "T3a", "T3b", "T3c", "T4a", "T4b", "T4", "T4c")       ~ "High",
    # wriie low before intermediate to compare with "not otherwise low risk
    PSA <= 10 &
      gleason == "<=6" &
      tstage %in% c("T1", "T1a", "T1b" , "T1c", "T2", "T2a")                    ~ "Low", 
    PSA <= 20 &
      gleason %in% c("<=6", "7") &
      tstage %in% c("T1", "T1a", "T1b" , "T1c", "T2", "T2a", "T2b", "T2c")      ~ "Intermediate", 
    TRUE                                                                        ~ NA_character_
  )) %>% 
  # https://www.ncbi.nlm.nih.gov/pubmed/27483464/ 
  mutate(CPG = case_when( 
    (PSA > 20 & isup == 4) |
      (PSA > 20 & tstage %in% c("T3", "T3a", "T3b", "T3c")) |
      (isup == 4 &tstage %in% c("T3", "T3a", "T3b", "T3c")) |
      isup == 5 |
      tstage %in% c("T4", "T4a", "T4b", "T4c") ~ "Very High",
    PSA > 20 | isup == 4 | tstage %in% c("T3", "T3a", "T3b", "T3c")             ~ "High",
    ((PSA >= 10 & PSA <= 20) & isup  == 2 &
       tstage %in% c("T1", "T1a", "T1b", "T1c", "T2", "T2a", "T2b", "T2c")) |
      (isup == 3 &
         tstage %in% c("T1", "T1a", "T1b", "T1c", "T2", "T2a", "T2b", "T2c"))   ~ "Intermediate Unfavorable",                                             
    isup == 2 | (PSA >= 10 & PSA <= 20) &
      tstage %in% c("T1", "T1a", "T1b", "T1c", "T2", "T2a", "T2b", "T2c")       ~ "Intermediate Favorable",
    PSA < 10 & isup == 1 &
      tstage %in% c("T1", "T1a", "T1b" , "T1c", "T2", "T2a", "T2b", "T2c")      ~ "Low",
    TRUE                                                                        ~ NA_character_
  )) %>%
  mutate(capra_score = rowSums(select(.,capra_psa:capra_age), na.rm = TRUE)) %>%
  mutate(capra = case_when(
    capra_score %in% c(0:2)                                                     ~ "Low",
    capra_score %in% c(3:5)                                                     ~ "Intermediate",
    capra_score %in% c(6:10)                                                    ~ "High",
    TRUE                                                                        ~ NA_character_
  ))

  ))

############################################################################## TCGA ############# 
# https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-entity-list&anchor=clinical

table(clinical_TCGA_333$CLINICAL_GLEASON)
table(clinical_TCGA_333$CLINICAL_GLEASON_CATEGORY)
table(clinical_TCGA_333$CLINICAL_GLEASON_SUM)


clinical_TCGA_333 <- clinical_TCGA_333 %>% 
  mutate(PSA = PREOPERATIVE_PSA) %>% 
  mutate(gleason = case_when(
    REVIEWED_GLEASON_SUM == 6                                                    ~ "<=6",
    REVIEWED_GLEASON_SUM == 7                                                    ~ "7",
    REVIEWED_GLEASON_SUM == 8                                                    ~ "8",
    REVIEWED_GLEASON_SUM %in% c(9, 10)                                           ~ "9-10"
  )) %>%
  # mutate(tstage = case_when(
  #   T_clinical == 1                                                              ~ "T1", # Not available for the 333 patients...
  #   T_clinical == 2                                                              ~ "T2",
  #   T_clinical == 3                                                              ~ "T3",
  #   T_clinical == 4                                                              ~ "T4"
  # )) %>%
  mutate(isup = case_when(
  REVIEWED_GLEASON %in% c("1+1", "1+2", "1+3", "2+1", "2+2", "2+3",
                          "3+1", "3+2", "3+3")                                   ~ 1,
  REVIEWED_GLEASON %in% c("3+5", "4+4", "4+5", "5+3", "5+4", "5+5")              ~ 4,
  REVIEWED_GLEASON == "3+4"                                                      ~ 2,
  REVIEWED_GLEASON == "4+3"                                                      ~ 3
  )) %>%
  mutate(capra_psa = case_when( # need to remove 323
    PSA <= 6                                                                     ~ 0,
    (PSA > 6 & PSA <= 10)                                                        ~ 1,
    (PSA > 10 & PSA <= 20)                                                       ~ 2,
    (PSA > 20 & PSA <= 30)                                                       ~ 3,
    PSA > 30                                                                     ~ 4,
    TRUE ~ NA_real_
  )) %>%
  mutate(capra_gleason = case_when(
    REVIEWED_GLEASON %in% c("1+1", "1+2", "1+3", "2+1", "2+2", "2+3",
                          "3+1", "3+2", "3+3")                                   ~ 0,
    REVIEWED_GLEASON %in% c("1+4", "1+5", "2+4", "2+5", 
                            "3+4", "3+5")                                        ~ 1,
    REVIEWED_GLEASON %in% c("4+1", "4+2", "4+3", "4+4", "4+5",
                            "5+1", "5+2", "5+3", "5+4", "5+5")                   ~ 3
  )) #%>%
  # mutate(capra_tstage = case_when(
  #   str_detect(tstage, "(T1|T2)")                                                ~ 0,
  #   str_detect(tstage, "(T3|T4)")                                                ~ 1
  # ))


############################################################################## MSKCC ############# 
clinical_MSKCC <- clinical_MSKCC
# no relevant info...

############################################################################## ICGC-FR #############
# ICGC doesn't have much clinical information (canada has only Tx stage)
clinical_ICGC_FR <- read_delim("/Users/colinccm/Downloads/donor.PRAD-FR.tsv", delim = "\t")

clinical_ICGC_FR <- clinical_ICGC_FR %>% 
  mutate(tstage = case_when(
    str_detect(donor_tumour_stage_at_diagnosis, "T1c")                          ~ "T1c",
    str_detect(donor_tumour_stage_at_diagnosis, "T2a")                          ~ "T2a",
    str_detect(donor_tumour_stage_at_diagnosis, "T2b")                          ~ "T2b",
    str_detect(donor_tumour_stage_at_diagnosis, "T2")                           ~ "T2"
  )) %>% 
  mutate(AGE = donor_age_at_diagnosis)

<<<<<<< HEAD
table(clinical_ICGC_FR$donor_tumour_stage_at_diagnosis)
clinical_CPC_GENE <- clinical_CPC_GENE %>%
  mutate(PSA = PSA_MOST_RECENT_RESULTS) %>% # Don't have anything else
  mutate(gleason = case_when(
    GLEASON_SCORE %in% c("3+3", "6")                                             ~ "<=6",
    GLEASON_SCORE %in% c("3+4", "3+4;5", "4+3","4+3;5", "7")                     ~ "7", 
    # Add "3+4;5" here (only 2 patients) and "4+3;5" (only 1 patient)
    GLEASON_SCORE %in% c("3+5", "4+4")                                           ~ "8"
  )) %>% 
  mutate(tstage = case_when(
    CLIN_T_STAGE == "T1B"                                                        ~ "T1b", # Found in dictionary but not in data...
    CLIN_T_STAGE == "T1C"                                                        ~ "T1c",
    CLIN_T_STAGE == "T2"                                                         ~ "T2",
    CLIN_T_STAGE == "T2A"                                                        ~ "T2a",
    CLIN_T_STAGE == "T2B"                                                        ~ "T2b",
    CLIN_T_STAGE == "T2C"                                                        ~ "T2c",
    CLIN_T_STAGE == "T3A"                                                        ~ "T3Aa",
    CLIN_T_STAGE == "T3B"                                                        ~ "T3Ab"
  )) %>%
  mutate(isup = case_when(
    GLEASON_SCORE %in% c("3+3", "6")                                             ~ 1,
    GLEASON_SCORE %in% c("3+4", "3+4;5") |
      (GLEASON_PATTERN_PRIMARY == "3" &
         GLEASON_PATTERN_SECONDARY == "4")                                       ~ 2,
    GLEASON_SCORE %in% c("4+3", "4+3;5") |
      (GLEASON_PATTERN_PRIMARY == "4" &
         GLEASON_PATTERN_SECONDARY == "3")                                       ~ 3
  )) %>%
  mutate(capra_psa = case_when( 
    PSA <= 6                                                                     ~ 0,
    (PSA > 6 & PSA <= 10)                                                        ~ 1,
    (PSA > 10 & PSA <= 20)                                                       ~ 2,
    (PSA > 20 & PSA <= 30)                                                       ~ 3,
    PSA > 30                                                                     ~ 4
  )) %>%
  mutate(capra_gleason = case_when(
    gleason == "<=6"                                                             ~ 0,
    gleason == "7"                                                               ~ 1,
    gleason ==  "8"                                                              ~ 3
  )) %>%
  mutate(capra_tstage = case_when(
    str_detect(tstage, "(T1|T2)")                                                ~ 0,
    str_detect(tstage, "(T3|T4)")                                                ~ 1
  )) %>%
  mutate(capra_age = case_when(
    AGE < 50                                                                     ~ 0,
    AGE >= 50                                                                    ~ 1
  ))
=======
clinical_ICGC_UK <- read_delim("/Users/colinccm/Downloads/donor.PRAD-UK.tsv", delim = "\t")
>>>>>>> Update loading and recoding

clinical_ICGC_UK <- clinical_ICGC_UK %>% 
  mutate(tstage = case_when(
    str_detect(donor_tumour_stage_at_diagnosis, "T1c")                          ~ "T1c",
    str_detect(donor_tumour_stage_at_diagnosis, "T2a")                          ~ "T2a",
    str_detect(donor_tumour_stage_at_diagnosis, "T2b")                          ~ "T2b",
    str_detect(donor_tumour_stage_at_diagnosis, "T2")                           ~ "T2"
  )) %>% 
  mutate(AGE = donor_age_at_diagnosis)
