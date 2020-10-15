library(magrittr)
library(dplyr)

initial_curated_df <- function(
  df_rownames,
  template_name
){
  # import template - hand written by ...? 
  template <- utils::read.csv(template_name, as.is=TRUE)
  output <- matrix(NA,
                   ncol = nrow(template),
                   nrow = length(df_rownames))
  colnames(output) <- template$col.name
  rownames(output) <- df_rownames
  output <- data.frame(output)
  for (i in 1:ncol(output)){
    class(output[,i]) <- template[i,"var.class"]
  }
  output$sample_name <- df_rownames
  return(output)
}

###############################################################################
#  _________  ________  ________  ________     
# |\___   ___\\   ____\|\   ____\|\   __  \    
# \|___ \  \_\ \  \___|\ \  \___|\ \  \|\  \   
#      \ \  \ \ \  \    \ \  \  __\ \   __  \  
#       \ \  \ \ \  \____\ \  \|\  \ \  \ \  \ 
#        \ \__\ \ \_______\ \_______\ \__\ \__\
#         \|__|  \|_______|\|_______|\|__|\|__|
#  
###############################################################################
# get the first file with more sample specfic information ---------------------
## We are downloading all 499 cases and not the analytic subset of 333
## Will download and organize all the data here and then subset in later stages

mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
uncurated <- cgdsr::getClinicalData(mycgds, caseList="prad_tcga_sequenced")


# create the curated object
curated <- initial_curated_df(
  df_rownames = rownames(uncurated),
  template_name="data-raw/template_prad.csv")

curated <- curated %>% 
  dplyr::mutate(study_name = "TCGA, provisional") %>%
  dplyr::mutate(sample_name = row.names(uncurated)) %>% 
  dplyr::mutate(frozen_ffpe = uncurated$IS_FFPE) %>%
  # do we want the same value for missing and no? 
  # do we want true NA or character NA?
  dplyr::mutate(frozen_ffpe = dplyr::case_when(
    frozen_ffpe %in% c("NO", "[Not Available]") ~ "NA",
    frozen_ffpe == "YES" ~ "ffpe",
    TRUE ~ frozen_ffpe
  )) %>%
  dplyr::mutate(sample_type = uncurated$SAMPLE_TYPE) %>% 
  dplyr::mutate(patient_id = stringr::str_sub(sample_name,1,12)) %>%
  dplyr::mutate(alt_sample_name = uncurated$OTHER_SAMPLE_ID) %>%
  dplyr::mutate(gleason_grade = uncurated$GLEASON_SCORE) %>%
  dplyr::mutate(gleason_major = uncurated$GLEASON_PATTERN_PRIMARY) %>%
  dplyr::mutate(gleason_minor = uncurated$GLEASON_PATTERN_SECONDARY) %>%
  dplyr::mutate(source_of_gleason = "biopsy") %>%
  dplyr::mutate(grade_group = dplyr::case_when(
    gleason_grade == 6 ~ "<=6",
    gleason_major == 3 & gleason_minor == 4 ~ "3+4",
    gleason_major == 4 & gleason_minor == 3 ~ "4+3",
    gleason_grade %in% 8:10 ~ ">=8",
    TRUE ~ "NA"
  )) %>% 
  dplyr::mutate(zone_of_origin = uncurated$PRIMARY_SITE) %>%
  dplyr::mutate(zone_of_origin = dplyr::case_when(
    zone_of_origin == "Peripheral Zone" ~ "peripheral",
    zone_of_origin == "Overlapping / Multiple Zones" ~ "mixed",
    zone_of_origin == "Central Zone" ~ "central",
    zone_of_origin == "Transition Zone" ~ "transitional",
    zone_of_origin == "[Not Available]" ~ "NA",
    TRUE ~ "NA"
  )) %>%
  dplyr::mutate(year_diagnosis = uncurated$INITIAL_PATHOLOGIC_DX_YEAR) %>%
  dplyr::mutate(overall_survival_status = uncurated$OS_STATUS) %>%
  dplyr::mutate(overall_survival_status = dplyr::case_when(
    is.na(overall_survival_status) ~ NA_real_,
    overall_survival_status == "DECEASED" ~ 1,
    overall_survival_status != "DECEASED" ~ 0
  )) %>% 
  dplyr::mutate(days_to_overall_survival = as.numeric(uncurated$OS_MONTHS) * 30.5) %>%
  dplyr::mutate(disease_specific_recurrence_status = uncurated$DFS_STATUS) %>% 
  dplyr::mutate(disease_specific_recurrence_status = dplyr::case_when(
    disease_specific_recurrence_status == "[Not Available]" ~ NA_real_,
    disease_specific_recurrence_status == "Recurred/Progressed" ~ 1,
    disease_specific_recurrence_status == "DiseaseFree" ~ 0
  )) %>% 
  dplyr::mutate(days_to_disease_specific_recurrence = as.numeric(uncurated$DFS_MONTHS) * 30.5) %>% 
  dplyr::mutate(psa = uncurated$PSA_MOST_RECENT_RESULTS) %>% 
  dplyr::mutate(age_at_initial_diagnosis = uncurated$AGE) %>%
  dplyr::mutate(M_stage = uncurated$CLIN_M_STAGE) %>% 
  dplyr::mutate(M_stage = dplyr::case_when(
    M_stage == "M0" ~ 0,
    M_stage %in% c("M1a", "M1b", "M1c") ~ 1,
    TRUE ~ NA_real_
  )) %>% 
  # stringr:: commands return true NA not character NA
  dplyr::mutate(M_substage = stringr::str_sub(uncurated$CLIN_M_STAGE, 3, 3)) %>% 
  dplyr::mutate(T_clinical = readr::parse_number(uncurated$CLIN_T_STAGE)) %>% 
  dplyr::mutate(T_substage_clinical = stringr::str_extract(uncurated$CLIN_T_STAGE, "[a-c]+")) %>%
  dplyr::mutate(T_pathological = readr::parse_number(uncurated$PATH_T_STAGE)) %>%
  dplyr::mutate(T_substage_pathological = stringr::str_extract(uncurated$PATH_T_STAGE, "[a-c]+")) %>%
  dplyr::mutate(race = uncurated$RACE) %>%
  dplyr::mutate(race = dplyr::case_when(
    race == "WHITE" ~ "caucasian",
    race == "ASIAN" ~ "asian",
    race == "BLACK OR AFRICAN AMERICAN" ~ "african_american",
    TRUE ~ "NA"
  )) %>% 
  dplyr::mutate(tumor_purity_pathology = uncurated$TUMOR_CONTENT)

clinical_tcga <- curated

save(clinical_tcga, file = "data-raw/clinical_tcga.RData")

###############################################################################
#  ________  ___  ___  ________           _______  _________            ________  ___              
# |\   ____\|\  \|\  \|\   ___  \        |\  ___ \|\___   ___\         |\   __  \|\  \             
# \ \  \___|\ \  \\\  \ \  \\ \  \       \ \   __/\|___ \  \_|         \ \  \|\  \ \  \            
#  \ \_____  \ \  \\\  \ \  \\ \  \       \ \  \_|/__  \ \  \           \ \   __  \ \  \           
#   \|____|\  \ \  \\\  \ \  \\ \  \       \ \  \_|\ \  \ \  \ ___       \ \  \ \  \ \  \____      
#     ____\_\  \ \_______\ \__\\ \__\       \ \_______\  \ \__\\__\       \ \__\ \__\ \_______\    
#    |\_________\|_______|\|__| \|__|        \|_______|   \|__\|__|        \|__|\|__|\|_______|    
#                                                                                
############################################################################### 

gse <- GEOquery::getGEO("GSE25136", GSEMatrix = TRUE)

uncurated <- Biobase::pData(gse[[1]])

curated <- initial_curated_df(
  df_rownames = rownames(uncurated),
  template_name="data-raw/template_prad.csv")

curated <- curated %>% 
  dplyr::mutate(study_name = "Sun, et al.") %>%
  dplyr::mutate(sample_name = row.names(uncurated)) %>% 
  dplyr::mutate(patient_id = row.names(uncurated)) %>%
  dplyr::mutate(alt_sample_name = stringr::str_remove(uncurated$title,
                                                      "Prostate cancer primary tumor ")) %>%
  dplyr::mutate(disease_specific_recurrence_status = stringr::str_remove(uncurated$characteristics_ch1.1,
                                                                         "recurrence status: ")) %>% 
  dplyr::mutate(disease_specific_recurrence_status = dplyr::case_when(
    disease_specific_recurrence_status == "Non-Recurrent" ~ 0,
    disease_specific_recurrence_status == "Recurrent" ~ 1,
    TRUE ~ NA_real_
  )) %>%
  dplyr::mutate(sample_type = stringr::str_remove(uncurated$`tissue:ch1`,
                                                  "Prostate cancer ")) %>% 
  dplyr::mutate(sample_type = stringr::str_remove(sample_type, " tumor"))

clinical_sun <- curated

save(clinical_sun, file = "data-raw/clinical_sun.RData")

###############################################################################
#  _________  ________      ___    ___ ___       ________  ________     
# |\___   ___\\   __  \    |\  \  /  /|\  \     |\   __  \|\   __  \    
# \|___ \  \_\ \  \|\  \   \ \  \/  / | \  \    \ \  \|\  \ \  \|\  \   
#      \ \  \ \ \   __  \   \ \    / / \ \  \    \ \  \\\  \ \   _  _\  
#       \ \  \ \ \  \ \  \   \/  /  /   \ \  \____\ \  \\\  \ \  \\  \| 
#        \ \__\ \ \__\ \__\__/  / /      \ \_______\ \_______\ \__\\ _\ 
#         \|__|  \|__|\|__|\___/ /        \|_______|\|_______|\|__|\|__|
#                         \|___|/                                       
#   
############################################################################### 

# GEOquery for GSE21032 is BUSTED 
gse <- GEOquery::getGEO("GSE21032", GSEMatrix = TRUE)
## TDL: Split fetch into GSE21034 (GEX) and GSE21035 (CNA) for making the map
gse_gex <- GEOquery::getGEO("GSE21034", GSEMatrix = TRUE)
gse_cna <- GEOquery::getGEO("GSE21035", GSEMatrix = TRUE)
gse_gex_1 <- gse_gex[[1]]@phenoData@data[,'sample id:ch1',drop=FALSE] # Exon part
gse_gex_2 <- gse_gex[[2]]@phenoData@data[,'sample id:ch1',drop=FALSE] # Transcript part
gse_cna_1 <- gse_cna[[1]]@phenoData@data[,'sample id:ch1',drop=FALSE]

map <- data.frame(
	assay = c(rep("cna", nrow(gse_cna_1)), rep("gex", nrow(gse_gex_1)+nrow(gse_gex_2))),
	colname = c(rownames(gse_cna_1), rownames(gse_gex_1), rownames(gse_gex_2)),
	primary = c(gse_cna_1[,1], gse_gex_1[,1], gse_gex_2[,1])
)
# grep down into patient samples only
map <- map[grep("PCA", map[,"primary"]),]


uncurated <- Biobase::pData(gse[[2]])
## TDL: Above seems to only pull the GSMs for aCGH (CNA) data

# mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
# uncurated <- cgdsr::getClinicalData(mycgds, caseList = "prad_mskcc_all")

curated <- initial_curated_df(
  df_rownames = rownames(uncurated),
  template_name="data-raw/template_prad.csv")

# uncurated_edited <- uncurated %>% 
#   mutate(sample_type = case_when(
#     stringr::str_starts(characteristics_ch1.3, "tumor type: ") ~ 
#       stringr::str_remove(characteristics_ch1.3, "tumor type: "),
#     stringr::str_starts(characteristics_ch1.4, "tumor type: ") ~ 
#       stringr::str_remove(characteristics_ch1.4, "tumor type: "),
#     TRUE ~ "Unknown"
#   )) 

curated <- curated %>% 
  dplyr::mutate(study_name = "Taylor, et al.") %>% 
  ## TDL: Trying to use the GSM-derived mapping of multi-platform GSM### to unique PCA### sample ids in the sampleMap instead of sample_name
  #dplyr::mutate(sample_name = row.names(uncurated)) %>% 
  dplyr::mutate(sample_name = uncurated$`sample id:ch1`) %>%
  ## End of changes
  dplyr::mutate(patient_id = uncurated$`sample id:ch1`) %>%
  dplyr::mutate(sample_type = tolower(uncurated$`tumor type:ch1`)) %>%
  dplyr::mutate(sample_type = dplyr::case_when(
    sample_type == "primary tumor" ~ "primary",
    sample_type == "cell line" ~ "cell.line",
    TRUE ~ sample_type
  )) %>% 
  dplyr::mutate(gleason_grade = as.numeric(uncurated$`biopsy_gleason_grade:ch1`)) %>%
  dplyr::mutate(T_clinical = stringr::str_sub(uncurated$`clint_stage:ch1`,2,2)) %>%
  dplyr::mutate(T_clinical = as.numeric(T_clinical)) %>% 
  dplyr::mutate(T_substage_clinical = tolower(stringr::str_sub(uncurated$`clint_stage:ch1`,3,3))) %>% 
  dplyr::mutate(T_substage_clinical = dplyr::case_when(
    T_substage_clinical == "" ~ NA_character_,
    TRUE ~ T_substage_clinical
  )) %>% 
  dplyr::mutate(T_pathological = as.numeric(stringr::str_sub(uncurated$`pathological_stage:ch1`,2,2))) %>% 
  dplyr::mutate(T_substage_pathological = tolower(stringr::str_sub(uncurated$`pathological_stage:ch1`,3,3))) %>% 
  dplyr::mutate(T_substage_pathological = dplyr::case_when(
    T_substage_pathological == "" ~ NA_character_,
    TRUE ~ T_substage_pathological
  )) 

clinical_taylor <- curated

save(clinical_taylor, file = "data-raw/clinical_taylor.RData")

###############################################################################
# ___  ___  ___  _______   ________  ________  ________       ___    ___ _____ ______   ___  ___  ________      
# |\  \|\  \|\  \|\  ___ \ |\   __  \|\   __  \|\   ___  \    |\  \  /  /|\   _ \  _   \|\  \|\  \|\   ____\     
# \ \  \\\  \ \  \ \   __/|\ \  \|\  \ \  \|\  \ \  \\ \  \   \ \  \/  / | \  \\\__\ \  \ \  \\\  \ \  \___|_    
#  \ \   __  \ \  \ \  \_|/_\ \   _  _\ \  \\\  \ \  \\ \  \   \ \    / / \ \  \\|__| \  \ \  \\\  \ \_____  \   
#   \ \  \ \  \ \  \ \  \_|\ \ \  \\  \\ \  \\\  \ \  \\ \  \   \/  /  /   \ \  \    \ \  \ \  \\\  \|____|\  \  
#    \ \__\ \__\ \__\ \_______\ \__\\ _\\ \_______\ \__\\ \__\__/  / /      \ \__\    \ \__\ \_______\____\_\  \ 
#     \|__|\|__|\|__|\|_______|\|__|\|__|\|_______|\|__| \|__|\___/ /        \|__|     \|__|\|_______|\_________\
#                                                           \|___|/                                 \|_________|
#   
###############################################################################  

gse <- GEOquery::getGEO("GSE54691", GSEMatrix = TRUE)

uncurated <- Biobase::pData(gse[[1]])

curated <- initial_curated_df(
  df_rownames = rownames(uncurated),
  template_name="data-raw/template_prad.csv")

curated <- curated %>% 
  dplyr::mutate(study_name = "Hieronymus, et al.") %>% 
  dplyr::mutate(sample_name = row.names(uncurated)) %>% 
  dplyr::mutate(patient_id = row.names(uncurated)) %>%
  dplyr::mutate(overall_survival_status = dplyr::case_when(
    is.na(uncurated$`survivalevent:ch1`) ~ 0,
    TRUE ~ 1
  )) %>%
  dplyr::mutate(days_to_overall_survival = 
                  as.numeric(uncurated$`survival_or_followup_time_months:ch1`)*30.5) %>% 
  dplyr::mutate(age_at_initial_diagnosis = 
                  as.numeric(uncurated$`dxage:ch1`)) %>% 
  dplyr::mutate(gleason_grade = as.numeric(uncurated$`pathggs:ch1`)) %>%
  dplyr::mutate(gleason_minor = as.numeric(uncurated$`pathgg2:ch1`)) %>%
  dplyr::mutate(gleason_major = as.numeric(uncurated$`pathgg1:ch1`)) %>%
  dplyr::mutate(grade_group = dplyr::case_when(
    gleason_grade == 6 ~ "<=6",
    gleason_grade > 8 ~ ">=8",
    gleason_grade == 7 & gleason_major == "3" ~ "3+4",
    gleason_grade == 7 & gleason_major == "4" ~ "4+3",
  )) %>% 
  dplyr::mutate(source_of_gleason = "prostatectomy") %>% 
  dplyr::mutate(T_pathological = readr::parse_number(uncurated$`pathstage:ch1`)) %>%
  dplyr::mutate(T_substage_pathological = stringr::str_extract(uncurated$`pathstage:ch1`,
                                                                "[a-c]+")) %>% 
  dplyr::mutate(T_clinical = readr::parse_number(uncurated$`clint_stage:ch1`)) %>% 
  dplyr::mutate(T_substage_clinical = stringr::str_extract(uncurated$`clint_stage:ch1`,
                                                           "[a-c]+")) %>%
  dplyr::mutate(metastasis_occurrence_status = dplyr::case_when(
    uncurated$`metsevent:ch1` == "no" ~ 0,
    uncurated$`metsevent:ch1` == "yes" ~ 1
  )) %>%
  dplyr::mutate(days_to_metastatic_occurrence = as.numeric(
    uncurated$`metsfreetime_months:ch1`
    )*30.5) %>%
  dplyr::mutate(psa = as.numeric(uncurated$`pretxpsa:ch1`)) %>%
  dplyr::mutate(extraprostatic_extension = dplyr::case_when(
    uncurated$`ece_binary:ch1` == "No" ~ 0,
    uncurated$`ece_binary:ch1` == "Yes" ~ 1
  )) %>% 
  dplyr::mutate(seminal_vesicle_invasion= case_when(
    uncurated$`svi:ch1` == "Negative" ~ 0,
    uncurated$`svi:ch1` == "Positive" ~ 1
  )) 
  
  
clinical_hieronymus <- curated

save(clinical_hieronymus, file = "data-raw/clinical_hieronymus.RData")

