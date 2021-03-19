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
  ## TDL: Below cases might've changed in cBio as they were no longer correct (prefixed with {"0:","1:"})
  #dplyr::mutate(overall_survival_status = dplyr::case_when(
  #  is.na(overall_survival_status) ~ NA_real_,
  #  overall_survival_status == "DECEASED" ~ 1,
  #  overall_survival_status != "DECEASED" ~ 0
  #)) %>% 
  dplyr::mutate(overall_survival_status = dplyr::case_when(
    is.na(overall_survival_status) ~ NA_real_,
    overall_survival_status == "1:DECEASED" ~ 1,
    overall_survival_status == "0:LIVING" ~ 0
  )) %>% 
  dplyr::mutate(days_to_overall_survival = as.numeric(uncurated$OS_MONTHS) * 30.5) %>%
  dplyr::mutate(disease_specific_recurrence_status = uncurated$DFS_STATUS) %>% 
  ## TDL: Below cases might've changed in cBio as they were no longer correct (prefixed with {"0:","1:"})
  #dplyr::mutate(disease_specific_recurrence_status = dplyr::case_when(
  #  disease_specific_recurrence_status == "[Not Available]" ~ NA_real_,
  #  disease_specific_recurrence_status == "Recurred/Progressed" ~ 1,
  #  disease_specific_recurrence_status == "DiseaseFree" ~ 0
  #)) %>% 
  dplyr::mutate(disease_specific_recurrence_status = dplyr::case_when(
    disease_specific_recurrence_status == "" ~ NA_real_,
    disease_specific_recurrence_status == "1:Recurred/Progressed" ~ 1,
    disease_specific_recurrence_status == "0:DiseaseFree" ~ 0
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
  )) 

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

## TDL: below script is saved only temporarily as it overlaps with Jordan's solution
if(FALSE){
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
}

uncurated <- Biobase::pData(gse[[2]])
gse_gex <- GEOquery::getGEO("GSE21034", GSEMatrix = TRUE)
uncurated_gex_exon <- Biobase::pData(gse_gex[[2]])
uncurated_gex_transcript <- Biobase::pData(gse_gex[[1]])

gse_cna <- GEOquery::getGEO("GSE21035", GSEMatrix = TRUE)
uncurated_cna <- Biobase::pData(gse_cna[[1]])

gse_mrna <- GEOquery::getGEO("GSE21036", GSEMatrix = TRUE)
uncurated_mrna <- Biobase::pData(gse_mrna[[1]])

mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
uncurated_cbio <- cgdsr::getClinicalData(mycgds, caseList = "prad_mskcc_all")

uncurated <- uncurated_cna %>% 
  dplyr::select(title, cna = geo_accession, dplyr::everything()) %>% 
  dplyr::mutate(title = stringr::str_remove(title, "Prostate tumor ")) %>% 
  dplyr::mutate(title = stringr::str_remove(title, " \\(aCGH\\)")) %>% 
  dplyr::left_join(uncurated_gex_exon %>% 
                     dplyr::select(title, gex_exon = geo_accession) %>%
                     dplyr::mutate(title = stringr::str_remove(title, "Prostate tumor ")) %>% 
                     dplyr::mutate(title = stringr::str_remove(title, " exon")),
                   by = "title") %>% 
  dplyr::left_join(uncurated_gex_transcript %>% 
                     dplyr::select(title, gex_transcript = geo_accession) %>%
                     dplyr::mutate(title = stringr::str_remove(title, "Prostate tumor ")) %>% 
                     dplyr::mutate(title = stringr::str_remove(title, " transcript")),
                   by = "title") %>% 
  dplyr::left_join(uncurated_mrna %>% 
                     dplyr::select(title, mrna = geo_accession) %>%
                     dplyr::mutate(title = stringr::str_remove(title, "Prostate tumor ")) %>% 
                     dplyr::mutate(title = stringr::str_remove(title, " miRNA")), 
                   by = "title") %>%
  dplyr::mutate(geo_accession = paste(paste0("cna: ", cna),
                                      paste0("gex_exon: ", gex_exon), 
                                      paste0("gex_transcript: ", gex_transcript), 
                                      paste0("mrna: ", mrna), sep = "|"))

rownames(uncurated) <- uncurated$title

curated <- initial_curated_df(
  df_rownames = rownames(uncurated),
  template_name="data-raw/template_prad.csv")

curated <- curated %>% 
  dplyr::mutate(study_name = "Taylor, et al.") %>% 
  dplyr::mutate(sample_name = uncurated$geo_accession) %>% 
  dplyr::mutate(patient_id = row.names(uncurated)) %>%
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

uncurated_cbio <- uncurated_cbio[match(curated$patient_id, row.names(uncurated_cbio)),]

curated <- curated %>%
  dplyr::mutate(gleason_major = uncurated_cbio$GLEASON_SCORE_1) %>% 
  dplyr::mutate(gleason_minor = uncurated_cbio$GLEASON_SCORE_2) %>% 
  dplyr::mutate(grade_group = dplyr::case_when(
    gleason_grade == 6 ~ "<=6",
    gleason_grade %in% 8:10 ~ ">=8",
    gleason_major == 3 & gleason_minor == 4 ~ "3+4",
    gleason_major == 4 & gleason_minor == 3 ~ "4+3",
  )) %>% 
  dplyr::mutate(ERG_fusion_GEX = uncurated_cbio$ERG_FUSION_GEX) %>% 
  dplyr::mutate(ERG_fusion_GEX = dplyr::case_when(
    ERG_fusion_GEX == "Negative" ~ 0,
    ERG_fusion_GEX == "Positive" ~ 1,
    TRUE ~ NA_real_
  )) %>% 
  dplyr::mutate(ERG_fusion_CNA = uncurated_cbio$ERG_FUSION_ACGH) %>% 
  dplyr::mutate(ERG_fusion_CNA =  dplyr::case_when(
    ERG_fusion_CNA == "Positive" ~ 1,
    ERG_fusion_CNA %in% c("Negative", "Flat") ~ 0,
    TRUE ~ NA_real_
  )) %>% 
  dplyr::mutate(disease_specific_recurrence_status = uncurated_cbio$DFS_STATUS) %>% 
  dplyr::mutate(disease_specific_recurrence_status = dplyr::case_when(
    disease_specific_recurrence_status == "Recurred" ~ 1,
    disease_specific_recurrence_status == "DiseaseFree" ~ 0,
    TRUE ~ NA_real_
  )) %>% 
  dplyr::mutate(days_to_disease_specific_recurrence == uncurated_cbio$DFS_MONTHS) 

curated <- curated %>%
  dplyr::filter(sample_type %in% c("metastasis", "primary"))

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


###############################################################################
# 
# ICGC clinical data curation scripts
#
###############################################################################

###
## CA (Canadian) ICGC dataset
###

# Generation script from generate.R, same function as is used for downloading and transforming the omics
uncurated <- curatedPCaData:::generate_icgc("PRAD_CA", "clinical")

# Format empty data frames according to the prad template
curated <- initial_curated_df(
  df_rownames = rownames(uncurated),
  template_name="data-raw/template_prad.csv")

# Mimic previous curation piping
curated <- curated %>% 
  dplyr::mutate(study_name = "ICGC_CA") %>% 
  dplyr::mutate(sample_name = uncurated$icgc_sample_id) %>% 
  dplyr::mutate(patient_id = uncurated$icgc_donor_id) %>%
  dplyr::mutate(overall_survival_status = dplyr::case_when(
    uncurated$donor_vital_status == "alive" ~ 0,
    uncurated$donor_vital_status == "deceased" ~ 1,
    uncurated$donor_vital_status == "" ~ NA_real_,
    is.na(uncurated$donor_vital_status) ~ NA_real_
  )) %>%
  dplyr::mutate(age_at_initial_diagnosis = uncurated$donor_age_at_diagnosis) %>%
  dplyr::mutate(days_to_overall_survival = (uncurated$donor_age_at_last_followup - uncurated$donor_age_at_diagnosis)*365) %>%
  dplyr::mutate(gleason_major = as.integer(stringr::str_sub(uncurated$tumour_grade,1,1))) %>%
  dplyr::mutate(gleason_minor = as.integer(stringr::str_sub(uncurated$tumour_grade,3,3))) %>%
  dplyr::mutate(grade_group = dplyr::case_when(
      stringr::str_sub(uncurated$tumour_grade,1,3) == "3+3" ~ "<=6",
      stringr::str_sub(uncurated$tumour_grade,1,3) == "3+4" ~ "3+4",
      stringr::str_sub(uncurated$tumour_grade,1,3) == "4+3" ~ "4+3",
      stringr::str_sub(uncurated$tumour_grade,1,3) %in% c("4+4", "4+5") ~ ">=8"
  )) %>%
  dplyr::mutate(gleason_grade = gleason_minor + gleason_major) %>%  
  ## TODO: Double-check: clinical or pathological?
  dplyr::mutate(T_clinical = readr::parse_number(uncurated$tumour_stage)) %>%   
  dplyr::mutate(T_substage_clinical = stringr::str_extract(uncurated$tumour_stage, "[a-c]+")) 
  
clinical_icgcca <- curated

save(clinical_icgcca, file = "data-raw/clinical_icgcca.RData")

###
## FR (French) ICGC dataset
###

# Generation script from generate.R, same function as is used for downloading and transforming the omics
uncurated <- curatedPCaData:::generate_icgc("PRAD_FR", "clinical")

# Format empty data frames according to the prad template
curated <- initial_curated_df(
  df_rownames = rownames(uncurated),
  template_name="data-raw/template_prad.csv")

# Mimic previous curation piping
curated <- curated %>% 
  dplyr::mutate(study_name = "ICGC_FR") %>% 
  dplyr::mutate(sample_name = uncurated$icgc_sample_id) %>% 
  dplyr::mutate(patient_id = uncurated$icgc_donor_id) %>%
  dplyr::mutate(age_at_initial_diagnosis = uncurated$donor_age_at_diagnosis) %>%
  dplyr::mutate(days_to_overall_survival = uncurated$donor_survival_time) %>%
  dplyr::mutate(overall_survival_status = dplyr::case_when(
    uncurated$donor_vital_status == "alive" ~ 0,
    uncurated$donor_vital_status == "deceased" ~ 1,
    is.na(uncurated$donor_vital_status) ~ NA_real_
  )) 
  # Source claims to use GLEASON_SCORE and PT_STAGE fields, yet the values in corresponding columns are NA?
  # TODO: Verify if these are at some other table 
   
clinical_icgcfr <- curated

save(clinical_icgcfr, file = "data-raw/clinical_icgcfr.RData")

###
## UK (United Kingdom) ICGC dataset
###

# Generation script from generate.R, same function as is used for downloading and transforming the omics
uncurated <- curatedPCaData:::generate_icgc("PRAD_UK", "clinical")

# Format empty data frames according to the prad template
curated <- initial_curated_df(
  df_rownames = rownames(uncurated),
  template_name="data-raw/template_prad.csv")

# Mimic previous curation piping
curated <- curated %>% 
  dplyr::mutate(study_name = "ICGC_UK") %>% 
  dplyr::mutate(sample_name = uncurated$icgc_sample_id) %>% 
  dplyr::mutate(patient_id = uncurated$icgc_donor_id) %>%
  dplyr::mutate(age_at_initial_diagnosis = uncurated$donor_age_at_diagnosis) %>%
  dplyr::mutate(days_to_overall_survival = uncurated$donor_survival_time) %>%
  dplyr::mutate(overall_survival_status = dplyr::case_when(
    uncurated$donor_vital_status == "alive" ~ 0,
    uncurated$donor_vital_status == "deceased" ~ 1,
    uncurated$donor_vital_status == "" ~ NA_real_,
    is.na(uncurated$donor_vital_status) ~ NA_real_
  )) %>%
  dplyr::mutate(gleason_major = as.integer(stringr::str_sub(uncurated$tumour_grade,1,1))) %>%
  dplyr::mutate(gleason_minor = as.integer(stringr::str_sub(uncurated$tumour_grade,3,3))) %>%
  dplyr::mutate(grade_group = dplyr::case_when(
      stringr::str_sub(uncurated$tumour_grade,1,3) == "3+3" ~ "<=6",
      stringr::str_sub(uncurated$tumour_grade,1,3) == "3+4" ~ "3+4",
      stringr::str_sub(uncurated$tumour_grade,1,3) == "4+3" ~ "4+3",
      stringr::str_sub(uncurated$tumour_grade,1,3) %in% c("4+4", "4+5") ~ ">=8"
  )) %>%
  dplyr::mutate(gleason_grade = gleason_minor + gleason_major) %>%  
  dplyr::mutate(metastasis_occurrence_status = dplyr::case_when(
    stringr::str_sub(uncurated$tumour_stage,6,7) == "M1" ~ 1,
    stringr::str_sub(uncurated$tumour_stage,6,7) == "M0" ~ 0,
    stringr::str_sub(uncurated$tumour_stage,6,7) %in% c("Mx", "") ~ NA_real_
  )) %>%
  dplyr::mutate(M_stage = stringr::str_sub(uncurated$tumour_stage,6,7))
  # clinical -> uncurated$donor_tumour_stage_at_diagnosis
  # pathology -> uncurated$tumour_stage
  ## contains a lot of information for the TNM
  #> table(uncurated$donor_tumour_stage_at_diagnosis)
  #
  #         pT2c pN0  pT3 pN0 pT3b pN1  pT4 pN1  T1bN0Mx  T1cN0M0  T1cN0Mx 
  #       5        5       10        6        4        2        4       18 
  # T1cNxM1  T1cNxMx    T2 Nx  T2aN0Mx  T2aNxMx  T2b pN1  T2bN0M0  T2bN0Mx 
  #       2      164        3       45       20        6        2       14 
  # T2bNxMx  T2cN0M0  T2cN0M1  T2cN0Mx  T2cN1Mx  T2cNxM0  T2cNxMx   T2N0Mx 
  #       8        8        2       50        2        2        8        2 
  #  T2NxMx    T3 Nx  T3aN0M0  T3aN0Mx  T3aNxMx  T3bN0Mx  T3bNxMx   T3N0Mx 
  #       6        5        2       20       23        2        6        2 
  #  T3NxM0    T4 Nx   TxNxMx 
  #       2       17       32 
  #> table(uncurated$tumour_stage)
  #
  #        T2aNxMx T2cN0M0 T2cN0M1 T2cN0Mx T2cNxMx T3aN0M0 T3aN0Mx T3aN1Mx T3aNxMx 
  #    273       6       8       2      20      44       6      65       2      37 
  #T3bN0Mx T3bN1M0 T3bN1Mx T3bNxMx T4xNxMx TxxNxMx 
  #   12       2      12       4       2      14
  #
  
clinical_icgcuk <- curated

save(clinical_icgcuk, file = "data-raw/clinical_icgcuk.RData")


#############################################################################
#
#
# Friedrich et at German samples.  This code downloads both gene expression 
# and clinical data -- just the curation of the clinical data is
# provided here 
#
#
##############################################################################

gset <- GEOquery::getGEO("GSE134051", GSEMatrix =TRUE, getGPL=TRUE)

uncurated <- Biobase::pData(gset[[1]]) 
curated <- initial_curated_df(
  df_rownames = rownames(uncurated),
  template_name="./template_prad.csv")

curated <- curated %>% 
  dplyr::mutate(study_name = "Friedrich et al.") %>%
  dplyr::mutate(sample_name = row.names(uncurated)) %>% 
  dplyr::mutate(patient_id = row.names(uncurated)) %>%
  dplyr::mutate(alt_sample_name = unlist(lapply(stringr::str_split(uncurated$title, '_'), function(x) x[3]))) %>%
  dplyr::mutate(gleason_grade = dplyr::case_when(uncurated$'gleason score:ch1' == 'None' ~ NA_character_, 
                                                 TRUE ~ uncurated$'gleason score:ch1')) %>%
  dplyr::mutate(gleason_grade = as.numeric(gleason_grade)) %>% 
  dplyr::mutate(race = 'caucasian') %>% 
  dplyr::mutate(tissue_source = 'prostatectomy') %>%
  dplyr::mutate(sample_type = dplyr::case_when(
                                               uncurated$'risk group:ch1' == 'C' ~ 'BPH', 
                                               uncurated$'risk group:ch1' == 'V' ~ 'primary',
                                               uncurated$'risk group:ch1' == 'Ms' ~ 'primary',
                                               uncurated$'risk group:ch1' == 'Md' ~ 'primary',
                                               uncurated$'risk group:ch1' == 'L' ~ 'primary',
                                               uncurated$'risk group:ch1' == 'H+st' ~ 'primary',
                                               uncurated$'risk group:ch1' == 'H+dt' ~ 'primary',
                                               uncurated$'risk group:ch1' == 'H-st' ~ 'primary',
                                               uncurated$'risk group:ch1' == 'H-dt' ~ 'primary',
                                               uncurated$'risk group:ch1' == 'H+sf' ~ 'adjacentnormal',
                                               uncurated$'risk group:ch1' == 'H+df' ~ 'adjacentnormal',
                                               uncurated$'risk group:ch1' == 'H-sf' ~ 'adjacentnormal',
                                               uncurated$'risk group:ch1' == 'H-df' ~ 'adjacentnormal'
                                               )) %>%
  dplyr::mutate(days_to_overall_survival = ceiling(as.numeric(uncurated$'follow-up time in months:ch1')*30.5)) %>%
  dplyr::mutate(frozen_ffpe = 'frozen') %>%
  dplyr::mutate(microdissected = 1) %>%
  dplyr::mutate(overall_survival_status = dplyr::case_when(
                                                           uncurated$'risk group:ch1' == 'C' ~ 0, 
                                                           uncurated$'risk group:ch1' == 'V' ~ 0,
                                                           uncurated$'risk group:ch1' == 'Ms' ~ 0,
                                                           uncurated$'risk group:ch1' == 'Md' ~ 1,
                                                           uncurated$'risk group:ch1' == 'L' ~ 0,
                                                           uncurated$'risk group:ch1' == 'H+st' ~ 0,
                                                           uncurated$'risk group:ch1' == 'H+dt' ~ 1,
                                                           uncurated$'risk group:ch1' == 'H-st' ~ 0,
                                                           uncurated$'risk group:ch1' == 'H-dt' ~ 1,
                                                           uncurated$'risk group:ch1' == 'H+sf' ~ 0,
                                                           uncurated$'risk group:ch1' == 'H+df' ~ 0,
                                                           uncurated$'risk group:ch1' == 'H-sf' ~ 0,
                                                           uncurated$'risk group:ch1' == 'H-df' ~ 0
                                                           )) 



clinical_friedrich <- curated

save(clinical_friedrich, file = "./clinical_friedrich.RData")

##################################################################
##################################################################
#
#
# Wallace et al. samples: This code downloads both gene expression
# BOTH clinical and GPL, but just the clinical data is used
#
#
##################################################################
##################################################################

library(GEOquery)
#library(limma)
#library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE6956", GSEMatrix =TRUE, getGPL=FALSE)


# clinical

library(magrittr)
library(dplyr)


uncurated <- Biobase::pData(gset[[1]])

unmatched_healty_tissue = c('GSM160418', 'GSM160419', 'GSM160420', 'GSM160421', 'GSM160422', 'GSM160430')

uncurated = uncurated[!is.element(uncurated$geo_accession, unmatched_healty_tissue), ]

curated <- initial_curated_df(
  df_rownames = rownames(uncurated),
  template_name="./template_prad.csv")

curated <- curated %>%
  dplyr::mutate(study_name = "Wallace et al.") %>%
  dplyr::mutate(sample_name = row.names(uncurated)) %>%
  dplyr::mutate(patient_id = stringr::str_extract(uncurated$title, 'patient \\d+')) %>%
  dplyr::mutate(alt_sample_name = stringr::str_extract(uncurated$title, 'patient \\d+')) %>%
  dplyr::mutate(gleason_grade = dplyr::case_when(
                                                 uncurated$'gleason sum:ch1' == 'NA' ~ NA_character_,
                                                 TRUE ~ uncurated$'gleason sum:ch1')) %>%
  dplyr::mutate(gleason_grade = as.numeric(gleason_grade)) %>%
  dplyr::mutate(race = dplyr::case_when(
                                        uncurated$"race:ch1" == 'African American' ~ 'african_american',
                                        uncurated$"race:ch1" == 'Caucasian' ~ 'caucasian',
                                        uncurated$"race:ch1" == 'NA' ~ NA_character_
                                        )) %>%
  dplyr::mutate(smoking_status = as.numeric(dplyr::case_when(
                                                 uncurated$"smoking status:ch1" == 'Current' ~ '1',
                                                 uncurated$"smoking status:ch1" == 'Past' ~ '1',
                                                 uncurated$"smoking status:ch1" == 'Never' ~ '0',
                                                 uncurated$"smoking status:ch1" == 'Unknown' ~ NA_character_,
                                                 uncurated$"smoking status:ch1" == 'NA' ~ NA_character_
                                                 ))) %>%
  dplyr::mutate(angiolymphatic_invasion = as.numeric(dplyr::case_when(
                                                           uncurated$"angio lymphatic invasion:ch1" == 'Yes' ~ '1',
                                                           uncurated$"angio lymphatic invasion:ch1" == 'No' ~ '0',
                                                           uncurated$"angio lymphatic invasion:ch1" == 'NA' ~ NA_character_
                                                           ))) %>%
  dplyr::mutate(seminal_vesicle_invasion = as.numeric(dplyr::case_when(
                                                           uncurated$"seminal vesicle invasion:ch1" == 'Yes' ~ '1',
                                                           uncurated$"seminal vesicle invasion:ch1" == 'No' ~ '0',
                                                           uncurated$"seminal vesicle invasion:ch1" == 'NA' ~ NA_character_
                                                           ))) %>%
  dplyr::mutate(perineural_invasion = as.numeric(dplyr::case_when(
                                                           uncurated$"perineural invasion:ch1" == 'Yes' ~ '1',
                                                           uncurated$"perineural invasion:ch1" == 'No' ~ '0',
                                                           uncurated$"perineural invasion:ch1" == 'NA' ~ NA_character_
                                                           ))) %>%
  dplyr::mutate(extraprostatic_extension = as.numeric(dplyr::case_when(
                                                           uncurated$"extraprostatic extension:ch1" == 'Focal' ~ '1',
                                                           uncurated$"extraprostatic extension:ch1" == 'Multifocal' ~ '1',
                                                           uncurated$"extraprostatic extension:ch1" == 'Established' ~ '1',
                                                           uncurated$"extraprostatic extension:ch1" == 'None' ~ '0',
                                                           uncurated$"perineural invasion:ch1" == 'NA' ~ NA_character_
                                                           ))) %>%
  dplyr::mutate(tumor_margins_positive = as.numeric(dplyr::case_when(
                                                           uncurated$"are surgical margins involved:ch1" == 'Tumor focal at margin' ~ '1',
                                                           uncurated$"are surgical margins involved:ch1" == 'Tumor widespread a surgical margins' ~ '1',
                                                           uncurated$"are surgical margins involved:ch1" == 'Tumor widespread at margin' ~ '1',
                                                           uncurated$"are surgical margins involved:ch1" == 'All surgical margins are free of tumor' ~ '0',
                                                           uncurated$"are surgical margins involved:ch1" == 'Unknown' ~ NA_character_,
                                                           uncurated$"are surgical margins involved:ch1" == 'NA' ~ NA_character_
                                                           ))) %>%
  dplyr::mutate(sample_type = dplyr::case_when(
                                               uncurated$source_name_ch1 == "Adenocarcinoma (NOS) of the prostate" ~ 'primary',
                                               uncurated$source_name_ch1 == "Normal prostate" ~ 'adjacentnormal'
                                               ))  %>%
  dplyr::mutate(frozen_ffpe = 'frozen') %>%
  dplyr::mutate(microdissected = 0) %>%
  dplyr::mutate(therapy_radiation_initial = 0) %>%
  dplyr::mutate(therapy_radiation_salvage = 0) %>%
  dplyr::mutate(therapy_surgery_initial = 0) %>%
  dplyr::mutate(therapy_hormonal_initial = 0)




clinical_wallace <- curated

save(clinical_wallace, file = "./clinical_wallace.RData")


###
#
# Chandran et al., BMC Cancer 2007
# Yu et al. J Clin Oncol 2004
#
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6919
# 
# Includes normal samples, tumor samples, as well as metastatic samples
#
###

gse <- GEOquery::getGEO("GSE6919", GSEMatrix = TRUE)

uncurated1 <- Biobase::pData(gse[[1]]) # Batch 1
uncurated2 <- Biobase::pData(gse[[2]]) # Batch 2
uncurated3 <- Biobase::pData(gse[[3]]) # Batch 3
# Bind metadata over the 3 GPLs
uncurated <- rbind(uncurated1, uncurated2, uncurated3)

# Base curation metadat
curated <- initial_curated_df(
  df_rownames = rownames(uncurated),
  template_name="./template_prad.csv")

# Pipe through the available fields
curated <- curated %>% 
	dplyr::mutate(study_name = "Chandran et al.") %>%
	dplyr::mutate(patient_id = uncurated$'title') %>%
	dplyr::mutate(sample_name = row.names(uncurated)) %>% 
	#dplyr::mutate(alt_sample_name = uncurated$'title') %>%
	dplyr::mutate(age_at_initial_diagnosis = uncurated$'Age:ch1') %>%
	dplyr::mutate(race = dplyr::case_when(
		uncurated$'Race:ch1' == 'Caucasian' ~ 'caucasian', 
		uncurated$'Race:ch1' == 'African American' ~ 'african_american'
	)) %>%
	dplyr::mutate(gleason_grade = uncurated$'Gleason Grade:ch1') %>%
	# TODO: Double-check if the T were determined at diagnosis or post-surgery
	# T_pathological & T_substage_pathological OR T_clinical & T_substage_clinical
	dplyr::mutate(T_pathological = dplyr::case_when(
		uncurated$'Tumor stage:ch1' %in% c('T2a','T2b') ~ 2, 
		uncurated$'Tumor stage:ch1' %in% c('T3a','T3b') ~ 3,
		uncurated$'Tumor stage:ch1' %in% c('T4', 'T4a') ~ 4,
		is.na(uncurated$'Tumor stage:ch1') ~ NA_real_,
	)) %>%
	dplyr::mutate(T_substage_pathological = dplyr::case_when(
		uncurated$'Tumor stage:ch1' %in% c('T2a','T3a', 'T4a') ~ 'a', 
		uncurated$'Tumor stage:ch1' %in% c('T2b','T3b') ~ 'b',
		uncurated$'Tumor stage:ch1' == "T4" ~ NA_character_,
		is.na(uncurated$'Tumor stage:ch1') ~ NA_character_,
	)) %>%
	dplyr::mutate(sample_type = dplyr::case_when(
		uncurated$'Tissue:ch1' %in% c('primary prostate tumor') ~ 'primary', 
		uncurated$'Tissue:ch1' %in% c('normal prostate tissue adjacent to tumor') ~ 'adjacentnormal', 
		uncurated$'Tissue:ch1' %in% c('normal prostate tissue free of any pathological alteration from brain-dead organ donor') ~ 'healthy', 
		uncurated$'Tissue:ch1' %in% c('metastases recurrent in prostate','prostate tumor metastases in adrenal gland','prostate tumor metastases in kidney','prostate tumor metastases in left inguinal lymph node','prostate tumor metastases in liver','prostate tumor metastases in lung','prostate tumor metastases in para aortic lymph node','prostate tumor metastases in para tracheal lymph node','prostate tumor metastases in paratracheal lymph node','prostate tumor metastases in retroperitoneal lymph node') ~ 'metastatic', 
		is.na(uncurated$'Tissue:ch1') ~ NA_character_,
	)) %>%
	dplyr::mutate(metastatic_site = dplyr::case_when(
		uncurated$'Tissue:ch1' %in% c('metastases recurrent in prostate') ~ 'prostate', 
		uncurated$'Tissue:ch1' %in% c('prostate tumor metastases in adrenal gland') ~ 'adrenal_gland',
		uncurated$'Tissue:ch1' %in% c('prostate tumor metastases in kidney') ~ 'kidney',
		uncurated$'Tissue:ch1' %in% c('prostate tumor metastases in liver') ~ 'liver',
		uncurated$'Tissue:ch1' %in% c('prostate tumor metastases in lung') ~ 'lung',
		uncurated$'Tissue:ch1' %in% c('prostate tumor metastases in left inguinal lymph node','prostate tumor metastases in para aortic lymph node','prostate tumor metastases in para tracheal lymph node','prostate tumor metastases in paratracheal lymph node','prostate tumor metastases in retroperitoneal lymph node') ~ 'lymph_node',
		uncurated$'Tissue:ch1' %in% c('primary prostate tumor','normal prostate tissue adjacent to tumor','normal prostate tissue free of any pathological alteration from brain-dead organ donor') ~ NA_character_, 
		is.na(uncurated$'Tissue:ch1') ~ NA_character_,
	))

rownames(curated) <- curated$sample_name
clinical_chandran <- curated

save(clinical_chandran, file = "./clinical_chandran.RData")

######################################################################
#cBioportal Barbieri Broad/Cornell Data
#####################################################################

mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
uncurated <- cgdsr::getClinicalData(mycgds, caseList="prad_broad_sequenced")
#mycancerstudy = cgdsr::getCancerStudies(mycgds)
#mycaselist = cgdsr::getCaseLists(mycgds,"prad_broad")

# create the curated object
curated <- initial_curated_df(
  df_rownames = rownames(uncurated),
  template_name="data-raw/template_prad.csv")

curated <- curated %>% 
  dplyr::mutate(study_name = "Barbieri") %>%
  dplyr::mutate(sample_name = row.names(uncurated)) %>%
  dplyr::mutate(patient_id = row.names(uncurated)) %>%
  dplyr::mutate(gleason_grade = uncurated$GLEASON_SCORE) %>%
  dplyr::mutate(gleason_major = as.integer(stringr::str_sub(uncurated$GLEASON_SCORE,1,1))) %>%
  dplyr::mutate(gleason_minor = as.integer(stringr::str_sub(uncurated$GLEASON_SCORE,3,3))) %>%
  dplyr::mutate(grade_group = dplyr::case_when(
    stringr::str_sub(uncurated$GLEASON_SCORE,1,3) == "3+3" ~ "<=6",
    stringr::str_sub(uncurated$GLEASON_SCORE,1,3) == "3+4" ~ "3+4",
    stringr::str_sub(uncurated$GLEASON_SCORE,1,3) == "4+3" ~ "4+3",
    stringr::str_sub(uncurated$GLEASON_SCORE,1,3) %in% c("4+4", "4+5") ~ ">=8"
  )) %>%
  dplyr::mutate(psa = uncurated$SERUM_PSA) %>% 
  dplyr::mutate(age_at_initial_diagnosis = uncurated$AGE) %>%
  dplyr::mutate(T_clinical = readr::parse_number(uncurated$TUMOR_STAGE)) %>% 
  dplyr::mutate(T_substage_clinical = stringr::str_extract(uncurated$TUMOR_STAGE, "[a-c]+")) %>%
  dplyr::mutate(TMPRSS2_ERG_FUSION_STATUS=uncurated$TMPRSS2_ERG_FUSION_STATUS) %>%
  dplyr::mutate(GLEASON_SCORE_PERCENT_4_AND_5=uncurated$GLEASON_SCORE_PERCENT_4_AND_5)

clinical_barbieri <- curated
save(clinical_barbieri, file = "data-raw/clinical_barbieri.RData")

#####################################################################################
#cBioportal Ren eururol 2017 Data
#####################################################################################

mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
uncurated <- cgdsr::getClinicalData(mycgds, caseList="prad_eururol_2017_sequenced")
#mycancerstudy = cgdsr::getCancerStudies(mycgds)
#mycaselist = cgdsr::getCaseLists(mycgds,"prad_eururol_2017")

# create the curated object
curated <- initial_curated_df(
  df_rownames = rownames(uncurated),
  template_name="data-raw/template_prad.csv")

curated <- curated %>% 
  dplyr::mutate(study_name = "Ren") %>%
  dplyr::mutate(sample_name = row.names(uncurated)) %>%
  dplyr::mutate(patient_id = row.names(uncurated)) %>%
  dplyr::mutate(psa = uncurated$PSA) %>%
  dplyr::mutate(gleason_grade = uncurated$GLEASON_SCORE) %>%
  dplyr::mutate(gleason_major = as.integer(stringr::str_sub(uncurated$GLEASON_SCORE,1,1))) %>%
  dplyr::mutate(gleason_minor = as.integer(stringr::str_sub(uncurated$GLEASON_SCORE,3,3))) %>%
  dplyr::mutate(grade_group = dplyr::case_when(
    stringr::str_sub(uncurated$GLEASON_SCORE,1,3) == "3+3" ~ "<=6",
    stringr::str_sub(uncurated$GLEASON_SCORE,1,3) == "3+4" ~ "3+4",
    stringr::str_sub(uncurated$GLEASON_SCORE,1,3) == "4+3" ~ "4+3",
    stringr::str_sub(uncurated$GLEASON_SCORE,1,3) %in% c("4+4", "4+5") ~ ">=8"
  )) %>%
  dplyr::mutate(age_at_initial_diagnosis = uncurated$AGE) %>%
  dplyr::mutate(T_clinical = readr::parse_number(uncurated$TNMSTAGE)) %>% 
  dplyr::mutate(T_substage_clinical = stringr::str_extract(uncurated$TNMSTAGE, "[a-c]+")) %>%
  dplyr::mutate(FPSA_PSA = uncurated$FPSA_PSA) %>%
  dplyr::mutate(Tumor_Purity = uncurated$TUMOR_PURITY) %>%
  dplyr::mutate(BLADDER_NECK_INVASION = uncurated$BLADDER_NECK_INVASION) %>%
  dplyr::mutate(SEMINAL_VESICLE_INVASION = uncurated$SEMINAL_VESICLE_INVASION) %>%
  dplyr::mutate(Fraction_genome_altered=uncurated$FRACTION_GENOME_ALTERED) %>%
  dplyr::mutate(LYMPH_NODE_METASTASIS=uncurated$LYMPH_NODE_METASTASIS) %>%
  dplyr::mutate(MUTATION_COUNT=uncurated$MUTATION_COUNT)


clinical_ren <- curated
save(clinical_ren, file = "data-raw/clinical_ren.RData")


#######################################################################
#Kim et al
######################################################################

gse <- GEOquery::getGEO("GSE119616", GSEMatrix = TRUE)

uncurated <- Biobase::pData(gse[[1]])

curated <- initial_curated_df(
  df_rownames = rownames(uncurated),
  template_name="data-raw/template_prad.csv")

curated <- curated %>% 
  dplyr::mutate(study_name = "Kim et al.") %>%
  dplyr::mutate(sample_name = row.names(uncurated)) %>% 
  dplyr::mutate(patient_id = stringr::str_remove(uncurated$title,
                                                      "prostate_cancer_biopsy_sample_pid_")) %>%
  dplyr::mutate(tissue_source = stringr::str_remove(uncurated$characteristics_ch1.10,"tissue:")) %>%
  dplyr::mutate(age_at_initial_diagnosis = stringr::str_remove(uncurated$characteristics_ch1.2,
                                                 "age:")) %>%
  dplyr::mutate(psa = stringr::str_remove(uncurated$characteristics_ch1.6,
                                                               "psa:")) %>%
  dplyr::mutate(gleason_major = stringr::str_remove(uncurated$characteristics_ch1.8,"primary gleason grade:")) %>%
  dplyr::mutate(gleason_minor = stringr::str_remove(uncurated$characteristics_ch1.9,"secondary gleason grade:")) %>%
  dplyr::mutate(gleason_grade = gsub(" ", "", paste(gleason_major,"+",gleason_minor))) %>%
  dplyr::mutate(grade_group = dplyr::case_when(
    gleason_grade == "3+3" ~ "<=6",
    gleason_grade == "3+4" ~ "3+4",
    gleason_grade == "4+3" ~ "4+3",
    gleason_grade %in% c("4+4", "4+5") ~ ">=8"
  )) %>%
  dplyr::mutate(T_clinical = readr::parse_number(uncurated$`tumor stage:ch1`)) %>% 
  dplyr::mutate(T_substage_clinical = stringr::str_extract(uncurated$`tumor stage:ch1`, "[a-c]+"))
  
clinical_kim <- curated

save(clinical_kim, file = "data-raw/clinical_kim.RData")
  

####
#
# Barwick et al.
#
####

gse <- GEOquery::getGEO("GSE18655", GSEMatrix = TRUE)

uncurated <- Biobase::pData(gse[[1]])

curated <- initial_curated_df(
  df_rownames = rownames(uncurated),
  template_name="data-raw/template_prad.csv")

curated <- curated %>% 
  dplyr::mutate(study_name = "Barwick et al.") %>%
  dplyr::mutate(patient_id = paste0("X", uncurated$'title')) %>%
  dplyr::mutate(age_at_initial_diagnosis = as.numeric(uncurated$'age:ch1')) %>%
  dplyr::mutate(psa = as.numeric(uncurated$'psa:ch1')) %>%
  dplyr::mutate(gleason_grade = as.numeric(uncurated$'gleason score:ch1')) %>%
  dplyr::mutate(days_to_disease_specific_recurrence = round(30.5*as.numeric(uncurated$'follow-up (months):ch1'),0)) %>%
  dplyr::mutate(disease_specific_recurrence_status = as.numeric(uncurated$'recurrence:ch1' == "Rec")) %>%
  dplyr::mutate(tumor_margins_positive = as.numeric(uncurated$'positive surgical margin:ch1' == "Positive Margin")) %>%
  dplyr::mutate(tissue_source = "prostatectomy") %>% 
  # Based on the publication text, ERG-fusion status was determined using over-expression of ERG-transcripts in gene expression and then experimentally validated
  dplyr::mutate(ERG_fusion_GEX = as.numeric(uncurated$'tmprss2:ch1' == "ERG Fusion: Fusion")) %>%
  # It appears the T staging was based on pathology
  dplyr::mutate(T_pathological = as.numeric(uncurated$'grade:ch1'))
  
clinical_barwick <- curated

save(clinical_barwick, file = "./data-raw/clinical_barwick.RData")

# Clinical mapping of Barwick replicates? The raw GEX file of Barwick contains 180 samples, but is still not equal to the 139 samples after 7 presumed replicated samples


##########################################################
#
# Kunderfranco et al. Italian data
#
#########################################################

library(GEOquery)
gset <- getGEO("GSE14206", GSEMatrix =TRUE, getGPL=TRUE)
labels = Biobase::fData(gset[[1]])


# clinical
uncurated <- Biobase::pData(gset[[1]]) 

curated <- initial_curated_df(
  df_rownames = rownames(uncurated),
  template_name="data-raw/template_prad.csv")

curated <- curated %>% 
  dplyr::mutate(study_name = "Kunderfranco et al.") %>%
  dplyr::mutate(sample_name = uncurated$geo_accession) %>% 
  dplyr::mutate(patient_id = stringr::str_sub(uncurated$description.1, 22, 30)) %>%
  dplyr::mutate(alt_sample_name = uncurated$title) %>%
  dplyr::mutate(gleason_major = as.numeric(stringr::str_sub(uncurated$"gleason grade:ch1", 1, 1))) %>%
  dplyr::mutate(gleason_minor = as.numeric(stringr::str_sub(uncurated$"gleason grade:ch1", 3, 3))) %>%
  dplyr::mutate(gleason_grade = gleason_major + gleason_minor) %>%
  dplyr::mutate(grade_group = dplyr::case_when(
    gleason_grade  %in%  4:6 ~ "<=6",
    gleason_major == 3 & gleason_minor == 4 ~ "3+4",
    gleason_major == 4 & gleason_minor == 3 ~ "4+3",
    gleason_grade %in% 8:10 ~ ">=8",
    TRUE ~ "NA"
  )) %>%
  dplyr::mutate(age_at_initial_diagnosis = uncurated$"age:ch1") %>%
  dplyr::mutate(sample_type = dplyr::case_when(
                                               uncurated$source_name_ch1 == "prostate cancer" ~ 'primary',
                                               uncurated$source_name_ch1 == "normal prostate" ~ 'BPH')) %>%
  dplyr::mutate(frozen_ffpe = 'FFPE') %>%
  dplyr::mutate(source_of_gleason = dplyr::case_when(
                                                     uncurated$source_name_ch1 == "prostate cancer" ~ 'prostatectomy',
                                                     uncurated$source_name_ch1 == "normal prostate" ~ 'biopsy')) %>%
  dplyr::mutate(microdissected = 0) %>%
  dplyr::mutate(tissue_source = dplyr::case_when(
                                               uncurated$source_name_ch1 == "prostate cancer" ~ 'prostatectomy',
                                               uncurated$source_name_ch1 == "normal prostate" ~ 'biopsy'))
 

 

clinical_kunderfranco <- curated

save(clinical_kunderfranco, file = "./data-raw/clinical_kunderfranco.RData")

############################################################################


