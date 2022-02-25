library(magrittr)
library(dplyr)

initial_curated_df <- function(
  df_rownames,
  template_name
){
  # import template - drafted by Jim & Svitlana 
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

curated <- initial_curated_df(
  df_rownames = rownames(uncurated),
  template_name="data-raw/template_prad.csv")

## Note: 
#
# Clinical data obtained from the N=501 in https://www.cbioportal.org/study/clinicalData?id=prad_tcga,
# then supplemented by additional information available for N=334 in https://www.cbioportal.org/study/clinicalData?id=prad_tcga_pub
# then subset to samples we ultimately accept.
#
curated <- curated %>% 
  # Not just TCGA provisional with TCGA firehose brought in for GEX
  dplyr::mutate(study_name = "TCGA") %>%
  dplyr::mutate(sample_name = row.names(uncurated)) %>% 
  dplyr::mutate(frozen_ffpe = uncurated$IS_FFPE) %>%
  dplyr::mutate(frozen_ffpe = dplyr::case_when(
    frozen_ffpe %in% c("NO", "[Not Available]") ~ "NA",
    frozen_ffpe == "YES" ~ "ffpe",
    TRUE ~ frozen_ffpe
  )) %>%
  # Somatic differences had paired samples available, GEX was unpaired
  dplyr::mutate(sample_paired = 1) %>%
  dplyr::mutate(sample_type = dplyr::case_when(
  	uncurated$SAMPLE_TYPE == "Metastasis" ~ "metastatic",
  	uncurated$SAMPLE_TYPE == "Primary" ~ "primary"
  )) %>% 
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
    overall_survival_status == "1:DECEASED" ~ 1,
    overall_survival_status == "0:LIVING" ~ 0
  )) %>% 
  dplyr::mutate(days_to_overall_survival = as.numeric(uncurated$OS_MONTHS) * 30.5) %>%
  dplyr::mutate(disease_specific_recurrence_status = uncurated$DFS_STATUS) %>% 
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
  # single instance '[Unknown]' will throw a warning for below
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
  	dplyr::mutate(therapy_radiation_initial = uncurated$RADIATION_TREATMENT_ADJUVANT) %>%
  	# Radiation treatment given at initial treatment
  	dplyr::mutate(therapy_radiation_initial = dplyr::case_when(
  		therapy_radiation_initial == "YES" ~ 1,
  		therapy_radiation_initial == "NO" ~ 0,
  		therapy_radiation_initial == "" ~ NA_real_
  	)) %>%
  	dplyr::mutate(other_treatment = uncurated$TARGETED_MOLECULAR_THERAPY) %>%
  	# Add additional treatments
	dplyr::mutate(other_treatment = dplyr::case_when(
		other_treatment == "YES" ~ 1,
		other_treatment == "NO" ~ 0,
		other_treatment == "" ~ NA_real_
	)) %>%
	# Fraction of genome altered
	dplyr::mutate(genome_altered = uncurated$FRACTION_GENOME_ALTERED) %>%
	# Nx, N0 or N1 if findings in lymph nodes
	dplyr::mutate(N_stage = uncurated$PATH_N_STAGE) %>%
	dplyr::mutate(N_stage = dplyr::case_when(
		N_stage == "N0" ~ 0,
		N_stage == "N1" ~ 1,
		N_stage == "" ~ NA_real_
	))

# cgdsr::getGeneticProfiles(mycgds,'tcga_prad')

# Append additional information from Firehose Legacy; loop over all 'omics profiles
uncurated2 <- do.call("rbind", lapply(cgdsr::getGeneticProfiles(mycgds,'prad_tcga_pub')[,1], FUN=function(x){
	cgdsr::getClinicalData(mycgds, caseList=x)
}))
# Differences in reported fields (already based on GEX)
# > length(uncurated)
# [1] 80
# > length(uncurated2)
# [1] 91
#
curated2 <- initial_curated_df(
  df_rownames = rownames(uncurated2),
  template_name="data-raw/template_prad.csv")

# Match rownames between subsets
uncurated2 <- uncurated2[match(rownames(uncurated), rownames(uncurated2)),]

curated2 <- uncurated2 %>% 
	# ERG fusions were determined using gene expression
	dplyr::mutate(ERG_fusion_GEX = dplyr::case_when(
		ERG_STATUS == "fusion" ~ 1,
		ERG_STATUS == "none" ~ 0,
		is.na(ERG_STATUS) ~ NA_real_
	)) %>%
	# Purity estimates from TCGA Firehose legacy
	dplyr::mutate(tumor_purity_absolute = uncurated2$ABSOLUTE_EXTRACT_PURITY) %>%
	dplyr::mutate(tumor_purity_demixt = uncurated2$DEMIX_PURITY) %>%
	dplyr::mutate(tumor_purity_pathology = uncurated2$TUMOR_CELLULARITY_PATHOLOGY_REVIEW) %>%
	# Precomputed AR-scores
	dplyr::mutate(AR_activity = uncurated2$AR_SCORE) %>%
	# Ethnicity
	dplyr::mutate(race = dplyr::case_when(
		RACE == "ASIAN" ~ "asian",
		RACE == "BLACK OR AFRICAN AMERICAN" ~ "african_american",
		RACE == "WHITE" ~ "caucasian",
		is.na(RACE) ~ NA_character_
	))

# Copy additional curated field from Firehose legacy
curated[,"ERG_fusion_GEX"] <- curated2[,"ERG_fusion_GEX"]
curated[,"tumor_purity_absolute"] <- curated2[,"tumor_purity_absolute"]
curated[,"tumor_purity_demixt"] <- curated2[,"tumor_purity_demixt"]
curated[,"tumor_purity_pathology"] <- curated2[,"tumor_purity_pathology"]
curated[,"AR_activity"] <- curated2[,"AR_activity"]
curated[,"race"] <- curated2[,"race"]

# Subset to analytical subset of 333 samples prior to removing low quality samples
curated <- curated[which(curated$sample_name %in% intersect(rownames(uncurated), rownames(uncurated2))),]
rownames(curated) <- curated$sample_name

# Append information for the normal samples from GDC (XenaBrowser) as well as other potentially interesting fields
uncurated3 <- curatedPCaData:::generate_xenabrowser(id="TCGA-PRAD", type="clinical")
curated3 <- initial_curated_df(
  df_rownames = rownames(uncurated3),
  template_name="data-raw/template_prad.csv")

# GDC portion of clinical metadata (especially normal samples)
curated3 <- curated3 %>%
	dplyr::mutate(study_name = "TCGA") %>%
	dplyr::mutate(patient_id = gsub("-", ".", uncurated3$X_PATIENT)) %>%
	dplyr::mutate(sample_name = gsub(".01A", ".01", gsub("-", ".", sample_name))) %>%
	dplyr::mutate(alt_sample_name = uncurated3$file_uuid) %>%
	dplyr::mutate(age_at_initial_diagnosis = uncurated3$age_at_initial_pathologic_diagnosis) %>%
	dplyr::mutate(gleason_grade = uncurated3$gleason_score) %>%
	dplyr::mutate(disease_specific_recurrence_status = dplyr::case_when(
		uncurated3$biochemical_recurrence == "YES" ~ 1,
		uncurated3$biochemical_recurrence == "NO" ~ 0,
		is.na(uncurated3$biochemical_recurrence) | uncurated3$biochemical_recurrence == "" ~ NA_real_
	)) %>%
	dplyr::mutate(days_to_disease_specific_recurrence = uncurated3$days_to_first_biochemical_recurrence) %>%
	dplyr::mutate(psa = uncurated3$psa_value) %>%
	dplyr::mutate(overall_survival_status = uncurated3$OS) %>%
	dplyr::mutate(days_to_overall_survival = uncurated3$OS.time) %>%
	dplyr::mutate(gleason_major = uncurated3$primary_pattern) %>%
	dplyr::mutate(gleason_minor = uncurated3$secondary_pattern) %>%
	dplyr::mutate(grade_group = dplyr::case_when(
		uncurated3$primary_pattern + uncurated3$secondary_pattern <= 6 ~ "<=6",
		uncurated3$primary_pattern == 3 & uncurated3$secondary_pattern == 4 ~ "3+4",
		uncurated3$primary_pattern == 4 & uncurated3$secondary_pattern == 3 ~ "4+3",
		uncurated3$primary_pattern + uncurated3$secondary_pattern >= 8 ~ ">=8"
	)) %>%
	dplyr::mutate(sample_type = dplyr::case_when(
		uncurated3$sample_type.samples == "Metastatic" ~ "metastatic",
		uncurated3$sample_type.samples == "Primary Tumor" ~ "primary",
		uncurated3$sample_type.samples == "Solid Tissue Normal" ~ "healthy",
		is.na(uncurated3$sample_type.samples) | uncurated3$sample_type.samples == "" ~ NA_character_
	)) %>%
	dplyr::mutate(year_diagnosis = uncurated3$year_of_diagnosis.diagnoses) %>%
	dplyr::mutate(T_pathological = as.numeric(substr(uncurated3$pathologic_T, 2, 2))) %>%
	dplyr::mutate(T_substage_pathological = substr(uncurated3$pathologic_T, 3, 3)) %>%
	dplyr::mutate(zone_of_origin = dplyr::case_when(
		uncurated3$zone_of_origin == "Central Zone" ~ "central",
		uncurated3$zone_of_origin == "Overlapping / Multiple Zones" ~ "mixed",
		uncurated3$zone_of_origin == "Peripheral Zone" ~ "peripheral",
		uncurated3$zone_of_origin == "Transition Zone" ~ "transition",
		is.na(uncurated3$zone_of_origin) | uncurated3$zone_of_origin == "" ~ NA_character_
	)) %>%
	dplyr::mutate(frozen_ffpe = ifelse(uncurated3$is_ffpe.samples == "True", "FFPE", "frozen")) %>%
	dplyr::mutate(N_stage = as.numeric(substr(uncurated3$pathologic_N,2,2))) %>%
	dplyr::mutate(race = dplyr::case_when(
		uncurated3$race.demographic == "asian" ~ "asian",
		uncurated3$race.demographic == "american indian or alaska native" ~ "other",
		uncurated3$race.demographic == "black or african american" ~ "african_american",		
		uncurated3$race.demographic == "white" ~ "caucasian",
		is.na(uncurated3$race.demographic) | uncurated3$race.demographic == "not reported" ~ NA_character_
	)) %>%
	dplyr::mutate(therapy_radiation_initial = dplyr::case_when(
		uncurated3$radiation_therapy == "NO" ~ 0,
		uncurated3$radiation_therapy == "YES" ~ 0,
		is.na(uncurated3$radiation_therapy) | uncurated3$radiation_therapy == "" ~ NA_real_
	)) %>%
	dplyr::mutate(therapy_hormonal_initial = dplyr::case_when(
		uncurated3$additional_pharmaceutical_therapy == "YES" ~ 1,
		uncurated3$additional_pharmaceutical_therapy == "NO" ~ 0,
		is.na(uncurated3$additional_pharmaceutical_therapy) | uncurated3$additional_pharmaceutical_therapy == "" ~ NA_real_
	)) %>%
	dplyr::mutate(batch = uncurated3$batch_number) %>%
	# Append custom features of interest
	dplyr::mutate(other_feature = 
		paste(
			paste("primary_therapy_outcome_success=", uncurated3$primary_therapy_outcome_success, sep=""), 
			paste("additional_radiation_therapy=", uncurated3$additional_radiation_therapy, sep=""), 
		sep=";")
	)
# Note regarding BCR: sometimes patient has "biochemical_recurrence" == "NO" but has das to first/second_biochemical_recurrence (?)
# uncurated3[,c("biochemical_recurrence", "days_to_first_biochemical_recurrence", "days_to_second_biochemical_recurrence")]

# Append normal samples with GDC metadata to curated data from cBio's metadata
curated <- rbind(curated, curated3[which(curated3$sample_type == "healthy"),])

# Utilize additional information of low quality samples from pathology review and RNA degradation
# Samples excluded in pathology review
PathReview <- c(
	"TCGA-G9-6332",
	"TCGA-HC-7738",
	"TCGA-HC-7745",
	"TCGA-HC-7819",
	"TCGA-HC-8259",
	"TCGA-EJ-A46B",
	"TCGA-EJ-A46E",
	"TCGA-EJ-A46H",
	"TCGA-H9-A6BX",
	"TCGA-HC-7817",
	"TCGA-KK-A8IM",
	"TCGA-V1-A8MK",
	"TCGA-EJ-A8FO",
	"TCGA-EJ-A8FP",
	"TCGA-J4-A83J",
	"TCGA-KK-A8I7"
)
# Samples indicated to suffer from RNA degradation
RNAdegrade <- c(
	"TCGA-J4-A67O-01A-11R-A30B-07",
	"TCGA-QU-A6IM-01A-11R-A31N-07",
	"TCGA-QU-A6IL-01A-11R-A31N-07",
	"TCGA-KK-A7AW-01A-11R-A32O-07",
	"TCGA-G9-6347-01A-11R-A31N-07",
	"TCGA-EJ-A46F-01A-31R-A250-07",
	"TCGA-QU-A6IO-01A-11R-A31N-07",
	"TCGA-J4-A67Q-01A-21R-A30B-07",
	"TCGA-QU-A6IN-01A-11R-A31N-07",
	"TCGA-G9-6379-01A-11R-A31N-07",
	"TCGA-KK-A59V-01A-11R-A29R-07",
	"TCGA-EJ-7325-01B-11R-A32O-07",
	"TCGA-HC-A6HX-01A-11R-A31N-07",
	"TCGA-G9-6362-01A-11R-1789-07",
	"TCGA-J4-A67N-01A-11R-A30B-07",
	"TCGA-FC-A6HD-01A-11R-A31N-07",
	"TCGA-KK-A7AY-01A-11R-A33R-07",
	"TCGA-G9-6498-01A-12R-A311-07",
	"TCGA-KK-A7B0-01A-11R-A32O-07",
	"TCGA-V1-A9OT-01A-11R-A41O-07",
	"TCGA-HC-A6AQ-01A-11R-A30B-07",
	"TCGA-EJ-A65D-01A-11R-A30B-07",
	"TCGA-J4-A67R-01A-21R-A30B-07",
	"TCGA-J4-A67M-01A-11R-A30B-07",
	"TCGA-EJ-A65M-01A-11R-A29R-07",
	"TCGA-MG-AAMC-01A-11R-A41O-07",
	"TCGA-M7-A723-01A-12R-A32O-07",
	"TCGA-EJ-A65B-01A-12R-A30B-07",
	"TCGA-KK-A6E3-01A-21R-A30B-07",
	"TCGA-H9-A6BY-01A-11R-A30B-07",
	"TCGA-EJ-AB20-01A-12R-A41O-07",
	"TCGA-J4-A67S-01A-11R-A30B-07",
	"TCGA-HC-A6AP-01A-11R-A30B-07",
	"TCGA-X4-A8KQ-01A-12R-A36G-07",
	"TCGA-HC-7752-01A-11R-2118-07",
	"TCGA-G9-6496-01A-11R-1789-07",
	"TCGA-HC-A6AN-01A-11R-A30B-07",
	"TCGA-FC-A66V-01A-21R-A30B-07",
	"TCGA-HC-A6AS-01A-11R-A30B-07",
	"TCGA-HC-A6AO-01A-11R-A30B-07",
	"TCGA-HI-7169-01A-11R-2118-07",
	"TCGA-J4-A67L-01A-11R-A30B-07",
	"TCGA-HC-A6HY-01A-11R-A31N-07",
	"TCGA-G9-6354-01A-11R-A311-07",
	"TCGA-J4-A67K-01A-21R-A30B-07",
	"TCGA-G9-6369-01A-21R-1965-07",
	"TCGA-KK-A6E4-01A-11R-A30B-07",
	"TCGA-KK-A6E7-01A-11R-A31N-07",
	"TCGA-G9-6339-01A-12R-A311-07",
	"TCGA-G9-6373-01A-11R-1789-07",
	"TCGA-KK-A7AQ-01A-11R-A33R-07",
	"TCGA-HC-A6AL-01A-11R-A30B-07",
	"TCGA-ZG-A9L9-01A-11R-A41O-07",
	"TCGA-KK-A6E8-01A-11R-A31N-07",
	"TCGA-2A-A8VX-01A-11R-A37L-07",
	"TCGA-G9-6343-01A-21R-1965-07",
	"TCGA-KC-A7F5-01A-11R-A33R-07",
	"TCGA-J9-A52E-01A-11R-A26U-07",
	"TCGA-J9-A52C-01A-11R-A26U-07",
	"TCGA-XJ-A9DX-01A-11R-A37L-07",
	"TCGA-G9-6338-01A-12R-1965-07",
	"TCGA-M7-A724-01A-12R-A32O-07",
	"TCGA-KK-A7B2-01A-12R-A32O-07",
	"TCGA-ZG-A9LB-01A-11R-A41O-07",
	"TCGA-V1-A8MM-01A-11R-A37L-07",
	"TCGA-M7-A71Z-01A-12R-A32O-07",
	"TCGA-XK-AAK1-01A-11R-A41O-07",
	"TCGA-VN-A88M-01A-11R-A352-07",
	"TCGA-KK-A7AZ-01A-12R-A32O-07",
	"TCGA-ZG-A9LM-01A-11R-A41O-07",
	"TCGA-XJ-A83F-01A-11R-A352-07",
	"TCGA-XJ-A9DQ-01A-11R-A37L-07",
	"TCGA-ZG-A9LU-01A-11R-A41O-07",
	"TCGA-G9-7525-01A-31R-2263-07",
	"TCGA-G9-6496-11A-01R-1789-07", # Normal sample
	"TCGA-G9-6362-11A-01R-1789-07", # Normal sample
	"TCGA-G9-6348-11A-01R-1789-07" # Normal sample
)
PathReview <- paste(gsub("-", ".", PathReview), ".01", sep="")
RNAdegrade <- gsub("-", ".", RNAdegrade)
RNAdegrade <- gsub(".01A", ".01", substr(RNAdegrade, start=0, stop=16))
# Exclusion of low quality samples (3 normals, after the larger data was subset to provisional tumor samples)
curated <- curated[which(!curated$sample_name %in% c(PathReview, RNAdegrade)),]

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
  dplyr::mutate(sample_paired = 0) %>%
  dplyr::mutate(sample_type = "primary") %>%
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

gse_all <- GEOquery::getGEO("GSE21032", GSEMatrix = TRUE)


gse <- rbind(
	# Transcript
	Biobase::pData(gse_all[[1]])[,c("sample id:ch1", "tissue:ch1", "tumor type:ch1", "biopsy_gleason_grade:ch1", "clint_stage:ch1", "disease status:ch1", "pathological_stage:ch1")],
	# aCGH
	Biobase::pData(gse_all[[2]])[,c("sample id:ch1", "tissue:ch1", "tumor type:ch1", "biopsy_gleason_grade:ch1", "clint_stage:ch1", "disease status:ch1", "pathological_stage:ch1")],
	# Exon
	Biobase::pData(gse_all[[3]])[,c("sample id:ch1", "tissue:ch1", "tumor type:ch1", "biopsy_gleason_grade:ch1", "clint_stage:ch1", "disease status:ch1", "pathological_stage:ch1")],
	# miRNA
	Biobase::pData(gse_all[[4]])[,c("sample id:ch1", "tissue:ch1", "tumor type:ch1", "biopsy_gleason_grade:ch1", "clint_stage:ch1", "disease status:ch1", "pathological_stage:ch1")]
)



#gse_gex <- GEOquery::getGEO("GSE21034", GSEMatrix = TRUE)
#uncurated_gex_exon <- Biobase::pData(gse_gex[[2]])
#uncurated_gex_transcript <- Biobase::pData(gse_gex[[1]])

#gse_cna <- GEOquery::getGEO("GSE21035", GSEMatrix = TRUE)
#uncurated_cna <- Biobase::pData(gse_cna[[1]])

#gse_mrna <- GEOquery::getGEO("GSE21036", GSEMatrix = TRUE)
#uncurated_mrna <- Biobase::pData(gse_mrna[[1]])

# mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
# uncurated_cbio <- cgdsr::getClinicalData(mycgds, caseList = "prad_mskcc_all")
mae <-cBioPortalData::cBioDataPack("prad_mskcc",ask = FALSE)
uncurated=colData(mae)
uncurated_cbio=as.data.frame(uncurated)
#rownames(uncurated)<-gsub("-",".",rownames(uncurated))

# Outdated names no longer in use (GSM######)
if(FALSE){
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
					      paste0("mrna: ", mrna),
					      # TDL: Mutation sample names come from cBio rather than GEO, but we don't know which ones will have mutation data
					      paste0("mut: ", uncurated$"sample id:ch1"), sep = "|"))

	rownames(uncurated) <- uncurated$title
}

curated <- initial_curated_df(
  df_rownames = gse$"sample id:ch1",
  template_name="data-raw/template_prad.csv")

curated <- curated %>% 
  dplyr::mutate(study_name = "Taylor, et al.") %>% 
  #dplyr::mutate(alt_sample_name = uncurated$geo_accession) %>% 
  #dplyr::mutate(sample_name = row.names(uncurated)) %>%
  dplyr::mutate(patient_id = gse$"sample id:ch1") %>%
  dplyr::mutate(sample_paired = 0) %>%
  dplyr::mutate(sample_type = tolower(gse$"tumor type:ch1")) %>%
  dplyr::mutate(sample_type = dplyr::case_when(
    sample_type == "primary tumor" ~ "primary",
    sample_type == "cell line" ~ "cell.line",
    sample_type == "xenograft" ~ "xenograft",
    sample_type == "metastatsis" ~ "metastasis",
    is.na(sample_type) ~ "normal",
    TRUE ~ sample_type
  )) %>%
  dplyr::mutate(gleason_grade = as.numeric(gse$"biopsy_gleason_grade:ch1")) %>%
  dplyr::mutate(T_clinical = stringr::str_sub(gse$"clint_stage:ch1",2,2)) %>%
  dplyr::mutate(T_clinical = as.numeric(T_clinical)) %>%
  dplyr::mutate(T_substage_clinical = tolower(stringr::str_sub(gse$"clint_stage:ch1",3,3))) %>%
  dplyr::mutate(T_substage_clinical = dplyr::case_when(
    T_substage_clinical == "" ~ NA_character_,
    TRUE ~ T_substage_clinical
  )) %>%
  dplyr::mutate(T_pathological = as.numeric(stringr::str_sub(gse$"pathological_stage:ch1",2,2))) %>%
  dplyr::mutate(T_substage_pathological = tolower(stringr::str_sub(gse$"pathological_stage:ch1",3,3))) %>%
  dplyr::mutate(T_substage_pathological = dplyr::case_when(
    T_substage_pathological == "" ~ NA_character_,
    TRUE ~ T_substage_pathological
  ))

# Order additional clinical information from cBio according to the order in current df
uncurated_cbio <- uncurated_cbio[match(curated$patient_id, row.names(uncurated_cbio)),]

# Additional clinical parameters as recored in cBioPortal
curated <- curated %>%
	# Gleason grades reported in GEO seem to differ from cBioPortal even for same IDs; using the ones provided by cBio:
	dplyr::mutate(gleason_grade = uncurated_cbio$GLEASON_SCORE) %>% 
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
		disease_specific_recurrence_status == "1:Recurred" ~ 1,
		disease_specific_recurrence_status == "0:DiseaseFree" ~ 0,
		TRUE ~ NA_real_
	)) %>% 
	dplyr::mutate(days_to_disease_specific_recurrence = round(uncurated_cbio$DFS_MONTHS*30.5,0)) %>%
	dplyr::mutate(genome_altered = uncurated_cbio$FRACTION_GENOME_ALTERED)

## Old filters when normals weren't requested
#curated <- curated %>%
#  dplyr::filter(sample_type %in% c("metastasis", "primary"))

# Leave out cell.line and xenograft samples
curated <- curated[which(curated$sample_type %in% c("primary", "metastasis", "normal")),]
curated <- curated[grep("PCA|PAN", curated$sample_name),]
curated <- curated[order(curated$sample_name),]
# Only include unique entries
curated <- curated[which(!duplicated(curated$patient_id)),]


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
  dplyr::mutate(sample_paired = 0) %>%
  dplyr::mutate(sample_type = "primary") %>%
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
  template_name="data-raw/template_prad.csv")

curated <- curated %>% 
  dplyr::mutate(study_name = "Friedrich et al.") %>%
  dplyr::mutate(sample_name = row.names(uncurated)) %>% 
  dplyr::mutate(sample_paired = 0) %>%
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

extra_cl = read.table('data-raw/GSE134051_Pheno_Data_ControlType.txt', h = T) # extra clinical from Friedrich et al. 
extra_cl[extra_cl$PathoState == 'None', 'PathoState'] = NA
extra_cl$PathoState = as.numeric(extra_cl$PathoState)
curated[extra_cl$Control_status == 'BLADDER_CANCER', 'tissue_source'] = 'cystoprostatectomy'
curated[extra_cl$Control_status == 'OTHER', 'tissue_source'] = 'TURP'
curated[, 'T_pathological'] = extra_cl$PathoState

clinical_friedrich <- curated

save(clinical_friedrich, file = "data-raw/clinical_friedrich.RData")

##################################################################
#
# Wallace et al. samples: This code downloads both gene expression
# BOTH clinical and GPL, but just the clinical data is used
#
##################################################################

gset <- getGEO("GSE6956", GSEMatrix =TRUE, getGPL=FALSE)

uncurated <- Biobase::pData(gset[[1]])

unmatched_healty_tissue = c('GSM160418', 'GSM160419', 'GSM160420', 'GSM160421', 'GSM160422', 'GSM160430')

uncurated = uncurated[!is.element(uncurated$geo_accession, unmatched_healty_tissue), ]

curated <- initial_curated_df(
  df_rownames = rownames(uncurated),
  template_name="data-raw/template_prad.csv")

curated <- curated %>%
  dplyr::mutate(study_name = "Wallace et al.") %>%
  dplyr::mutate(sample_name = row.names(uncurated)) %>%
  dplyr::mutate(sample_paired = 0) %>%
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

save(clinical_wallace, file = "data-raw/clinical_wallace.RData")


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
  template_name="data-raw/template_prad.csv")

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

save(clinical_chandran, file = "data-raw/clinical_chandran.RData")

######################################################################
#cBioportal Barbieri Broad/Cornell Data
#####################################################################
mae <-cBioPortalData::cBioDataPack("prad_broad",ask = FALSE)
uncurated=colData(mae)
uncurated=as.data.frame(uncurated)
rownames(uncurated)<-gsub("-",".",rownames(uncurated))
# mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
# uncurated <- cgdsr::getClinicalData(mycgds, caseList="prad_broad_sequenced")

# create the curated object
curated <- initial_curated_df(
  df_rownames = rownames(uncurated),
  template_name="data-raw/template_prad.csv")

curated <- curated %>% 
  dplyr::mutate(study_name = "Barbieri et al.") %>%
  dplyr::mutate(sample_name = row.names(uncurated)) %>%
  dplyr::mutate(patient_id = row.names(uncurated)) %>%
  # from the publication: Primary treatment naive radical prostatectomy specimen from American and Australian patients
  dplyr::mutate(sample_paired = 1) %>%
  dplyr::mutate(sample_type = "primary") %>%
  dplyr::mutate(race = "caucasian") %>%
  dplyr::mutate(tissue_source = "prostatectomy") %>%
  dplyr::mutate(gleason_major = as.integer(stringr::str_sub(uncurated$GLEASON_SCORE,1,1))) %>%
  dplyr::mutate(gleason_minor = as.integer(stringr::str_sub(uncurated$GLEASON_SCORE,3,3))) %>%
  dplyr::mutate(gleason_grade = gleason_major+gleason_minor) %>%
  #dplyr::mutate(gleason_grade = gsub(";.*","",uncurated$GLEASON_SCORE)) %>%
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
  # from the publication: TMPRSS2-ERG fusion status assessed by fluorescence in situ hybridization (FISH)
  #dplyr::mutate(TMPRSS2_ERG_FUSION_STATUS=uncurated$TMPRSS2_ERG_FUSION_STATUS) %>%
  dplyr::mutate(ERG_fusion_IHC = dplyr::case_when(
  	uncurated$TMPRSS2_ERG_FUSION_STATUS == "Negative" ~ 0, # Negative -> 0
  	TRUE ~ 1 # Else -> 1
  )) %>%
  dplyr::mutate(GLEASON_SCORE_PERCENT_4_AND_5=uncurated$GLEASON_SCORE_PERCENT_4_AND_5)

clinical_barbieri <- curated

save(clinical_barbieri, file = "data-raw/clinical_barbieri.RData")

#####################################################################################
#cBioportal Ren eururol 2017 Data
#####################################################################################

# mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
# uncurated <- cgdsr::getClinicalData(mycgds, caseList="prad_eururol_2017_sequenced")
mae <-cBioPortalData::cBioDataPack("prad_eururol_2017",ask = FALSE)
uncurated=colData(mae)
uncurated=as.data.frame(uncurated)
#rownames(uncurated)<-gsub("-",".",rownames(uncurated))

# create the curated object
curated <- initial_curated_df(
  df_rownames = rownames(uncurated),
  template_name="data-raw/template_prad.csv")

curated <- curated %>% 
  dplyr::mutate(study_name = "Ren et al.") %>%
  dplyr::mutate(sample_name = row.names(uncurated)) %>%
  # From the publication: "The study sequenced whole-genome and transcriptome of tumor-benign paired tissues from 65 treatment-naive Chinese PCa patients"
  dplyr::mutate(sample_paired = 1) %>%
  dplyr::mutate(sample_type = "primary") %>%
  dplyr::mutate(race = "asian") %>%
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
  # Tumor purity was estimated using ASCAT (Allele-Specific Copy number Analysis of Tumours)
  dplyr::mutate(ASCAT_Tumor_Purity = uncurated$TUMOR_PURITY) %>%
  dplyr::mutate(BLADDER_NECK_INVASION = uncurated$BLADDER_NECK_INVASION) %>%
  dplyr::mutate(seminal_vesicle_invasion = dplyr::case_when(
  	uncurated$SEMINAL_VESICLE_INVASION == "No" ~ 0,
  	uncurated$SEMINAL_VESICLE_INVASION == "Yes" ~ 1,
  	TRUE ~ NA_real_  	
  )) %>%
  dplyr::mutate(genome_altered=uncurated$FRACTION_GENOME_ALTERED) %>%
  dplyr::mutate(N_stage = dplyr::case_when(
  	uncurated$LYMPH_NODE_METASTASIS == "No" ~ 0,
  	uncurated$LYMPH_NODE_METASTASIS == "Yes" ~ 1,
  	uncurated$LYMPH_NODE_METASTASIS == "" ~ NA_real_
  )) %>%
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
  dplyr::mutate(tissue_source = gsub("prostate cancer biopsy", "biopsy", stringr::str_remove(uncurated$characteristics_ch1.10,"tissue:"))) %>%
  dplyr::mutate(age_at_initial_diagnosis = floor(as.numeric(uncurated$`age:ch1`))) %>%
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
  dplyr::mutate(T_substage_clinical = stringr::str_extract(uncurated$`tumor stage:ch1`, "[a-c]+")) %>%
  dplyr::mutate(sample_type = "primary")
  
clinical_kim <- curated

save(clinical_kim, file = "data-raw/clinical_kim.RData")
  

#######################################################################
#
# Abida et al. 2019
#
#######################################################################
mae <-cBioPortalData::cBioDataPack("prad_su2c_2019",ask = FALSE)
uncurated=colData(mae)
uncurated=as.data.frame(uncurated)
uncurated$SAMPLE_ID<-gsub("-",".",uncurated$SAMPLE_ID)

# mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
# uncurated <- cgdsr::getClinicalData(mycgds, caseList="prad_su2c_2019_cnaseq")
# gp = cgdsr::getGeneticProfiles(mycgds,"prad_su2c_2019")

curated <- initial_curated_df(
  df_rownames = uncurated$SAMPLE_ID,
  template_name="data-raw/template_prad.csv")

curated <- curated %>% 
  dplyr::mutate(study_name = "Abida et al.") %>%
  dplyr::mutate(alt_sample_name = uncurated$OTHER_SAMPLE_ID) %>%
  dplyr::mutate(patient_id = uncurated$PATIENT_ID) %>%
  # Samples were paired for analyses
  dplyr::mutate(sample_paired = 1) %>%
  dplyr::mutate(sample_name = uncurated$SAMPLE_ID) %>%
  ## TDL: See template_prad.csv, tissue_source does not mean tissue extraction site but rather method of extraction - follow template_prad.csv instructions for all fields
  #dplyr::mutate(tissue_source = uncurated$TISSUE_SITE) %>%
  dplyr::mutate(age_at_initial_diagnosis = floor(uncurated$AGE_AT_DIAGNOSIS)) %>%
  dplyr::mutate(age_at_procurement = floor(uncurated$AGE_AT_PROCUREMENT)) %>%
  # All samples are metastatic
  dplyr::mutate(metastasis_occurrence_status = 1) %>%
  dplyr::mutate(sample_type = "metastatic") %>%
  dplyr::mutate(metastatic_site = dplyr::case_when(
  	uncurated$TISSUE_SITE == "Adrenal" ~ "adrenal_gland",
  	uncurated$TISSUE_SITE == "Bone" ~ "bone",
  	uncurated$TISSUE_SITE == "Brain" ~ "brain",
  	uncurated$TISSUE_SITE == "Liver" ~ "liver",
  	uncurated$TISSUE_SITE == "LN" ~ "lymph_node",
  	uncurated$TISSUE_SITE == "Lung" ~ "lung",
  	uncurated$TISSUE_SITE == "Other Soft tissue" ~ "soft_tissue",
  	uncurated$TISSUE_SITE == "Prostate" ~ "prostate",
  	uncurated$TISSUE_SITE == "Unknown" ~ "other")) %>%
  dplyr::mutate(psa = uncurated$PSA) %>%
  #dplyr::mutate(gleason_score=uncurated$GLEASON_SCORE) %>%
  dplyr::mutate(gleason_grade=as.numeric(uncurated$GLEASON_SCORE)) %>%
  dplyr::mutate(AR_score=uncurated$AR_SCORE) %>%
  dplyr::mutate(genome_altered = uncurated$FRACTION_GENOME_ALTERED) %>%
  dplyr::mutate(NEPC_score=uncurated$NEPC_SCORE) %>%
  dplyr::mutate(ABI_ENZA_exposure_status=uncurated$ABI_ENZA_EXPOSURE_STATUS) %>%
  dplyr::mutate(ETS_FUSION_DETAILS=uncurated$ETS_FUSION_DETAILS) %>%
  #dplyr::mutate(TMPRSS2_ERG_FUSION_STATUS = dplyr::case_when(
  ## Fusions were determined from gene expression
  dplyr::mutate(ERG_fusion_GEX = dplyr::case_when(
    uncurated$ETS_FUSION_DETAILS == "TMPRSS2-ERG" ~ 1,
    uncurated$ETS_FUSION_DETAILS != "TMPRSS2-ERG" ~ 0)) %>%
  dplyr::mutate(overall_survival_status = dplyr::case_when(
    uncurated$OS_STATUS == "1:DECEASED" ~ 1,
    uncurated$OS_STATUS == "0:LIVING" ~ 0)) %>%
  dplyr::mutate(days_to_overall_survival=as.numeric(uncurated$OS_MONTHS)*30.5) %>%
  ## Possibly TODO: AGE_AT_PROCUREMENT may be an indicator for time from non-metastatic to metastatic disease (difference in years)
  dplyr::mutate(other_treatment = uncurated$CHEMO_REGIMEN_CATEGORY)

clinical_abida <- curated
save(clinical_abida, file = "data-raw/clinical_abida.RData")


####
#
# Barwick et al.
#
####

gse <- GEOquery::getGEO("GSE18655", GSEMatrix = TRUE)

uncurated <- Biobase::pData(gse[[1]])

curated <- initial_curated_df(
	df_rownames = rownames(uncurated),
	template_name="data-raw/template_prad.csv"
)

curated <- curated %>% 
	dplyr::mutate(study_name = "Barwick et al.") %>%
	dplyr::mutate(patient_id = paste0("X", uncurated$'title')) %>%
	dplyr::mutate(sample_name = paste0("X", uncurated$'title')) %>%
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
	dplyr::mutate(T_pathological = as.numeric(uncurated$'grade:ch1')) %>% 
	dplyr::mutate(sample_type = "primary") %>%
	dplyr::mutate(sample_paired = 0)

# Replicates provided in the raw gene expression matrix:
#> grep("rep", colnames(gex), value=TRUE)
# [1] "X1_rep1"   "X1_rep2"   "X100_rep1" "X100_rep2" "X105_rep1" "X105_rep2" "X151_rep1" "X151_rep2" "X172_rep1" "X172_rep2" "X61_rep1"  "X61_rep2"  "X77_rep1"  "X77_rep2"  "X9_rep1"   "X9_rep2"
#> unique(unlist(lapply(grep("rep", colnames(gex), value=TRUE), FUN=function(x) { strsplit(x, "_")[[1]][1]})))
#[1] "X1"   "X100" "X105" "X151" "X172" "X61"  "X77"  "X9"

# Append correct "_rep1", "_rep2" suffixes to correct samples
curated <- rbind(curated, curated[which(curated$patient_id %in% c("X1", "X100", "X105", "X151", "X172", "X61", "X77", "X9")),])
curated[curated$sample_name %in% c("X1", "X100", "X105", "X151", "X172", "X61", "X77", "X9"),"sample_name"] <- paste0(curated$sample_name[curated$sample_name %in% c("X1", "X100", "X105", "X151", "X172", "X61", "X77", "X9")], rep(c("_rep1", "_rep2"), each=7))
rownames(curated) <- NULL
  
clinical_barwick <- curated

save(clinical_barwick, file = "./data-raw/clinical_barwick.RData")

##########################################################
#
# Kunderfranco et al. Italian data
#
#########################################################

library(GEOquery)
gset <- getGEO("GSE14206", GSEMatrix =TRUE, getGPL=TRUE)
#labels = Biobase::fData(gset[[1]])


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
  dplyr::mutate(sample_paired = 1) %>%
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

###########################################################################################
#
#
#  True et al. Proc Natl Acad Sci U S A 2006 Jul 18;103(29):10991-6.
#
#  this is an exceedingly old dataset, where the expression data is relative and not absolute
#  which makes the data not directly comparable with the other dataset
#  
#     
#
###########################################################################################

library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE5132", GSEMatrix =TRUE, getGPL=TRUE)

# for reasons unknown the clinical is split in a set of 31 samples and a set of 1 sample
uncurated1 <- Biobase::pData(gset[[1]]) 
uncurated2 <- Biobase::pData(gset[[2]]) 

curated1 <- initial_curated_df(
  df_rownames = rownames(uncurated1),
  template_name="data-raw/template_prad.csv")

curated2 <- initial_curated_df(
  df_rownames = rownames(uncurated2),
  template_name="data-raw/template_prad.csv")

#  additional issue is that, because the data seems to have been put in at random, most of it will have
# to be pulled 'manually' rather than algorithmically
curated1 <- curated1 %>% 
  dplyr::mutate(study_name = "True et al.") %>%
  dplyr::mutate(sample_name = uncurated1$geo_accession) %>% 
  dplyr::mutate(sample_paired = 1) %>%
  dplyr::mutate(sample_type = "primary") %>%
  dplyr::mutate(patient_id = unlist(lapply(stringr::str_split(uncurated1$title, '_'), function(x) x[2]))) %>%
  dplyr::mutate(gleason_major = as.numeric(substr(unlist(lapply(stringr::str_split(uncurated1$description, ' '), function(x) x[9])), 1, 1))) %>%
  dplyr::mutate(gleason_minor = as.numeric(substr(unlist(lapply(stringr::str_split(uncurated1$description, ' '), function(x) x[9])), 3, 3))) %>%
  dplyr::mutate(gleason_grade = gleason_major + gleason_minor) %>%
  dplyr::mutate(microdissected = 1) %>%
  dplyr::mutate(gleason_group = dplyr::case_when(
    gleason_grade  %in%  4:6 ~ "<=6",
    gleason_major == 3 & gleason_minor == 4 ~ "3+4",
    gleason_major == 4 & gleason_minor == 3 ~ "4+3",
    gleason_grade %in% 8:10 ~ ">=8",
    TRUE ~ "NA"
  ))
curated1[1,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[1,'Prostate Cancer patient 03-138A Gleason_Score:ch1'], ' '))[5], 5, 8))
curated1[1, 'tumor_margins_positive'] = 0
curated1[1, 'other_treatment'] = substr(unlist(stringr::str_split(uncurated1[1,'Prostate Cancer patient 03-138A Gleason_Score:ch1'], ' '))[8], 11, 20)
curated1[1, 'therapy_radiation_initial'] = 0
curated1[1, 'therapy_radiation_salvage'] = 0
curated1[1, 'therapy_hormonal_initial'] = 0
curated1[1, 'other_feature'] = paste(
                                     unlist(stringr::str_split(uncurated1[1,'Prostate Cancer patient 03-138A Gleason_Score:ch1'], ' '))[4],
                                     unlist(stringr::str_split(uncurated1[1,'Prostate Cancer patient 03-138A Gleason_Score:ch1'], ' '))[6],
                                     unlist(stringr::str_split(uncurated1[1,'Prostate Cancer patient 03-138A Gleason_Score:ch1'], ' '))[2],
                                     sep = '|')
curated[1, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[1,'Prostate Cancer patient 03-138A Gleason_Score:ch1'], ' '))[3], 17, 17))

curated1[2,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[2,'Prostate Cancer patient 03-135C Gleason_Score:ch1'], ' '))[5], 5, 8))
curated1[2, 'tumor_margins_positive'] = 0
curated1[2, 'therapy_radiation_initial'] = 0
curated1[2, 'therapy_radiation_salvage'] = 0
curated1[2, 'therapy_hormonal_initial'] = 0
curated1[2, 'other_feature'] = paste(
                                     unlist(stringr::str_split(uncurated1[2,'Prostate Cancer patient 03-135C Gleason_Score:ch1'], ' '))[4],
                                     unlist(stringr::str_split(uncurated1[2,'Prostate Cancer patient 03-135C Gleason_Score:ch1'], ' '))[6],
                                     unlist(stringr::str_split(uncurated1[2,'Prostate Cancer patient 03-135C Gleason_Score:ch1'], ' '))[2],
                                     sep = '|')
curated1[2, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[2,'Prostate Cancer patient 03-135C Gleason_Score:ch1'], ' '))[3], 17, 17))

curated1[3,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[3,'Prostate Cancer patient 03-158F Gleason_Score:ch2'], ' '))[5], 5, 8))
curated1[3, 'tumor_margins_positive'] = 0
curated1[3, 'therapy_radiation_initial'] = 0
curated1[3, 'therapy_radiation_salvage'] = 0
curated1[3, 'therapy_hormonal_initial'] = 0
curated1[3, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[3,'Prostate Cancer patient 03-158F Gleason_Score:ch2'], ' '))[4],
                                     unlist(stringr::str_split(uncurated1[3,'Prostate Cancer patient 03-158F Gleason_Score:ch2'], ' '))[6],
                                     unlist(stringr::str_split(uncurated1[3,'Prostate Cancer patient 03-158F Gleason_Score:ch2'], ' '))[2],
                                     sep = '|')
curated1[3, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[3,'Prostate Cancer patient 03-158F Gleason_Score:ch2'], ' '))[3], 17, 17))


curated1[4,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[4,' LCM_Gleason_Pattern:ch2'], ' '))[4], 5, 8))
curated1[4, 'tumor_margins_positive'] = 0
curated1[4, 'other_treatment'] = paste(substr(unlist(stringr::str_split(uncurated1[4,' LCM_Gleason_Pattern:ch2'], ' '))[7], 11, 20), unlist(stringr::str_split(uncurated1[4,' LCM_Gleason_Pattern:ch2'], ' '))[8], sep = '_')
curated1[4, 'therapy_radiation_initial'] = 0
curated1[4, 'therapy_radiation_salvage'] = 0
curated1[4, 'therapy_hormonal_initial'] = 0
curated1[4, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[4,' LCM_Gleason_Pattern:ch2'], ' '))[3],
                                     unlist(stringr::str_split(uncurated1[4,' LCM_Gleason_Pattern:ch2'], ' '))[5],
                                     unlist(stringr::str_split(uncurated1[4, 'characteristics_ch2.1'], ' '))[2],
                                     sep = '|')
curated1[4, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[4,' LCM_Gleason_Pattern:ch2'], ' '))[2], 17, 17))


curated1[5,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[5,' LCM_Gleason_Pattern:ch1'], ' '))[4], 5, 8))
curated1[5, 'tumor_margins_positive'] = 0
curated1[5, 'therapy_radiation_initial'] = 0
curated1[5, 'therapy_radiation_salvage'] = 0
curated1[5, 'therapy_hormonal_initial'] = 0
curated1[5, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[5,' LCM_Gleason_Pattern:ch1'], ' '))[3],
                                     unlist(stringr::str_split(uncurated1[5,' LCM_Gleason_Pattern:ch1'], ' '))[5],
                                     unlist(stringr::str_split(uncurated1[5,'characteristics_ch1.1'], ' '))[2],
                                     sep = '|')
curated1[6, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[5,' LCM_Gleason_Pattern:ch1'], ' '))[2], 17,17))


curated1[6,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[6,'Prostate Cancer patient 03-060A Gleason_Score:ch1'], ' '))[5], 5, 8))
curated1[6, 'tumor_margins_positive'] = 1
curated1[6, 'therapy_radiation_initial'] = 0
curated1[6, 'therapy_radiation_salvage'] = 0
curated1[6, 'therapy_hormonal_initial'] = 0
curated1[6, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[6,'Prostate Cancer patient 03-060A Gleason_Score:ch1'], ' '))[4], 
                                     unlist(stringr::str_split(uncurated1[6,'Prostate Cancer patient 03-060A Gleason_Score:ch1'], ' '))[6],
                                     unlist(stringr::str_split(uncurated1[6,'Prostate Cancer patient 03-060A Gleason_Score:ch1'], ' '))[2],
                                     sep = '|')
curated1[6, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[6,'Prostate Cancer patient 03-060A Gleason_Score:ch1'], ' '))[3], 17, 17))

curated1[7,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[7,' LCM_Gleason_Pattern:ch2'], ' '))[4], 5, 8))
curated1[7, 'tumor_margins_positive'] = 1
curated1[7, 'therapy_radiation_initial'] = 0
curated1[7, 'therapy_radiation_salvage'] = 0
curated1[7, 'therapy_hormonal_initial'] = 0
curated1[7, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[7,' LCM_Gleason_Pattern:ch2'], ' '))[3], 
                                     unlist(stringr::str_split(uncurated1[7,' LCM_Gleason_Pattern:ch2'], ' '))[5], 
                                     unlist(stringr::str_split(uncurated1[7,'characteristics_ch2.1'], ' '))[2], 
                                     sep = '|')
curated1[7, 'grade_group'] =  as.numeric(substr(unlist(stringr::str_split(uncurated1[7,'characteristics_ch2.1'], ' '))[3], 17, 17))

curated1[8,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[8,'Prostate Cancer patient 03-066C Gleason_Score:ch2'], ' '))[8], 5, 8))
curated1[8, 'tumor_margins_positive'] = 1
curated1[8, 'therapy_radiation_initial'] = 0
curated1[8, 'therapy_radiation_salvage'] = 0
curated1[8, 'therapy_hormonal_initial'] = 0
curated1[8, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[8,'Prostate Cancer patient 03-066C Gleason_Score:ch2'], ' '))[7],
                                     unlist(stringr::str_split(uncurated1[8,'Prostate Cancer patient 03-066C Gleason_Score:ch2'], ' '))[9],
                                     unlist(stringr::str_split(uncurated1[8,'Prostate Cancer patient 03-066C Gleason_Score:ch2'], ' '))[5],
                                     sep = '|')
curated1[8, 'grade_group'] = as.numeric(substr( unlist(stringr::str_split(uncurated1[8,'Prostate Cancer patient 03-066C Gleason_Score:ch2'], ' '))[6], 17, 17))

curated1[9,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[9,'Prostate Cancer patient 03-159 Gleason_Score:ch1'], ' '))[5], 5, 8))
curated1[9, 'tumor_margins_positive'] = 1
curated1[9, 'therapy_radiation_initial'] = 0
curated1[9, 'therapy_radiation_salvage'] = 0
curated1[9, 'therapy_hormonal_initial'] = 0
curated1[9, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[9,'Prostate Cancer patient 03-159 Gleason_Score:ch1'], ' '))[4],
                                     unlist(stringr::str_split(uncurated1[9,'Prostate Cancer patient 03-159 Gleason_Score:ch1'], ' '))[6], 
                                     unlist(stringr::str_split(uncurated1[9,'Prostate Cancer patient 03-159 Gleason_Score:ch1'], ' '))[2], 
                                     sep ='|')
curated1[9, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[9,'Prostate Cancer patient 03-159 Gleason_Score:ch1'], ' '))[3] , 17, 17))

curated1[10,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[10,'Prostate Cancer patient 03-021F Gleason_Score:ch1'], ' '))[5], 5, 8))
curated1[10, 'tumor_margins_positive'] = 1
curated1[10, 'other_treatment'] = substr(unlist(stringr::str_split(uncurated1[10,'Prostate Cancer patient 03-021F Gleason_Score:ch1'], ' '))[8], 11, 20)
curated1[10, 'therapy_radiation_initial'] = 0
curated1[10, 'therapy_radiation_salvage'] = 0
curated1[10, 'therapy_hormonal_initial'] = 0
curated1[10, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[10,'Prostate Cancer patient 03-021F Gleason_Score:ch1'], ' '))[4], 
                                      unlist(stringr::str_split(uncurated1[10,'Prostate Cancer patient 03-021F Gleason_Score:ch1'], ' '))[6],
                                      unlist(stringr::str_split(uncurated1[10,'Prostate Cancer patient 03-021F Gleason_Score:ch1'], ' '))[2],
                                       sep = '|')
curated1[10, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[10,'Prostate Cancer patient 03-021F Gleason_Score:ch1'], ' '))[3], 17, 17))

curated1[11,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[11,'Prostate Cancer patient 02-003E Gleason_Score:ch2'], ' '))[5], 5, 8))
curated1[11, 'tumor_margins_positive'] = 0
curated1[11, 'therapy_radiation_initial'] = 0
curated1[11, 'therapy_radiation_salvage'] = 0
curated1[11, 'therapy_hormonal_initial'] = 0
curated1[11, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[11,'Prostate Cancer patient 02-003E Gleason_Score:ch2'], ' '))[4], 
                                      unlist(stringr::str_split(uncurated1[11,'Prostate Cancer patient 02-003E Gleason_Score:ch2'], ' '))[6], 
                                      unlist(stringr::str_split(uncurated1[11,'Prostate Cancer patient 02-003E Gleason_Score:ch2'], ' '))[2], 
                                      sep = '|')
curated1[11, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[11,'Prostate Cancer patient 02-003E Gleason_Score:ch2'], ' '))[3], 17, 17))

curated1[12,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[12,'Prostate Cancer patient 03-207C Gleason_Score:ch1'], ' '))[8], 5, 8))
curated1[12, 'tumor_margins_positive'] = 0
curated1[12, 'other_treatment'] = paste(substr(unlist(stringr::str_split(uncurated1[12,'Prostate Cancer patient 03-207C Gleason_Score:ch1'], ' '))[11], 11, 15), unlist(stringr::str_split(uncurated1[12,'Prostate Cancer patient 03-207C Gleason_Score:ch1'], ' '))[12], sep = '_')
curated1[12, 'therapy_radiation_initial'] = 0
curated1[12, 'therapy_radiation_salvage'] = 0
curated1[12, 'therapy_hormonal_initial'] = 0
curated1[12, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[12,'Prostate Cancer patient 03-207C Gleason_Score:ch1'], ' '))[7], 
                                      unlist(stringr::str_split(uncurated1[12,'Prostate Cancer patient 03-207C Gleason_Score:ch1'], ' '))[9], 
                                      unlist(stringr::str_split(uncurated1[12,'Prostate Cancer patient 03-207C Gleason_Score:ch1'], ' '))[5], 
                                      sep = '|')
curated1[12, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[12,'Prostate Cancer patient 03-207C Gleason_Score:ch1'], ' '))[6], 17, 17)) 

curated1[13,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[13,'LCM_Gleason_Pattern:ch2'], ' '))[4], 5, 8))
curated1[13, 'tumor_margins_positive'] = 1
curated1[13, 'therapy_radiation_initial'] = 0
curated1[13, 'therapy_radiation_salvage'] = 0
curated1[13, 'therapy_hormonal_initial'] = 0
curated1[13, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[13,'LCM_Gleason_Pattern:ch2'], ' '))[3],
                                      unlist(stringr::str_split(uncurated1[13,'LCM_Gleason_Pattern:ch2'], ' '))[5], 
                                      unlist(stringr::str_split(uncurated1[13,'characteristics_ch2.1'], ' '))[1], 
                                       sep = '|')
curated1[13, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[13,'LCM_Gleason_Pattern:ch2'], ' '))[2], 17, 17))

curated1[14,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[14,' LCM_Gleason_Pattern:ch2'], ' '))[4], 5, 8))
curated1[14, 'tumor_margins_positive'] = 0
curated1[14, 'therapy_radiation_initial'] = 0
curated1[14, 'therapy_radiation_salvage'] = 0
curated1[14, 'therapy_hormonal_initial'] = 0
curated1[14, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[14,' LCM_Gleason_Pattern:ch2'], ' '))[3], 
                                      unlist(stringr::str_split(uncurated1[14,' LCM_Gleason_Pattern:ch2'], ' '))[5], 
                                      unlist(stringr::str_split(uncurated1[14,'characteristics_ch2.1'], ' '))[2], 
                                      sep = '|')
curated1[14, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[14,' LCM_Gleason_Pattern:ch2'], ' '))[2], 17, 17))

curated1[15,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[15,' LCM_Gleason_Pattern:ch1'], ' '))[4], 5, 8))
curated1[15, 'tumor_margins_positive'] = 1
curated1[15, 'therapy_radiation_initial'] = 0
curated1[15, 'therapy_radiation_salvage'] = 0
curated1[15, 'therapy_hormonal_initial'] = 0
curated1[15, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[15,' LCM_Gleason_Pattern:ch1'], ' '))[3], 
                                      unlist(stringr::str_split(uncurated1[15,' LCM_Gleason_Pattern:ch1'], ' '))[5],
                                      unlist(stringr::str_split(uncurated1[15,'characteristics_ch1.1'], ' '))[2],
                                      sep = '|')
curated1[15, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[15,' LCM_Gleason_Pattern:ch1'], ' '))[2], 17, 17))

curated1[16,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[16,'Prostate Cancer patient 03-068D Gleason_Score:ch2'], ' '))[8], 5, 8))
curated1[16, 'tumor_margins_positive'] = 1
curated1[16, 'therapy_radiation_initial'] = 0
curated1[16, 'therapy_radiation_salvage'] = 0
curated1[16, 'therapy_hormonal_initial'] = 0
curated1[16, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[16,'Prostate Cancer patient 03-068D Gleason_Score:ch2'], ' '))[7],
                                     unlist(stringr::str_split(uncurated1[16,'Prostate Cancer patient 03-068D Gleason_Score:ch2'], ' '))[9],
                                     unlist(stringr::str_split(uncurated1[16,'Prostate Cancer patient 03-068D Gleason_Score:ch2'], ' '))[5],
                                     sep = '|')
curated1[16, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[16,'Prostate Cancer patient 03-068D Gleason_Score:ch2'], ' '))[6], 17, 17))

curated1[17,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[17,'Prostate Cancer patient 03-141C Gleason_Score:ch2'], ' '))[5], 5, 8))
curated1[17, 'tumor_margins_positive'] = 1
curated1[17, 'therapy_radiation_initial'] = 0
curated1[17, 'therapy_radiation_salvage'] = 0
curated1[17, 'therapy_hormonal_initial'] = 0
curated1[17, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[17,'Prostate Cancer patient 03-141C Gleason_Score:ch2'], ' '))[4], 
                                      unlist(stringr::str_split(uncurated1[17,'Prostate Cancer patient 03-141C Gleason_Score:ch2'], ' '))[6],
                                      unlist(stringr::str_split(uncurated1[17,'Prostate Cancer patient 03-141C Gleason_Score:ch2'], ' '))[2],
                                      sep = '|')
curated1[17, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[17,'Prostate Cancer patient 03-141C Gleason_Score:ch2'], ' '))[3], 17, 17))

curated1[18,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[18,'Prostate Cancer patient 02-053C Gleason_Score:ch1'], ' '))[5], 5, 8))
curated1[18, 'tumor_margins_positive'] = 0
curated1[18, 'therapy_radiation_initial'] = 0
curated1[18, 'therapy_radiation_salvage'] = 0
curated1[18, 'therapy_hormonal_initial'] = 0
curated1[18, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[18,'Prostate Cancer patient 02-053C Gleason_Score:ch1'], ' '))[4], 
                                      unlist(stringr::str_split(uncurated1[18,'Prostate Cancer patient 02-053C Gleason_Score:ch1'], ' '))[6],
                                      unlist(stringr::str_split(uncurated1[18,'Prostate Cancer patient 02-053C Gleason_Score:ch1'], ' '))[2],
                                      sep = '|')
curated1[18, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[18,'Prostate Cancer patient 02-053C Gleason_Score:ch1'], ' '))[3], 17, 17))

curated1[19,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[19,'Prostate Cancer patient 04-030A Gleason_Score:ch2'], ' '))[5], 5, 8))
curated1[19, 'tumor_margins_positive'] = 0
curated1[19, 'therapy_radiation_initial'] = 0
curated1[19, 'therapy_radiation_salvage'] = 0
curated1[19, 'therapy_hormonal_initial'] = 0
curated1[19, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[19,'Prostate Cancer patient 04-030A Gleason_Score:ch2'], ' '))[4], 
                                      unlist(stringr::str_split(uncurated1[19,'Prostate Cancer patient 04-030A Gleason_Score:ch2'], ' '))[6],
                                      unlist(stringr::str_split(uncurated1[19,'Prostate Cancer patient 04-030A Gleason_Score:ch2'], ' '))[2],
                                      sep = '|')
curated1[19, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[19,'Prostate Cancer patient 04-030A Gleason_Score:ch2'], ' '))[3], 17, 17))

curated1[20,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[20,' LCM_Gleason_Pattern:ch1'], ' '))[4], 5, 8))
curated1[20, 'tumor_margins_positive'] = 1
curated1[20, 'therapy_radiation_initial'] = 0
curated1[20, 'therapy_radiation_salvage'] = 0
curated1[20, 'therapy_hormonal_initial'] = 0
curated1[20, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[20,' LCM_Gleason_Pattern:ch1'], ' '))[3], 
                                      unlist(stringr::str_split(uncurated1[20,' LCM_Gleason_Pattern:ch1'], ' '))[5], 
                                      unlist(stringr::str_split(uncurated1[20,'characteristics_ch1.1'], ' '))[2], 
                                      sep = '|')
curated1[20, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[20,' LCM_Gleason_Pattern:ch1'], ' '))[2] , 17, 17))

curated1[21,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[21,'LCM_Gleason_Pattern:ch2'], ' '))[4], 5, 8))
curated1[21, 'tumor_margins_positive'] = 0
curated1[21, 'other_treatment'] = substr(unlist(stringr::str_split(uncurated1[21,'LCM_Gleason_Pattern:ch2'], ' '))[7], 11, 20)
curated1[21, 'therapy_radiation_initial'] = 0
curated1[21, 'therapy_radiation_salvage'] = 0
curated1[21, 'therapy_hormonal_initial'] = 0
curated1[21, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[21,'LCM_Gleason_Pattern:ch2'], ' '))[3], 
                                      unlist(stringr::str_split(uncurated1[21,'LCM_Gleason_Pattern:ch2'], ' '))[5], 
                                      unlist(stringr::str_split(uncurated1[21,'characteristics_ch2.1'], ' '))[1], 
                                      sep = '|')
curated1[21, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[21,'LCM_Gleason_Pattern:ch2'], ' '))[2], 17, 17))

curated1[22,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[22,'Prostate Cancer patient 03-184C Gleason_Score:ch2'], ' '))[5], 5, 8))
curated1[22, 'tumor_margins_positive'] = 0
curated1[22, 'therapy_radiation_initial'] = 0
curated1[22, 'therapy_radiation_salvage'] = 0
curated1[22, 'therapy_hormonal_initial'] = 0
curated1[22, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[22,'Prostate Cancer patient 03-184C Gleason_Score:ch2'], ' '))[4],
                                      unlist(stringr::str_split(uncurated1[22,'Prostate Cancer patient 03-184C Gleason_Score:ch2'], ' '))[6],
                                      unlist(stringr::str_split(uncurated1[22,'Prostate Cancer patient 03-184C Gleason_Score:ch2'], ' '))[2],
                                      sep = '|')
curated1[22, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[22,'Prostate Cancer patient 03-184C Gleason_Score:ch2'], ' '))[3], 17, 17))

curated1[23,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[23,'Prostate Cancer patient 03-152A Gleason_Score:ch2'], ' '))[5], 5, 8))
curated1[23, 'tumor_margins_positive'] = 1
curated1[23, 'other_treatment'] = paste(substr(unlist(stringr::str_split(uncurated1[23,'Prostate Cancer patient 03-152A Gleason_Score:ch2'], ' '))[8], 11, 20), unlist(stringr::str_split(uncurated1[23,'Prostate Cancer patient 03-152A Gleason_Score:ch2'], ' '))[9], sep = '_')
curated1[23, 'therapy_radiation_initial'] = 0
curated1[23, 'therapy_radiation_salvage'] = 0
curated1[23, 'therapy_hormonal_initial'] = 0
curated1[23, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[23,'Prostate Cancer patient 03-152A Gleason_Score:ch2'], ' '))[4], 
                                      unlist(stringr::str_split(uncurated1[23,'Prostate Cancer patient 03-152A Gleason_Score:ch2'], ' '))[6],
                                      unlist(stringr::str_split(uncurated1[23,'Prostate Cancer patient 03-152A Gleason_Score:ch2'], ' '))[2],
                                      sep = '|')
curated1[23, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[23,'Prostate Cancer patient 03-152A Gleason_Score:ch2'], ' '))[3], 17, 17))

curated1[24,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[24,'Prostate Cancer patient 03-029B Gleason_Score:ch1'], ' '))[5], 5, 8))
curated1[24, 'tumor_margins_positive'] = 0
curated1[24, 'therapy_radiation_initial'] = 0
curated1[24, 'therapy_radiation_salvage'] = 0
curated1[24, 'therapy_hormonal_initial'] = 0
curated1[24, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[24,'Prostate Cancer patient 03-029B Gleason_Score:ch1'], ' '))[4], 
                                      unlist(stringr::str_split(uncurated1[24,'Prostate Cancer patient 03-029B Gleason_Score:ch1'], ' '))[6],
                                      unlist(stringr::str_split(uncurated1[24,'Prostate Cancer patient 03-029B Gleason_Score:ch1'], ' '))[2],
                                      sep = '|')
curated1[24, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[24,'Prostate Cancer patient 03-029B Gleason_Score:ch1'], ' '))[3], 17, 17))

curated1[25,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[25,'Prostate Cancer patient 03-015F Gleason_Score:ch2'], ' '))[5], 5, 8))
curated1[25, 'tumor_margins_positive'] = 1
curated1[25, 'therapy_radiation_initial'] = 0
curated1[25, 'therapy_radiation_salvage'] = 0
curated1[25, 'therapy_hormonal_initial'] = 0
curated1[25, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[25,'Prostate Cancer patient 03-015F Gleason_Score:ch2'], ' '))[4],
                                      unlist(stringr::str_split(uncurated1[25,'Prostate Cancer patient 03-015F Gleason_Score:ch2'], ' '))[6],
                                      unlist(stringr::str_split(uncurated1[25,'Prostate Cancer patient 03-015F Gleason_Score:ch2'], ' '))[2],
                                      sep = '|')
curated1[25, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[25,'Prostate Cancer patient 03-015F Gleason_Score:ch2'], ' '))[3], 17, 17))

curated1[26,'psa'] = NA
curated1[26, 'tumor_margins_positive'] = 1
curated1[26, 'therapy_radiation_initial'] = 0
curated1[26, 'therapy_radiation_salvage'] = 0
curated1[26, 'therapy_hormonal_initial'] = 0
curated1[26, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[26,'Prostate Cancer patient 03-119D Gleason_Score:ch1'], ' '))[4], 
                                      unlist(stringr::str_split(uncurated1[26,'Prostate Cancer patient 03-119D Gleason_Score:ch1'], ' '))[6],
                                      unlist(stringr::str_split(uncurated1[26,'Prostate Cancer patient 03-119D Gleason_Score:ch1'], ' '))[2],
                                      sep = '|')
curated1[26, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[26,'Prostate Cancer patient 03-119D Gleason_Score:ch1'], ' '))[3], 17, 17))

curated1[27,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[27,'Prostate Cancer patient 03-140B Gleason_Score:ch1'], ' '))[5], 5, 8))
curated1[27, 'tumor_margins_positive'] = 0
curated1[27, 'other_treatment'] = paste(substr(unlist(stringr::str_split(uncurated1[27,'Prostate Cancer patient 03-140B Gleason_Score:ch1'], ' '))[8], 11, 20), unlist(stringr::str_split(uncurated1[27,'Prostate Cancer patient 03-140B Gleason_Score:ch1'], ' '))[9], sep = '_')
curated1[27, 'therapy_radiation_initial'] = 0
curated1[27, 'therapy_radiation_salvage'] = 0
curated1[27, 'therapy_hormonal_initial'] = 0
curated1[27, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[27,'Prostate Cancer patient 03-140B Gleason_Score:ch1'], ' '))[4],
                                      unlist(stringr::str_split(uncurated1[27,'Prostate Cancer patient 03-140B Gleason_Score:ch1'], ' '))[6],
                                      unlist(stringr::str_split(uncurated1[27,'Prostate Cancer patient 03-140B Gleason_Score:ch1'], ' '))[2],
                                      sep = '|')
curated1[27, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[27,'Prostate Cancer patient 03-140B Gleason_Score:ch1'], ' '))[3], 17, 17))

curated1[28,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[28,'LCM_Gleason_Pattern:ch1'], ' '))[4], 5, 8))
curated1[28, 'tumor_margins_positive'] = 0
curated1[28, 'therapy_radiation_initial'] = 0
curated1[28, 'therapy_radiation_salvage'] = 0
curated1[28, 'therapy_hormonal_initial'] = 0
curated1[28, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[28,'LCM_Gleason_Pattern:ch1'], ' '))[3],
                                      unlist(stringr::str_split(uncurated1[28,'LCM_Gleason_Pattern:ch1'], ' '))[5],
                                      unlist(stringr::str_split(uncurated1[28,'characteristics_ch1.1'], ' '))[1],
                                      sep = '|')
curated1[28, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[28,'LCM_Gleason_Pattern:ch1'], ' '))[2], 17, 17))

curated1[29,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[29,' LCM_Gleason_Pattern:ch2'], ' '))[4], 5, 8))
curated1[29, 'tumor_margins_positive'] = 0
curated1[29, 'therapy_radiation_initial'] = 0
curated1[29, 'therapy_radiation_salvage'] = 0
curated1[29, 'therapy_hormonal_initial'] = 0
curated1[29, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[29,' LCM_Gleason_Pattern:ch2'], ' '))[3], 
                                      unlist(stringr::str_split(uncurated1[29,' LCM_Gleason_Pattern:ch2'], ' '))[5],
                                      unlist(stringr::str_split(uncurated1[29,'characteristics_ch2.1'], ' '))[2],
                                      sep = '|')
curated1[29, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[29,' LCM_Gleason_Pattern:ch2'], ' '))[2], 17, 17))

curated1[30,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[30,'Prostate Cancer patient 03-063A Gleason_Score:ch2'], ' '))[5], 5, 8))
curated1[30, 'tumor_margins_positive'] = 0
curated1[30, 'therapy_radiation_initial'] = 0
curated1[30, 'therapy_radiation_salvage'] = 0
curated1[30, 'therapy_hormonal_initial'] = 0
curated1[30, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[30,'Prostate Cancer patient 03-063A Gleason_Score:ch2'], ' '))[4], 
                                      unlist(stringr::str_split(uncurated1[30,'Prostate Cancer patient 03-063A Gleason_Score:ch2'], ' '))[6],
                                      unlist(stringr::str_split(uncurated1[30,'Prostate Cancer patient 03-063A Gleason_Score:ch2'], ' '))[2],
                                      sep = '|')
curated1[30, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[30,'Prostate Cancer patient 03-063A Gleason_Score:ch2'], ' '))[3],17, 17))

curated1[31,'psa'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[31,'Prostate Cancer patient 03-055H Gleason_Score:ch1'], ' '))[5], 5, 8))
curated1[31, 'tumor_margins_positive'] = 0
curated1[31, 'other_treatment'] = paste(substr(unlist(stringr::str_split(uncurated1[31,'Prostate Cancer patient 03-055H Gleason_Score:ch1'], ' '))[8], 11, 20), unlist(stringr::str_split(uncurated1[31,'Prostate Cancer patient 03-055H Gleason_Score:ch1'], ' '))[9], sep = '_')
curated1[31, 'therapy_radiation_initial'] = 0
curated1[31, 'therapy_radiation_salvage'] = 0
curated1[31, 'therapy_hormonal_initial'] = 0
curated1[31, 'other_feature'] = paste(unlist(stringr::str_split(uncurated1[31,'Prostate Cancer patient 03-055H Gleason_Score:ch1'], ' '))[4],
                                      unlist(stringr::str_split(uncurated1[31,'Prostate Cancer patient 03-055H Gleason_Score:ch1'], ' '))[6],
                                      unlist(stringr::str_split(uncurated1[31,'Prostate Cancer patient 03-055H Gleason_Score:ch1'], ' '))[2],
                                      sep = '|')
curated1[31, 'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated1[31,'Prostate Cancer patient 03-055H Gleason_Score:ch1'], ' '))[3], 17, 17))

###########################
###########################

# being only one line the second dataset does not cause such inconvenients

curated2 <- curated2 %>% 
  dplyr::mutate(study_name = "True et al.") %>%
  dplyr::mutate(sample_name = uncurated2$geo_accession) %>% 
  dplyr::mutate(sample_paired = 1) %>%
  dplyr::mutate(sample_type = "primary") %>%
  dplyr::mutate(patient_id = unlist(lapply(stringr::str_split(uncurated2$title, '_'), function(x) x[2]))) %>%
  dplyr::mutate(gleason_major = as.numeric(substr(unlist(lapply(stringr::str_split(uncurated2$description, ' '), function(x) x[9])), 1, 1))) %>%
  dplyr::mutate(gleason_minor = as.numeric(substr(unlist(lapply(stringr::str_split(uncurated2$description, ' '), function(x) x[9])), 3, 3))) %>%
  dplyr::mutate(gleason_grade = gleason_major + gleason_minor) %>%
  dplyr::mutate(gleason_group = dplyr::case_when(
    gleason_grade  %in%  4:6 ~ "<=6",
    gleason_major == 3 & gleason_minor == 4 ~ "3+4",
    gleason_major == 4 & gleason_minor == 3 ~ "4+3",
    gleason_grade %in% 8:10 ~ ">=8",
    TRUE ~ "NA"
  )) %>%
  dplyr::mutate(psa = as.numeric(substr(unlist(stringr::str_split(uncurated2$characteristics_ch1, ' '))[13], 5, 8))) %>%
  dplyr::mutate(tumor_margins_positive = 0) %>%
  dplyr::mutate(therapy_radiation_initial = 0) %>%
  dplyr::mutate(therapy_radiation_salvage = 0) %>%
  dplyr::mutate(therapy_hormonal_initial = 0) %>%
  dplyr::mutate(microdissected = 1) %>%
  dplyr::mutate(other_feature = paste(
                                      unlist(stringr::str_split(uncurated2$characteristics_ch1, ' '))[12],
                                      unlist(stringr::str_split(uncurated2$characteristics_ch1, ' '))[14],
                                      unlist(stringr::str_split(uncurated2$characteristics_ch1, ' '))[10],
                                      paste('tertiary_pattern', stringr::str_split(uncurated2$characteristics_ch1, ' ')[[1]][8], sep = ':'),
                                      sep = '|' )) 

curated2[1,'grade_group'] = as.numeric(substr(unlist(stringr::str_split(uncurated2[1,"Prostate Cancer patient 02-209C Gleason_Score:ch1"], ' '))[7], 17, 17))

# the datastes are merged in a way that preserves the order of the GEO sample indexes

clinical_true = rbind(curated1[1:10, ], curated2[1, ], curated1[11:31, ])

save(clinical_true, file =  "./data-raw//clinical_true.RData")


#############################################################################
#############################################################################

##  WALLACE ET AL.

#############################################################################
#############################################################################

# load series and platform data from GEO

gset <- getGEO("GSE6956", GSEMatrix =TRUE, getGPL=FALSE)

# clinical
uncurated <- Biobase::pData(gset[[1]]) 

unmatched_healty_tissue = c('GSM160418', 'GSM160419', 'GSM160420', 'GSM160421', 'GSM160422', 'GSM160430') # as determined by the GSE6956 metadata

uncurated = uncurated[!is.element(uncurated$geo_accession, unmatched_healty_tissue), ] 

curated <- initial_curated_df(
  df_rownames = rownames(uncurated),
  template_name="data-raw/template_prad.csv")


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
  dplyr::mutate(sample_paired = 0) %>%                                                           
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

save(clinical_wallace, file = "data-raw/clinical_wallace.RData")


########################################################################
########################################################################
#
# Weiner et al.
#
########################################################################
########################################################################

gset <- getGEO("GSE157548", GSEMatrix =TRUE, getGPL=TRUE)


uncurated <- Biobase::pData(gset[[1]]) 

curated <- initial_curated_df(
  df_rownames = rownames(uncurated),
  template_name="data-raw/template_prad.csv")

curated <- curated %>% 
  dplyr::mutate(study_name = "Weiner et al.") %>%
  dplyr::mutate(sample_name = uncurated$geo_accession) %>% 
  dplyr::mutate(patient_id = uncurated$geo_accession) %>%
  dplyr::mutate(race = dplyr::case_when(
                                        is.na(uncurated$'race:ch1') ~ NA_character_, 
                                        uncurated$'race:ch1' == 'Black' ~ 'african_american',
                                        uncurated$'race:ch1' == 'White' ~ 'caucasian')) %>% 
  dplyr::mutate(psa = dplyr::case_when(uncurated$"psa:ch1" == "NA" ~ NA_character_,
                                       TRUE ~ uncurated$"psa:ch1"
                                       )) %>% 
  dplyr::mutate(psa = as.numeric(psa)) %>%
  dplyr::mutate(tissue_source = 'prostatectomy') %>%
  dplyr::mutate(age_at_initial_diagnosis = dplyr::case_when(
                                       is.na(uncurated$'age:ch1') ~ NA_character_, 
                                       TRUE ~ uncurated$'age:ch1')) %>% 
  dplyr::mutate(age_at_initial_diagnosis = as.numeric(age_at_initial_diagnosis)) %>%
  dplyr::mutate(grade_group = dplyr::case_when(
                                               uncurated$"grade group:ch1" == "1" ~ "<=6",
                                               uncurated$"grade group:ch1" == "2" ~ "3+4",
                                               uncurated$"grade group:ch1" == "3" ~ "4+3",
                                               uncurated$"grade group:ch1" == "4" ~ ">=8",
                                               uncurated$"grade group:ch1" == "5" ~ ">=8",
                                               )) %>%
  dplyr::mutate(frozen_ffpe = 'FFPE') %>%
  dplyr::mutate(batch = dplyr::case_when(
                                           is.na(uncurated$"race:ch1") ~ 'Durham Veterans Affairs Hospital',
                                           uncurated$"race:ch1" == 'Black' ~ 'Johns Hopkins Medical Institute',
                                           uncurated$"race:ch1" == 'White' ~ 'Johns Hopkins Medical Institute'
                                           )) %>%
  dplyr::mutate(therapy_radiation_initial = 0) %>%
  dplyr::mutate(therapy_radiation_salvage = 0) %>%
  dplyr::mutate(therapy_surgery_initial = 0) %>%
  dplyr::mutate(therapy_hormonal_initial = 0) %>%
  dplyr::mutate(sample_paired = 0) %>%
  dplyr::mutate(sample_type = 'primary') %>%
  dplyr::mutate(microdissected = 0) %>%
  dplyr::mutate(grade_group1 = uncurated$"grade group:ch1") %>%
  dplyr::mutate(grade_group1 = as.numeric(grade_group1))

clinical_weiner <- curated

save(clinical_weiner, file = 'data-raw/clinical_weiner.RData')



#######################################################################
#Wang et al
######################################################################

gse <- GEOquery::getGEO("GSE8218", GSEMatrix = TRUE)

uncurated <- Biobase::pData(gse[[1]])

curated <- initial_curated_df(
  df_rownames = rownames(uncurated),
  template_name="data-raw/template_prad.csv")

curated <- curated %>% 
  dplyr::mutate(study_name = "Wang et al.") %>%
  dplyr::mutate(patient_id = rownames(uncurated)) %>%
  dplyr::mutate(sample_name = row.names(uncurated)) %>% 
  dplyr::mutate(sample_paired = 0) %>%
  dplyr::mutate(sample_type = "primary") %>%
  ## TDL: Tumor purity of percentage is a key field rather than specific subsets
  # From GEO description: "The percentage of different cell types vary considerably among samples and were determined by pathologist."
  # Despite reported as 'prostate cancer sample(s)' in GEO, the samples are quite mixed.
  dplyr::mutate(tumor_purity_pathology = as.numeric(gsub("%", "", uncurated$"Percentage of Tumor:ch1"))/100)
  #dplyr::mutate(Percentage_of_Atrophic_Gland = uncurated$`Percentage of Atrophic Gland:ch1`) %>%
  #dplyr::mutate(Percentage_of_BPH = uncurated$`Percentage of BPH:ch1`) %>%
  #dplyr::mutate(Percentage_of_Stroma = uncurated$`Percentage of Stroma:ch1`) %>%
  #dplyr::mutate(Percentage_of_Tumor = uncurated$`Percentage of Tumor:ch1`)
  
clinical_wang <- curated

save(clinical_wang, file = "data-raw/clinical_wang.RData")

#######################################################################
# IGC - GSE2109
######################################################################
gse <- GEOquery::getGEO("GSE2109", GSEMatrix = TRUE)

uncurated <- Biobase::pData(gse[[1]])
uncurated <- uncurated[grep("Prostate", uncurated$title), ]

curated <- initial_curated_df(
  df_rownames = rownames(uncurated),
  template_name="data-raw/template_prad.csv")
#uncurated <- uncurated$title

# Grep row-wise interesting fields
uncurated$PathologicalM <- unlist(apply(uncurated, MARGIN=1, FUN=function(x) { tmp <- gsub("Pathological M: ", "", grep("Pathological M: ", x, value=TRUE)); ifelse(length(tmp)>0, tmp, NA) }))
uncurated$PathologicalN <- unlist(apply(uncurated, MARGIN=1, FUN=function(x) { tmp <- gsub("Pathological N: ", "", grep("Pathological N: ", x, value=TRUE)); ifelse(length(tmp)>0, tmp, NA) }))
uncurated$PrimarySite <- unlist(apply(uncurated, MARGIN=1, FUN=function(x) { tmp <- gsub("Primary Site: ", "", grep("Primary Site: ", x, value=TRUE)); ifelse(length(tmp)>0, tmp, NA) }))
uncurated$Histology <- unlist(apply(uncurated, MARGIN=1, FUN=function(x) { tmp <- gsub("Histology: ", "", grep("Histology: ", x, value=TRUE)); ifelse(length(tmp)>0, tmp, NA) }))
uncurated$PathologicalGleasonScore <- unlist(apply(uncurated, MARGIN=1, FUN=function(x) { tmp <- gsub("Pathological Gleason Score: ", "", grep("Pathological Gleason Score: ", x, value=TRUE)); ifelse(length(tmp)>0, tmp, NA) }))
uncurated$ClinicalGleasonScore <- unlist(apply(uncurated, MARGIN=1, FUN=function(x) { tmp <- gsub("Clinical Gleason Score: ", "", grep("Clinical Gleason Score: ", x, value=TRUE)); ifelse(length(tmp)>0, tmp, NA) }))

curated <- curated %>% 
  dplyr::mutate(study_name = "igc") %>%
  dplyr::mutate(sample_name = rownames(uncurated)) %>% 
  dplyr::mutate(sample_paired = 0) %>%
  dplyr::mutate(patient_id = rownames(uncurated)) %>%
  dplyr::mutate(age_at_initial_diagnosis = stringr::str_remove(uncurated$description.1,"Patient Age:")) %>%
  dplyr::mutate(race = stringr::str_remove(uncurated$description.3,"Ethnic Background:")) %>%
  dplyr::mutate(description.6 = uncurated$description.6)%>%
  dplyr::mutate(smoking_status = dplyr::case_when(
    description.6 %in% c('Type of Tobacco Use: Cigarettes','Type of Tobacco Use: Cigars') ~ 1,
    TRUE ~ 0)) %>%
  
  dplyr::mutate(psa_category = dplyr::case_when(
    uncurated$description.6 %in% 'PSA: Normal' ~ 'Normal',
    uncurated$description.6 %in% 'PSA: Elevated' ~ 'Elevated',
    uncurated$description.7 %in% 'PSA: Normal' ~ 'Normal',
    uncurated$description.7 %in% 'PSA: Elevated' ~ 'Elevated',
    uncurated$description.8 %in% 'PSA: Normal' ~ 'Normal',
    uncurated$description.8 %in% 'PSA: Elevated' ~ 'Elevated',
    uncurated$description.9 %in% 'PSA: Normal' ~ 'Normal',
    uncurated$description.9 %in% 'PSA: Elevated' ~ 'Elevated',
    uncurated$description.10 %in% 'PSA: Elevated' ~ 'Elevated',
    uncurated$description.10 %in% 'PSA: Normal' ~ 'Normal',
    uncurated$description.11 %in% 'PSA: Normal' ~ 'Normal',
    uncurated$description.11 %in% 'PSA: Elevated' ~ 'Elevated',
    uncurated$description.12 %in% 'PSA: Elevated' ~ 'Elevated',
    uncurated$description.12 %in% 'PSA: Normal' ~ 'Normal',
    uncurated$description.13 %in% 'PSA: Normal' ~ 'Normal',
    uncurated$description.13 %in% 'PSA: Elevated' ~ 'Elevated',
    uncurated$description.15 %in% 'PSA: Normal' ~ 'Normal',
    uncurated$description.15 %in% 'PSA: Elevated' ~ 'Elevated'
    )) %>%
  dplyr::mutate(T_pathological = dplyr::case_when(
    uncurated$description.7 %in% c('Pathological T: 2a','Pathological T: 2b', 'Pathological T: 2c') ~ 2, 
    uncurated$description.7 %in% c('Pathological T: 3a','Pathological T: 3b', 'Pathological T: 3c') ~ 3,
    uncurated$description.7 %in% c('Pathological T: 4', 'Pathological T: 4a') ~ 4,
    uncurated$description.8 %in% c('Pathological T: 2a','Pathological T: 2b', 'Pathological T: 2c') ~ 2,
    uncurated$description.8 %in% c('Pathological T: 3a','Pathological T: 3b', 'Pathological T: 3c') ~ 3,
    uncurated$description.8 %in% c('Pathological T: 4', 'Pathological T: 4a') ~ 4,
    uncurated$description.9 %in% c('Pathological T: 2a','Pathological T: 2b', 'Pathological T: 2c') ~ 2,
    uncurated$description.9 %in% c('Pathological T: 3a','Pathological T: 3b', 'Pathological T: 3c') ~ 3,
    uncurated$description.9 %in% c('Pathological T: 4', 'Pathological T: 4a') ~ 4,
    uncurated$description.10 %in% c('Pathological T: 2a','Pathological T: 2b', 'Pathological T: 2c') ~ 2,
    uncurated$description.10 %in% c('Pathological T: 3a','Pathological T: 3b', 'Pathological T: 3c') ~ 3,
    uncurated$description.10 %in% c('Pathological T: 4', 'Pathological T: 4a') ~ 4,
    uncurated$description.11 %in% c('Pathological T: 2a','Pathological T: 2b', 'Pathological T: 2c') ~ 2,
    uncurated$description.11 %in% c('Pathological T: 3a','Pathological T: 3b', 'Pathological T: 3c') ~ 3,
    uncurated$description.11 %in% c('Pathological T: 4', 'Pathological T: 4a') ~ 4,
    uncurated$description.12 %in% c('Pathological T: 2a','Pathological T: 2b', 'Pathological T: 2c') ~ 2,
    uncurated$description.12 %in% c('Pathological T: 3a','Pathological T: 3b', 'Pathological T: 3c') ~ 3,
    uncurated$description.12 %in% c('Pathological T: 4', 'Pathological T: 4a') ~ 4,
    uncurated$description.13 %in% c('Pathological T: 2a','Pathological T: 2b', 'Pathological T: 2c') ~ 2,
    uncurated$description.13 %in% c('Pathological T: 3a','Pathological T: 3b', 'Pathological T: 3c') ~ 3,
    uncurated$description.13 %in% c('Pathological T: 4', 'Pathological T: 4a') ~ 4,
    uncurated$description.14 %in% c('Pathological T: 2a','Pathological T: 2b', 'Pathological T: 2c') ~ 2,
    uncurated$description.14 %in% c('Pathological T: 3a','Pathological T: 3b', 'Pathological T: 3c') ~ 3,
    uncurated$description.14 %in% c('Pathological T: 4', 'Pathological T: 4a') ~ 4,
    uncurated$description.15 %in% c('Pathological T: 2a','Pathological T: 2b', 'Pathological T: 2c') ~ 2,
    uncurated$description.15 %in% c('Pathological T: 3a','Pathological T: 3b', 'Pathological T: 3c') ~ 3,
    uncurated$description.15 %in% c('Pathological T: 4', 'Pathological T: 4a') ~ 4,
    uncurated$description.17 %in% c('Pathological T: 2a','Pathological T: 2b', 'Pathological T: 2c') ~ 2,
    uncurated$description.17 %in% c('Pathological T: 3a','Pathological T: 3b', 'Pathological T: 3c') ~ 3,
    uncurated$description.17 %in% c('Pathological T: 4', 'Pathological T: 4a') ~ 4,
    uncurated$description.19 %in% c('Pathological T: 2a','Pathological T: 2b', 'Pathological T: 2c') ~ 2,
    uncurated$description.19 %in% c('Pathological T: 3a','Pathological T: 3b', 'Pathological T: 3c') ~ 3,
    uncurated$description.19 %in% c('Pathological T: 4', 'Pathological T: 4a') ~ 4 )) %>%
  
  dplyr::mutate(T_substage_pathological = dplyr::case_when(
    uncurated$description.7 %in% c('Pathological T: 2a','Pathological T: 3a', 'Pathological T: 4a') ~ 'a',
    uncurated$description.8 %in% c('Pathological T: 2a','Pathological T: 3a', 'Pathological T: 4a') ~ 'a',
    uncurated$description.9 %in% c('Pathological T: 2a','Pathological T: 3a', 'Pathological T: 4a') ~ 'a',
    uncurated$description.10 %in% c('Pathological T: 2a','Pathological T: 3a', 'Pathological T: 4a') ~ 'a',
    uncurated$description.11 %in% c('Pathological T: 2a','Pathological T: 3a', 'Pathological T: 4a') ~ 'a',
    uncurated$description.12 %in% c('Pathological T: 2a','Pathological T: 3a', 'Pathological T: 4a') ~ 'a',
    uncurated$description.13 %in% c('Pathological T: 2a','Pathological T: 3a', 'Pathological T: 4a') ~ 'a',
    uncurated$description.14 %in% c('Pathological T: 2a','Pathological T: 3a', 'Pathological T: 4a') ~ 'a',
    uncurated$description.15 %in% c('Pathological T: 2a','Pathological T: 3a', 'Pathological T: 4a') ~ 'a',
    uncurated$description.17 %in% c('Pathological T: 2a','Pathological T: 3a', 'Pathological T: 4a') ~ 'a',
    uncurated$description.19 %in% c('Pathological T: 2a','Pathological T: 3a', 'Pathological T: 4a') ~ 'a',
    uncurated$description.7 %in% c('Pathological T: 2b','Pathological T: 3b', 'Pathological T: 4b') ~ 'b',
    uncurated$description.8 %in% c('Pathological T: 2b','Pathological T: 3b', 'Pathological T: 4b') ~ 'b',
    uncurated$description.9 %in% c('Pathological T: 2b','Pathological T: 3b', 'Pathological T: 4b') ~ 'b',
    uncurated$description.10 %in% c('Pathological T: 2b','Pathological T: 3b', 'Pathological T: 4b') ~ 'b',
    uncurated$description.11 %in% c('Pathological T: 2b','Pathological T: 3b', 'Pathological T: 4b') ~ 'b',
    uncurated$description.12 %in% c('Pathological T: 2b','Pathological T: 3b', 'Pathological T: 4b') ~ 'b',
    uncurated$description.13 %in% c('Pathological T: 2b','Pathological T: 3b', 'Pathological T: 4b') ~ 'b',
    uncurated$description.14 %in% c('Pathological T: 2b','Pathological T: 3b', 'Pathological T: 4b') ~ 'b',
    uncurated$description.15 %in% c('Pathological T: 2b','Pathological T: 3b', 'Pathological T: 4b') ~ 'b',
    uncurated$description.17 %in% c('Pathological T: 2b','Pathological T: 3b', 'Pathological T: 4b') ~ 'b',
    uncurated$description.19 %in% c('Pathological T: 2b','Pathological T: 3b', 'Pathological T: 4b') ~ 'b',
    uncurated$description.7 %in% c('Pathological T: 2c','Pathological T: 3c', 'Pathological T: 4c') ~ 'c',
    uncurated$description.8 %in% c('Pathological T: 2c','Pathological T: 3c', 'Pathological T: 4c') ~ 'c',
    uncurated$description.9 %in% c('Pathological T: 2c','Pathological T: 3c', 'Pathological T: 4c') ~ 'c',
    uncurated$description.10 %in% c('Pathological T: 2c','Pathological T: 3c', 'Pathological T: 4c') ~ 'c',
    uncurated$description.11 %in% c('Pathological T: 2c','Pathological T: 3c', 'Pathological T: 4c') ~ 'c',
    uncurated$description.12 %in% c('Pathological T: 2c','Pathological T: 3c', 'Pathological T: 4c') ~ 'c',
    uncurated$description.13 %in% c('Pathological T: 2c','Pathological T: 3c', 'Pathological T: 4c') ~ 'c',
    uncurated$description.14 %in% c('Pathological T: 2c','Pathological T: 3c', 'Pathological T: 4c') ~ 'c',
    uncurated$description.15 %in% c('Pathological T: 2c','Pathological T: 3c', 'Pathological T: 4c') ~ 'c',
    uncurated$description.17 %in% c('Pathological T: 2c','Pathological T: 3c', 'Pathological T: 4c') ~ 'c',
    uncurated$description.19 %in% c('Pathological T: 2c','Pathological T: 3c', 'Pathological T: 4c') ~ 'c')) %>%
  
  dplyr::mutate(gleason_grade = dplyr::case_when(
    uncurated$description.14 %in% 'Clinical Gleason Score: 7' ~ '7',
    uncurated$description.14 %in% 'Clinical Gleason Score: 8-10' ~ '8',
    uncurated$description.14 %in% 'Clinical Gleason Score: 5-6' ~ '6',
    uncurated$description.14 %in% 'Pathological Gleason Score: 7' ~ '7',
    uncurated$description.14 %in% 'Pathological Gleason Score: 8-10' ~ '8',
    uncurated$description.14 %in% 'Pathological Gleason Score: 5-6' ~ '6',
    uncurated$description.15 %in% 'Clinical Gleason Score: 7' ~ '7',
    uncurated$description.15 %in% 'Clinical Gleason Score: 8-10' ~ '8',
    uncurated$description.15 %in% 'Clinical Gleason Score: 5-6' ~ '6',
    uncurated$description.15 %in% 'Pathological Gleason Score: 7' ~ '7',
    uncurated$description.15 %in% 'Pathological Gleason Score: 8-10' ~ '8',
    uncurated$description.15 %in% 'Pathological Gleason Score: 5-6' ~ '6',
    uncurated$description.16 %in% 'Clinical Gleason Score: 7' ~ '7',
    uncurated$description.16 %in% 'Clinical Gleason Score: 8-10' ~ '8',
    uncurated$description.16 %in% 'Clinical Gleason Score: 5-6' ~ '6',
    uncurated$description.16 %in% 'Pathological Gleason Score: 7' ~ '7',
    uncurated$description.16 %in% 'Pathological Gleason Score: 8-10' ~ '8',
    uncurated$description.16 %in% 'Pathological Gleason Score: 5-6' ~ '6',
    uncurated$description.17 %in% 'Clinical Gleason Score: 7' ~ '7',
    uncurated$description.17 %in% 'Clinical Gleason Score: 8-10' ~ '8',
    uncurated$description.17 %in% 'Clinical Gleason Score: 5-6' ~ '6',
    uncurated$description.17 %in% 'Pathological Gleason Score: 7' ~ '7',
    uncurated$description.17 %in% 'Pathological Gleason Score: 8-10' ~ '8',
    uncurated$description.17 %in% 'Pathological Gleason Score: 5-6' ~ '6',
    uncurated$description.18 %in% 'Clinical Gleason Score: 7' ~ '7',
    uncurated$description.18 %in% 'Clinical Gleason Score: 8-10' ~ '8',
    uncurated$description.18 %in% 'Clinical Gleason Score: 5-6' ~ '6',
    uncurated$description.18 %in% 'Pathological Gleason Score: 7' ~ '7',
    uncurated$description.18 %in% 'Pathological Gleason Score: 8-10' ~ '8',
    uncurated$description.18 %in% 'Pathological Gleason Score: 5-6' ~ '6',
    uncurated$description.19 %in% 'Clinical Gleason Score: 7' ~ '7',
    uncurated$description.19 %in% 'Clinical Gleason Score: 8-10' ~ '8',
    uncurated$description.19 %in% 'Clinical Gleason Score: 5-6' ~ '6',
    uncurated$description.19 %in% 'Pathological Gleason Score: 7' ~ '7',
    uncurated$description.19 %in% 'Pathological Gleason Score: 8-10' ~ '8',
    uncurated$description.19 %in% 'Pathological Gleason Score: 5-6' ~ '6',
    uncurated$description.20 %in% 'Clinical Gleason Score: 7' ~ '7',
    uncurated$description.20 %in% 'Clinical Gleason Score: 8-10' ~ '8',
    uncurated$description.20 %in% 'Clinical Gleason Score: 5-6' ~ '6',
    uncurated$description.20 %in% 'Pathological Gleason Score: 7' ~ '7',
    uncurated$description.20 %in% 'Pathological Gleason Score: 8-10' ~ '8',
    uncurated$description.20 %in% 'Pathological Gleason Score: 5-6' ~ '6',
    uncurated$description.21 %in% 'Clinical Gleason Score: 7' ~ '7',
    uncurated$description.21 %in% 'Clinical Gleason Score: 8-10' ~ '8',
    uncurated$description.21 %in% 'Clinical Gleason Score: 5-6' ~ '6',
    uncurated$description.21 %in% 'Pathological Gleason Score: 7' ~ '7',
    uncurated$description.21 %in% 'Pathological Gleason Score: 8-10' ~ '8',
    uncurated$description.21 %in% 'Pathological Gleason Score: 5-6' ~ '6',
    uncurated$description.22 %in% 'Clinical Gleason Score: 7' ~ '7',
    uncurated$description.22 %in% 'Clinical Gleason Score: 8-10' ~ '8',
    uncurated$description.22 %in% 'Clinical Gleason Score: 5-6' ~ '6',
    uncurated$description.22 %in% 'Pathological Gleason Score: 7' ~ '7',
    uncurated$description.22 %in% 'Pathological Gleason Score: 8-10' ~ '8',
    uncurated$description.22 %in% 'Pathological Gleason Score: 5-6' ~ '6',
    uncurated$description.23 %in% 'Clinical Gleason Score: 7' ~ '7',
    uncurated$description.23 %in% 'Clinical Gleason Score: 8-10' ~ '8',
    uncurated$description.23 %in% 'Clinical Gleason Score: 5-6' ~ '6',
    uncurated$description.23 %in% 'Pathological Gleason Score: 7' ~ '7',
    uncurated$description.23 %in% 'Pathological Gleason Score: 8-10' ~ '8',
    uncurated$description.23 %in% 'Pathological Gleason Score: 5-6' ~ '6')) %>%
    
  dplyr::mutate(T_clinical = dplyr::case_when(
    uncurated$description.9 %in% c('Clinical T: 2a','Clinical T: 2b', 'Clinical T: 2c') ~ 2,
    uncurated$description.9 %in% c('Clinical T: 3a','Clinical T: 3b', 'Clinical T: 3c') ~ 3,
    uncurated$description.9 %in% c('Clinical T: 4', 'Clinical T: 4a') ~ 4,
    uncurated$description.10 %in% c('Clinical T: 2a','Clinical T: 2b', 'Clinical T: 2c') ~ 2,
    uncurated$description.10 %in% c('Clinical T: 3a','Clinical T: 3b', 'Clinical T: 3c') ~ 3,
    uncurated$description.10 %in% c('Clinical T: 4', 'Clinical T: 4a') ~ 4,
    uncurated$description.11 %in% c('Clinical T: 2a','Clinical T: 2b', 'Clinical T: 2c') ~ 2,
    uncurated$description.11 %in% c('Clinical T: 3a','Clinical T: 3b', 'Clinical T: 3c') ~ 3,
    uncurated$description.11 %in% c('Clinical T: 4', 'Clinical T: 4a') ~ 4,
    uncurated$description.12 %in% c('Clinical T: 2a','Clinical T: 2b', 'Clinical T: 2c') ~ 2,
    uncurated$description.12 %in% c('Clinical T: 3a','Clinical T: 3b', 'Clinical T: 3c') ~ 3,
    uncurated$description.12 %in% c('Clinical T: 4', 'Clinical T: 4a') ~ 4,
    uncurated$description.13 %in% c('Clinical T: 2a','Clinical T: 2b', 'Clinical T: 2c') ~ 2,
    uncurated$description.13 %in% c('Clinical T: 3a','Clinical T: 3b', 'Clinical T: 3c') ~ 3,
    uncurated$description.13 %in% c('Clinical T: 4', 'Clinical T: 4a') ~ 4,
    uncurated$description.14 %in% c('Clinical T: 2a','Clinical T: 2b', 'Clinical T: 2c') ~ 2,
    uncurated$description.14 %in% c('Clinical T: 3a','Clinical T: 3b', 'Clinical T: 3c') ~ 3,
    uncurated$description.14 %in% c('Clinical T: 4', 'Clinical T: 4a') ~ 4,
    uncurated$description.15 %in% c('Clinical T: 2a','Clinical T: 2b', 'Clinical T: 2c') ~ 2,
    uncurated$description.15 %in% c('Clinical T: 3a','Clinical T: 3b', 'Clinical T: 3c') ~ 3,
    uncurated$description.15 %in% c('Clinical T: 4', 'Clinical T: 4a') ~ 4,
    uncurated$description.16 %in% c('Clinical T: 2a','Clinical T: 2b', 'Clinical T: 2c') ~ 2,
    uncurated$description.16 %in% c('Clinical T: 3a','Clinical T: 3b', 'Clinical T: 3c') ~ 3,
    uncurated$description.16 %in% c('Clinical T: 4', 'Clinical T: 4a') ~ 4)) %>%
  
  dplyr::mutate(T_substage_clinical = dplyr::case_when(
    uncurated$description.9 %in% c('Clinical T: 2a','Clinical T: 3a', 'Clinical T: 4a') ~ 'a',
    uncurated$description.10 %in% c('Clinical T: 2a','Clinical T: 3a', 'Clinical T: 4a') ~ 'a',
    uncurated$description.11 %in% c('Clinical T: 2a','Clinical T: 3a', 'Clinical T: 4a') ~ 'a',
    uncurated$description.12 %in% c('Clinical T: 2a','Clinical T: 3a', 'Clinical T: 4a') ~ 'a',
    uncurated$description.13 %in% c('Clinical T: 2a','Clinical T: 3a', 'Clinical T: 4a') ~ 'a',
    uncurated$description.14 %in% c('Clinical T: 2a','Clinical T: 3a', 'Clinical T: 4a') ~ 'a',
    uncurated$description.15 %in% c('Clinical T: 2a','Clinical T: 3a', 'Clinical T: 4a') ~ 'a',
    uncurated$description.16 %in% c('Clinical T: 2a','Clinical T: 3a', 'Clinical T: 4a') ~ 'a',
    uncurated$description.9 %in% c('Clinical T: 2b','Clinical T: 3b', 'Clinical T: 4b') ~ 'b',
    uncurated$description.10 %in% c('Clinical T: 2b','Clinical T: 3b', 'Clinical T: 4b') ~ 'b',
    uncurated$description.11 %in% c('Clinical T: 2b','Clinical T: 3b', 'Clinical T: 4b') ~ 'b',
    uncurated$description.12 %in% c('Clinical T: 2b','Clinical T: 3b', 'Clinical T: 4b') ~ 'b',
    uncurated$description.13 %in% c('Clinical T: 2b','Clinical T: 3b', 'Clinical T: 4b') ~ 'b',
    uncurated$description.14 %in% c('Clinical T: 2b','Clinical T: 3b', 'Clinical T: 4b') ~ 'b',
    uncurated$description.15 %in% c('Clinical T: 2b','Clinical T: 3b', 'Clinical T: 4b') ~ 'b',
    uncurated$description.16 %in% c('Clinical T: 2b','Clinical T: 3b', 'Clinical T: 4b') ~ 'b',
    uncurated$description.9 %in% c('Clinical T: 2c','Clinical T: 3c', 'Clinical T: 4c') ~ 'c',
    uncurated$description.10 %in% c('Clinical T: 2c','Clinical T: 3c', 'Clinical T: 4c') ~ 'c',
    uncurated$description.11 %in% c('Clinical T: 2c','Clinical T: 3c', 'Clinical T: 4c') ~ 'c',
    uncurated$description.12 %in% c('Clinical T: 2c','Clinical T: 3c', 'Clinical T: 4c') ~ 'c',
    uncurated$description.13 %in% c('Clinical T: 2c','Clinical T: 3c', 'Clinical T: 4c') ~ 'c',
    uncurated$description.14 %in% c('Clinical T: 2c','Clinical T: 3c', 'Clinical T: 4c') ~ 'c',
    uncurated$description.15 %in% c('Clinical T: 2c','Clinical T: 3c', 'Clinical T: 4c') ~ 'c',
    uncurated$description.16 %in% c('Clinical T: 2c','Clinical T: 3c', 'Clinical T: 4c') ~ 'c')) %>%

   # Pathological N -> Lymph node spread
  dplyr::mutate(N_stage = uncurated$PathologicalN) %>%
   # Pathological M -> Metastatic disease
  dplyr::mutate(M_stage = uncurated$PathologicalM) %>%    
   # Primary site was 'Prostate' in all, no metastases
  dplyr::mutate(sample_type = "primary")   
   # Histology  

curated <- dplyr::select(curated, -(description.6))

clinical_igc <- curated

save(clinical_igc, file = "data-raw/clinical_igc.RData")   
    
 # "Type of Tobacco Use: Cigarettes" %in% uncurated$description.6 ~ 1,
  #"Type of Tobacco Use: Cigars" %in% uncurated$description.6 ~ 1),  
#  is.na(uncurated$description.6) ~ 0))

#####################################################################################
#cBioportal BACA
#####################################################################################

mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
uncurated <- cgdsr::getClinicalData(mycgds, caseList="prad_broad_2013_all")
#mycancerstudy = cgdsr::getCancerStudies(mycgds)
mycaselistren = cgdsr::getCaseLists(mycgds,"prad_broad_2013")

# create the curated object
curated <- initial_curated_df(
  df_rownames = rownames(uncurated),
  template_name="data-raw/template_prad.csv")

curated <- curated %>% 
  dplyr::mutate(study_name = "Baca et al.") %>%
  dplyr::mutate(patient_id = rownames(uncurated)) %>%
  dplyr::mutate(age_at_initial_diagnosis = uncurated$AGE) %>%
  ## TDL: Should follow allowed values assigned in 'template_prad.csv'
  #dplyr::mutate(sample_type=uncurated$SAMPLE_TYPE)%>%
  dplyr::mutate(psa=uncurated$SERUM_PSA)%>%
  dplyr::mutate(gleason_major = as.integer(stringr::str_sub(uncurated$GLEASON_SCORE,1,1))) %>%
  dplyr::mutate(gleason_minor = as.integer(stringr::str_sub(uncurated$GLEASON_SCORE,3,3))) %>%
  dplyr::mutate(T_pathological = readr::parse_number(uncurated$PATH_T_STAGE)) %>%
  dplyr::mutate(T_substage_pathological = dplyr::case_when(
    uncurated$PATH_T_STAGE %in% c("Metastatic") ~ "NA",
    TRUE ~ stringr::str_extract(uncurated$PATH_T_STAGE, "[a-c]+")
  )) %>%
  dplyr::mutate(gleason_grade = gleason_major+gleason_minor) %>%
  #dplyr::mutate(gleason_grade = gsub(";.*","",uncurated$GLEASON_SCORE)) %>%
  dplyr::mutate(grade_group = dplyr::case_when(
    stringr::str_sub(uncurated$GLEASON_SCORE,1,3) == "3+3" ~ "<=6",
    stringr::str_sub(uncurated$GLEASON_SCORE,1,3) == "3+4" ~ "3+4",
    stringr::str_sub(uncurated$GLEASON_SCORE,1,3) == "4+3" ~ "4+3",
    stringr::str_sub(uncurated$GLEASON_SCORE,1,3) %in% c("4+4", "4+5") ~ ">=8"
  )) %>%
  dplyr::mutate(sample_type = dplyr::case_when(
    uncurated$SAMPLE_TYPE == "Metastasis" ~ "metastatic",
    uncurated$SAMPLE_TYPE == "Primary" ~ "primary"
  )) %>%
  # The samples were analyzed in a paired manner to normal tissue
  dplyr::mutate(sample_paired = 1)
  
clinical_baca <- curated

save(clinical_baca, file = "data-raw/clinical_baca.RData")   
  

    
  
