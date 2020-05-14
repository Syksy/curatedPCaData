rm(list=ls())
source("functions.R")

#get the first file with more sample specfic information
uncurated <- read.delim("../DFCI_version/uncurated_copy/prad_tcga/data_bcr_clinical_data_sample.txt",as.is=TRUE,row.names=2,sep="\t", skip=4)
colnames(uncurated) <- as.character(read.delim("../DFCI_version/uncurated_copy/prad_tcga/data_bcr_clinical_data_sample.txt",as.is=TRUE,header=T,row.names=2,sep="\t",skip=3)[1,])

# the last row in the sample file is an unannotated metastatic case
# just remove it since there is no information for it anyways
uncurated <- uncurated[-dim(uncurated)[1],]

# create the curated object
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_prad.csv")

##study name
curated$study_name <- rep("TCGA, provisional", dim(curated)[1])

## sample name
curated$sample_name <- row.names(uncurated)

## FFPE information
tmp <- uncurated$IS_FFPE
tmp <- gsub("NO", "NA", tmp)
tmp <- gsub("[Not Available]", "NA", fixed = T, tmp)
tmp <- gsub("YES", "ffpe", tmp)
curated$frozen_ffpe <- tmp

## sample type
curated$sample_type <- uncurated$SAMPLE_TYPE

## patient id
curated$patient_id <- uncurated$PATIENT_ID

## alternate patient id
curated$alt_sample_name <- uncurated$OTHER_SAMPLE_ID

###### patient file
uncurated <- read.delim("../DFCI_version/uncurated_copy/prad_tcga/data_bcr_clinical_data_patient.txt",as.is=TRUE,row.names=2,sep="\t", skip=4)
colnames(uncurated) <- as.character(read.delim("../DFCI_version/uncurated_copy/prad_tcga/data_bcr_clinical_data_patient.txt",as.is=TRUE,header=T,row.names=2,sep="\t",skip=3)[1,])
uncurated <- uncurated[curated$patient_id,]

#gleason grade
curated$gleason_grade <- uncurated$GLEASON_SCORE

#gleason major
curated$gleason_major <- uncurated$GLEASON_PATTERN_PRIMARY

#gleason minor
curated$gleason_minor <- uncurated$GLEASON_PATTERN_SECONDARY

#gleason source
curated$source_of_gleason <- rep("biopsy", dim(curated)[1])

#gleason group
tmp <- uncurated$GLEASON_SCORE
tmp <- gsub("6", "<=6", tmp)
tmp <- gsub("8", ">=8", tmp)
tmp <- gsub("9", ">=8", tmp)
tmp <- gsub("10", ">=8", tmp)
tmp[curated$gleason_major == 3 & curated$gleason_minor == 4] <- "3+4"
tmp[curated$gleason_major == 4 & curated$gleason_minor == 3] <- "4+3"
curated$grade_group <- tmp

##primary site of the sample
tmp <- uncurated$PRIMARY_SITE
tmp <- gsub("Peripheral Zone", "peripheral", tmp)
tmp <- gsub("Overlapping / Multiple Zones", "mixed", tmp)
tmp <- gsub("Central Zone", "central", tmp)
tmp <- gsub("Transition Zone", "transitional", tmp)
tmp <- gsub("[Not Available]", "NA", fixed = T, tmp)
curated$zone_of_origin <- tmp

##year of diagnosis
curated$year_diagnosis <- uncurated$INITIAL_PATHOLOGIC_DX_YEAR

##overall survival status
tmp <- uncurated$OS_STATUS
curated$overall_survival_status <- ifelse(grepl("DECEASED",tmp),1,0)
curated$overall_survival_status[is.na(tmp)] <- NA

##overall survival days
curated$days_to_overall_survival <- as.numeric(uncurated$OS_MONTHS) * 30.5

##disease free survival status
tmp <- uncurated$DFS_STATUS
tmp <- gsub("[Not Available]", "NA", fixed = T, tmp)
tmp <- gsub("Recurred/Progressed", "1",  tmp)
tmp <- gsub("DiseaseFree", "0",  tmp)
curated$disease_specific_recurrence_status <- as.numeric(tmp)

##disease free survival days
curated$days_to_disease_specific_recurrence <- as.numeric(uncurated$DFS_MONTHS) * 30.5

##PSA levels
tmp <- uncurated$PSA_MOST_RECENT_RESULTS
tmp <- gsub("[Not Available]", "NA", fixed = T, tmp)
curated$psa <- as.numeric(tmp)

#age at diagnosis
curated$age_at_initial_diagnosis <- uncurated$AGE 

#clinical M stage
tmp <- uncurated$CLIN_M_STAGE
tmp <- gsub("M", "",tmp)
tmp <- gsub("[abc]","",tmp,perl=TRUE)
curated$M_stage <- as.numeric(tmp)

#clinical M substage
tmp <- uncurated$CLIN_M_STAGE
tmp <- gsub("[Not Available]", "NA", fixed = T, tmp)
tmp <- gsub("M", "",tmp)
tmp <- gsub("[^abc]","",tmp,perl=TRUE)
tmp[tmp == ''] <- "NA"
curated$M_substage <- tmp

#clinical T stage
tmp <- uncurated$CLIN_T_STAGE
tmp <- gsub("T", "",tmp)
tmp <- gsub("[abc]","",tmp,perl=TRUE)
curated$T_clinical <- as.numeric(tmp)

#clinical T substage
tmp <- uncurated$CLIN_T_STAGE
tmp <- gsub("[Not Available]", "NA", fixed = T, tmp)
tmp <- gsub("T", "",tmp)
tmp <- gsub("[^abc]","",tmp,perl=TRUE)
tmp[tmp == ''] <- "NA"
curated$T_substage_clinical <- tmp

#pathological T stage
tmp <- uncurated$PATH_T_STAGE
tmp <- gsub("[pT]", "",tmp)
tmp <- gsub("[abc]","",tmp,perl=TRUE)
curated$T_pathological <- as.numeric(tmp)

#pathological T substage
tmp <- uncurated$PATH_T_STAGE
tmp <- gsub("[Not Available]", "NA", fixed = T, tmp)
tmp <- gsub("T", "",tmp)
tmp <- gsub("[^abc]","",tmp,perl=TRUE)
tmp[tmp == ''] <- "NA"
curated$T_substage_pathological <- tmp

#Race
tmp <- uncurated$RACE
tmp <- gsub("[Not Available]", "NA", fixed = T, tmp)
tmp <- gsub("WHITE", "caucasian", tmp)
tmp <- gsub("ASIAN", "asian", tmp)
tmp <- gsub("BLACK OR AFRICAN AMERICAN", "african_american", tmp)
curated$race <- tmp

##tumor purity
curated$tumor_purity_pathology <- uncurated$TUMOR_CONTENT

write.table(curated, row.names = FALSE, file="../DFCI_version/curated/prad_tcga_curated_pdata.txt",sep="\t")

