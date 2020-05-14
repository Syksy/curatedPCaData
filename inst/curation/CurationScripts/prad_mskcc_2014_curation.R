rm(list=ls())
source("functions.R")

uncurated <- read.delim("../DFCI_version/uncurated_copy/prad_mskcc_2014/data_clinical_patient.txt",as.is=TRUE,row.names=1,sep="\t", skip=4)
colnames(uncurated) <- as.character(read.delim("../DFCI_version/uncurated_copy/prad_mskcc_2014/data_clinical_patient.txt",as.is=TRUE,header=T,row.names=1,sep="\t",skip=3)[1,])

##initial creation of curated dataframe
##curated <- data.frame(sample_name=rownames(uncurated))#,row.names=rownames(uncurated))
##initial creation of curated dataframe
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_prad.csv")

##--------------------
##start the mappings
##--------------------

####################################
# some treatment information might be included


##study name
curated$study_name <- rep("Hieronymus, et al.", dim(curated)[1])

##patient id
curated$patient_id <- row.names(uncurated)

#age at diagnosis
curated$age_at_initial_diagnosis <- uncurated$AGE 

##overall survival status
tmp <- uncurated$OS_STATUS
curated$overall_survival_status <- ifelse(grepl("DECEASED",tmp,ignore.case=TRUE),1,0)
curated$overall_survival_status[is.na(tmp)] <- NA

##overall survival days
tmp <- uncurated$OS_MONTHS
tmp <- tmp * 30.5
curated$days_to_overall_survival <- tmp 

##metastasis survival status
tmp <- uncurated$METASTATIC_DISEASE
curated$metastasis_occurrence_status <- ifelse(grepl("POSITIVE",tmp,ignore.case=TRUE),1,0)
curated$metastasis_occurrence_status[is.na(tmp)] <- NA

##metastasis survival days
tmp <- uncurated$METS_FREE_TIME_MONTHS
tmp <- tmp * 30.5
curated$days_to_metastatic_occurrence <- tmp 

##disease survival status
tmp <- uncurated$DFS_STATUS
curated$disease_specific_recurrence_status <- ifelse(grepl("Recurred",tmp,ignore.case=TRUE),1,0)
curated$disease_specific_recurrence_status[is.na(tmp)] <- NA

##disease survival days
tmp <- uncurated$DFS_MONTHS
tmp <- tmp * 30.5
curated$days_to_disease_specific_recurrence <- tmp 

##PSA
curated$psa <- as.numeric(uncurated$PSA)

#clinical T stage
tmp <- uncurated$CLIN_T_STAGE
tmp <- gsub("T", "",tmp)
tmp <- gsub("[abc]","",tmp,perl=TRUE)
curated$T_clinical <- tmp

#clinical T substage
tmp <- uncurated$CLIN_T_STAGE
tmp <- gsub("T", "",tmp)
tmp <- gsub("[^abc]","",tmp,perl=TRUE)
curated$T_substage_clinical <- tmp

#pathological T stage
tmp <- uncurated$PATH_T_STAGE
tmp <- gsub("[pT]", "",tmp)
tmp <- gsub("[abc]","",tmp,perl=TRUE)
curated$T_pathological <- tmp

#pathological T substage
tmp <- uncurated$PATH_T_STAGE
tmp <- gsub("[pT]", "",tmp)
tmp <- gsub("[^abc]","",tmp,perl=TRUE)
curated$T_substage_pathological <- tmp


#get the second file with more sample specfic information
uncurated <- read.delim("../DFCI_version/uncurated_copy/prad_mskcc_2014/data_clinical_sample.txt",as.is=TRUE,row.names=1,sep="\t", skip=4)
colnames(uncurated) <- as.character(read.delim("../DFCI_version/uncurated_copy/prad_mskcc_2014/data_clinical_sample.txt",as.is=TRUE,header=T,row.names=1,sep="\t",skip=3)[1,])
uncurated <- uncurated[row.names(curated),]

#tumor sample site
tmp <- uncurated$SAMPLE_TYPE
tmp <- gsub("Primary", "primary",tmp)
curated$sample_type <- tmp

#gleason grade
curated$gleason_grade <- uncurated$BIOPSY_GLEASON_SCORE

#gleason major
curated$gleason_major <- uncurated$BIOPSY_GLEASON_SCORE_1

#gleason minor
curated$gleason_minor <- uncurated$BIOPSY_GLEASON_SCORE_2

#gleason source
curated$source_of_gleason <- rep("biopsy", dim(curated)[1])

#gleason group
tmp <- uncurated$BIOPSY_GLEASON_SCORE
tmp <- gsub("6", "<=6", tmp)
tmp <- gsub("8", ">=8", tmp)
tmp <- gsub("9", ">=8", tmp)
tmp <- gsub("10", ">=8", tmp)
tmp[curated$gleason_major == 3 & curated$gleason_minor == 4] <- "3+4"
tmp[curated$gleason_major == 4 & curated$gleason_minor == 3] <- "4+3"
curated$grade_group <- tmp


write.table(curated, row.names = FALSE, file="../DFCI_version/curated/prad_mskcc_2014_curated_pdata.txt",sep="\t")
