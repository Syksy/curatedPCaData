rm(list=ls())
source("functions.R")

#get the first file with more sample specfic information
uncurated <- read.delim("../DFCI_version/uncurated_copy/prad_mskcc/data_clinical_sample.txt",as.is=TRUE,row.names=1,sep="\t", skip=4)
colnames(uncurated) <- as.character(read.delim("../DFCI_version/uncurated_copy/prad_mskcc/data_clinical_sample.txt",as.is=TRUE,header=T,row.names=1,sep="\t",skip=3)[1,])

# there are sample ids at the end of the file that dont have any other information, so just remove
uncurated <- uncurated[uncurated$SAMPLE_CLASS != '',]

# create the curated object
curated <- initialCuratedDF(rownames(uncurated),template.filename="template_prad.csv")

##study name
curated$study_name <- rep("Taylor,et al", dim(curated)[1])

## sample name
curated$sample_name <- row.names(uncurated)
curated$patient_id <- row.names(uncurated)

## sample type
tmp <- uncurated$SAMPLE_CLASS
tmp <- gsub("PRIMARY", "primary", tmp)
tmp <- gsub("Metastatic", "metastatic", tmp)
tmp <- gsub("CELL LINE", "cell.line", tmp)
tmp <- gsub("XENOGRAFT", "xenograft", tmp)
curated$sample_type <- tmp

#gleason grade
curated$gleason_grade <- as.numeric(uncurated$GLEASON_SCORE)

#gleason major
curated$gleason_major <- as.numeric(uncurated$GLEASON_SCORE_1)

#gleason minor
curated$gleason_minor <- as.numeric(uncurated$GLEASON_SCORE_2)

#gleason group
tmp <- uncurated$GLEASON_SCORE
tmp <- gsub("6", "<=6", tmp)
tmp <- gsub("8", ">=8", tmp)
tmp <- gsub("9", ">=8", tmp)
tmp <- gsub("10", ">=8", tmp)
tmp[curated$gleason_major == 3 & curated$gleason_minor == 4] <- "3+4"
tmp[curated$gleason_major == 4 & curated$gleason_minor == 3] <- "4+3"
curated$grade_group <- tmp

#ERG fusion GEX
tmp <- uncurated$ERG_FUSION_GEX
tmp <- gsub("positive", "1",  tmp)
tmp <- gsub("not_assessed", "NA",  tmp)
tmp <- gsub("negative", "0",  tmp)
curated$ERG_fusion_GEX <- as.numeric(tmp)

#ERG fusion CNA
tmp <- uncurated$ERG_FUSION_ACGH
tmp <- gsub("positive", "1",  tmp)
tmp <- gsub("not_assessed", "NA",  tmp)
tmp <- gsub("negative", "0",  tmp)
tmp <- gsub("flat", "0",  tmp)
curated$ERG_fusion_CNA <- as.numeric(tmp)

###### patient file
uncurated <- read.delim("../DFCI_version/uncurated_copy/prad_mskcc/data_clinical_patient.txt",as.is=TRUE,row.names=1,sep="\t", skip=4)
colnames(uncurated) <- as.character(read.delim("../DFCI_version/uncurated_copy/prad_mskcc/data_clinical_patient.txt",as.is=TRUE,header=T,row.names=1,sep="\t",skip=3)[1,])
uncurated <- uncurated[curated$patient_id,]

##disease free survival status
tmp <- uncurated$DFS_STATUS
tmp <- gsub("Recurred", "1",  tmp)
tmp <- gsub("DiseaseFree", "0",  tmp)
curated$disease_specific_recurrence_status <- as.numeric(tmp)

##disease free survival days
curated$days_to_disease_specific_recurrence <- as.numeric(uncurated$DFS_MONTHS) * 30.5

#clinical T stage
tmp <- uncurated$CLIN_T_STAGE
tmp <- gsub("T", "",tmp)
tmp <- gsub("[ABC]","",tmp,perl=TRUE)
curated$T_clinical <- as.numeric(tmp)

#clinical T substage
tmp <- uncurated$CLIN_T_STAGE
tmp <- gsub("T", "",tmp)
tmp <- gsub("[^ABC]","",tmp,perl=TRUE)
tmp[tmp == ''] <- "NA"
tmp <- gsub("A", "a",tmp)
tmp <- gsub("B", "b",tmp)
tmp <- gsub("C", "c",tmp)
tmp <- gsub("Na", "NA", tmp)
curated$T_substage_clinical <- tmp

#pathological T stage
tmp <- uncurated$PATH_T_STAGE
tmp <- gsub("T", "",tmp)
tmp <- gsub("[ABC]","",tmp,perl=TRUE)
curated$T_pathological <- as.numeric(tmp)

#pathological T substage
tmp <- uncurated$PATH_T_STAGE
tmp <- gsub("T", "",tmp)
tmp <- gsub("[^ABC]","",tmp,perl=TRUE)
tmp[tmp == ''] <- "NA"
tmp <- gsub("A", "a",tmp)
tmp <- gsub("B", "b",tmp)
tmp <- gsub("C", "c",tmp)
tmp <- gsub("Na", "NA",tmp)
curated$T_substage_pathological <- tmp

write.table(curated, row.names = FALSE, file="../DFCI_version/curated/prad_mskcc_curated_pdata.txt",sep="\t")

