mapstage <- function(x){
  output <- x
  output[output=="I"] <- 1
  output[output=="II"] <- 2
  output[output=="III"] <- 3
  output[output=="IV"] <- 4
  output <- as.integer(output)
  return(output)
}

initialCuratedDF <- function(DF.rownames,template.filename){
  template <- read.csv(template.filename,as.is=TRUE)
  output <- matrix(NA,
                   ncol=nrow(template),
                   nrow=length(DF.rownames))
  colnames(output) <- template$col.name
  rownames(output) <- DF.rownames
  output <- data.frame(output)
  for (i in 1:ncol(output)){
    class(output[,i]) <- template[i,"var.class"]
  }
  output$sample_name <- DF.rownames
  return(output)
}

curate_grade <- function(characteristic, word="grade: ") {
 tmp <- sub(word, "", characteristic, fixed=TRUE)

 tmp[tmp=="I"] <- "low"
 tmp[tmp=="II"] <- "low"
 tmp[tmp=="III"] <- "high"
 tmp
}

##ASK LEVI ABOUT THIS:
##curate_title <- function(uncurated){
##tmp <- uncurated$title
##curated$alt_sample_name <- tmp
##}

curate_age <- function(characteristic, word="Age: ") {
 tmp <- sub(word, "", characteristic, fixed=TRUE)
}

#note, only works when the row w/ "stage" values contains 1 integer
#and only works, I think, when the row w/ "stage" begins w/ "Tumor stage:", although, I could probably put an OR operator...
curate_stage <- function(characteristic, word ="Tumor stage: ") {
tmp <- sub(word, "", characteristic, fixed=TRUE)
tmp <- gsub("[^\\d]","",tmp,perl=TRUE)
}

#Note, has same limitations as curate_stage
curate_substage <- function(characteristic, word ="Tumor stage: T") {
tmp <- sub(word, "", characteristic, fixed=TRUE)
tmp <- gsub("[^abcd]","",tmp,perl=TRUE)
#tmp[tmp==""] <- NA
}

curate_gleasongrade <- function(characteristic, word ="Gleason Grade:") {
tmp <- sub(word, "", characteristic, fixed=TRUE)
}

##ASK LEVI: WHY DOESN'T THIS WORK?
#curate_summarygrade <- function(curated)  {
#tmp <- curated$gleasongrade
#tmp[tmp=="2"] <- "low"
#tmp[tmp=="3"] <- "low"
#tmp[tmp=="4"] <- "low"
#tmp[tmp=="5"] <- "low"
#tmp[tmp=="6"] <- "low"
#tmp[tmp=="7"] <- "intermediate"
#tmp[tmp=="8"] <- "high"
#tmp[tmp=="9"] <- "high"
#tmp[tmp=="10"] <- "high" 
#}

##ASK LEVI: WHY DOESN'T THIS WORK?
#curate_summarygrade <- function(characteristic, word ="Gleason Grade: ") {
#tmp <- sub(word, "", characteristic, fixed=TRUE)
# tmp[tmp=="2"] <- "low"
# tmp[tmp=="3"] <- "low"
# tmp[tmp=="4"] <- "low"
# tmp[tmp=="5"] <- "low"
# tmp[tmp=="6"] <- "low"
# tmp[tmp=="7"] <- "intermediate"
# tmp[tmp=="8"] <- "high"
# tmp[tmp=="9"] <- "high"
# tmp[tmp=="10"] <- "high" 
#}

#postProcess <- function(curated){
#  tmp <- curated$days_to_death < 6*30 & curated$vital_status=="deceased"
#  tmp[is.na(tmp)] <- FALSE
#  curated$inferred_chemo_response[tmp] <- "refractory"
#  tmp <- curated$days_to_death >12*30 & curated$vital_status=="living"
#  tmp[is.na(tmp)] <- FALSE
#  curated$inferred_chemo_response[tmp] <- "sensitive"
#  return(curated)
#}
