getProfileDataWrapper <- function(
  x, 
  genes, 
  delay = 0.05, # For Sys.sleep to not fetch too fast from cBio API
  splitsize = 100, # How many genes are fetched at one time - max 1000
  verb = TRUE, # If call should be verbose; 0 = silent, 1 = info
  ...
){
  genesplit <- rep(1:ceiling(length(genes$hgnc)/splitsize), each=splitsize)[1:length(genes$hgnc)]
  splitgenes <- split(genes$hgnc, f=genesplit)
  # Fetch split gene name lists as separate calls
  pb <- progress::progress_bar$new(total = length(splitgenes))
  as.matrix(do.call("cbind", lapply(1:length(splitgenes), FUN=function(z){
    # if(verb == TRUE) print(paste("Downloading gene set", z, "of", length(splitgenes), "from cBioportal"))
    if(verb == TRUE) pb$tick()
    # Sleep if necessary to avoid API call overflow
    Sys.sleep(delay)
    # Fetch each split gene name list from the URL-based API, essentially a wrapper for cgdsr's own function
    cgdsr::getProfileData(x=x, genes=splitgenes[[z]], geneticProfiles=geneticProfiles, caseList=caseList)
  })))	
}