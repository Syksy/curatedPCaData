#' Download gene expression from GEO using study specific id
generate_gex_geo <- function(
  file_directory, 
  geo_code = "GSE25136", # code for Sun et al. (Taylor et al. - GSE21032)
  cleanup = TRUE, 
  collapseFUN = function(z) {apply(z, MARGIN = 2, FUN = stats::median)}, # Function to collapse probe(s) or select a probe, e.g. mean, median, or function that picks a probe with high variance
  ...
){
  if(!missing(file_directory)) here::set_here(file_directory)
  # Supplementary files include the raw CEL files
  supfiles <- GEOquery::getGEOSuppFiles(geo_code)

  # Open the tarball(s)
  utils::untar(tarfile = rownames(supfiles))
  
  ##
  # Sun et al.
  ##
  if(geo_code == "GSE25136"){
    # Make sure to function in a working directory where the are no other tarballs present
    gz_files <- list.files()
    gz_files <- gz_files[grep(".gz", gz_files)]
  
    # Read Affymetrix MA
    Sun <- affy::ReadAffy()
    colnames(affy::exprs(Sun)) <- gsub(".gz|.CEL", "", colnames(Sun))
    
    # Careful not to mask 'rma' from 'affy' by the 'rma' from 'oligo'
    gex <- affy::rma(Sun)
    
    # Removing .CEL and packaging names from the GEO-compatible sample names
    colnames(gex) <- gsub(".CEL.gz", "", colnames(affy::exprs(gex)))
  
    keys <- AnnotationDbi::mappedkeys(hgu133a.db::hgu133aGENENAME)
    nam <- names(as.character(hgu133a.db::hgu133aALIAS2PROBE)[match(rownames(gex),
                                                                    as.character(hgu133a.db::hgu133aALIAS2PROBE))])
    nam[is.na(nam)] <- "NA"
    gex <- do.call("rbind", by(as.matrix(affy::exprs(gex)), INDICES=nam, FUN=collapseFUN))
    
  }
  ##
  # Taylor et al.
  ##
  else if (geo_code == "GSE21032") {
    
    # breaks here -----
    # Error: vector memory exhausted (limit reached?)
    # cels <- oligo::read.celfiles(affy::list.celfiles(), pkgname='pd.huex.1.0.st.v2')	
    
  }
  ##
  # Unknown, throw an R error
  ##
  else{
  	stop("Unknown GEO id, see allowed parameter values for geo_code")
  }

  # Remove downloaded files
  if(cleanup){
    # First GEO download
    file.remove(rownames(supfiles))
    # Tarballs
    file.remove(gz_files)
    # Remove empty folder
    file.remove(paste0(here::here(), "/", geo_code))
  }
  # Return numeric matrix
  as.matrix(gex)
}

#' Download generic 'omics data from cBioPortal using dataset specific query
generate_cbioportal <- function(
  genes, # List of gene symbols to iterate over
  geneticProfiles, # for cgdsr calls, platform and dataset specific string
  caseList, # for cgdsr calls, platform and dataset specific string
  delay = 0.05, # For Sys.sleep to not fetch too fast from cBio API
  splitsize = 100, # How many genes are fetched at one time - max 1000
  verb = TRUE,
  ...
){
  mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
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
    cgdsr::getProfileData(mycgds, genes=splitgenes[[z]], geneticProfiles=geneticProfiles, caseList=caseList)
  })))
}