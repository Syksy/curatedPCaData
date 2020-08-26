#' Download gene expression from GEO using study specific id
generate_gex_geo <- function(
  file_directory, 
  ## Allowed GEO ids:
  # "GSE25136" : Sun et al.
  # "GSE21032" : Taylor et al.
  geo_code = "GSE25136", # code for Sun et al. (Taylor et al. - GSE21032)
  cleanup = TRUE, 
  collapse_probes = function(z) {apply(z, MARGIN = 2, FUN = stats::median)}, # Function to collapse probe(s) or select a probe, e.g. mean, median, or function that picks a probe with high variance
  # Function for cleaning rows/cols where cBio returned NaN or similar non-finite values only
  clean_columns = janitor::remove_empty,
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
  
    # keys <- AnnotationDbi::mappedkeys(hgu133a.db::hgu133aGENENAME)
    nam <- names(as.character(hgu133a.db::hgu133aALIAS2PROBE)[match(rownames(gex),
                                                                    as.character(hgu133a.db::hgu133aALIAS2PROBE))])
    nam[is.na(nam)] <- "NA"
    gex <- do.call("rbind", by(as.matrix(affy::exprs(gex)), INDICES=nam, FUN=collapse_probes))
    
    compare_names <- data.frame(original = row.names(gex),
                                current = limma::alias2SymbolTable(row.names(gex),
                                                                   species="Hs"))
    duplicated_hugo_symbols <- compare_names[duplicated(compare_names$current),]$current
    
    compare_names <- compare_names %>% 
      mutate(new_names = case_when(
        current %in% duplicated_hugo_symbols ~ original,
        TRUE ~ current
      ))
    
    row.names(gex) <- compare_names$new_names
    
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
  as.matrix(clean_columns(gex))
  
  # clean names to match conventions? 
  # gex <- gex %>% janitor::clean_names()
  
}

#' Download generic 'omics data from cBioPortal using dataset specific query
generate_cbioportal <- function(
  # By default used the gene symbol from package data
  genes = sort(curatedPCaData:::tcga_gene_names$hgnc), # List of gene symbols to iterate over
  # geneticProfiles; Allowed platform/data combinations:
  # "prad_tcga_pub_rna_seq_v2_mrna" : TCGA GEX
  # "prad_tcga_pub_gistic" : TCGA CNA (GISTIC)
  # "prad_mskcc_mrna_median_Zscores" : Taylor et al. GEX (z-score normalized)
  # "prad_mskcc_cna" : Taylor et al. CNA
  geneticProfiles = "prad_tcga_pub_rna_seq_v2_mrna", # for cgdsr calls, platform and dataset specific string
  # caseList; Allowed platform/data combinations:
  # "prad_tcga_pub_sequenced" : TCGA
  # "prad_mskcc_sequenced" : Taylor et al. 
  caseList = "prad_tcga_pub_sequenced", # for cgdsr calls, platform and dataset specific string
  # For Sys.sleep to not fetch too fast from cBio API
  delay = 0.05, 
  # Amount of genes fetched at a single API call - max 1000
  splitsize = 100, 
  # Function for cleaning rows/cols where cBio returned NaN or similar non-finite values only
  clean_columns = janitor::remove_empty,
  # Verbosity
  verb = TRUE,
  ...
){
  mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
  genesplit <- rep(1:ceiling(length(genes)/splitsize), each=splitsize)[1:length(genes)]
  #genesplit <- rep(1:ceiling(length(genes$hgnc)/splitsize), each=splitsize)[1:length(genes$hgnc)]
  splitgenes <- split(genes, f=genesplit)
  #splitgenes <- split(genes$hgnc, f=genesplit)
  # Fetch split gene name lists as separate calls
  pb <- progress::progress_bar$new(total = length(splitgenes))
  # Transpose to have genes as row and samples as columns
  gex <- clean_columns(t(
    # Bind the API calls as per columns
    as.matrix(do.call("cbind", lapply(1:length(splitgenes), FUN=function(z){
      # if(verb == TRUE) print(paste("Downloading gene set", z, "of", length(splitgenes), "from cBioportal"))
      if(verb == TRUE) pb$tick()
      # Sleep if necessary to avoid API call overflow
      Sys.sleep(delay)
      # Fetch each split gene name list from the URL-based API, essentially a wrapper for cgdsr's own function
      cgdsr::getProfileData(mycgds, genes=splitgenes[[z]], geneticProfiles=geneticProfiles, caseList=caseList)
    }))
  )))  
  
  # if cleaning column names 
  # gex <- gex %>% janitor::clean_names()
  
}
