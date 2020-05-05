generate_gex_geo <- function(
  file_directory, 
  geo_code = "GSE25136", # code for Sun et al. 
  cleanup = TRUE, 
  collapseFUN = function(z) {apply(z, MARGIN = 2, FUN = median)}, # Function to collapse probe(s) or select a probe, e.g. mean, median, or function that picks a probe with high variance
  ...
){
  if(!missing(file_directory)) here::set_here(file_directory)
  # Supplementary files include the raw CEL files
  supfiles <- GEOquery::getGEOSuppFiles(geo_code)

  # Open the tarball(s)
  utils::untar(tarfile = rownames(supfiles))
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

  # Remove downloaded files
  if(cleanup){
    # First GEO download
    file.remove(rownames(supfiles))
    # Tarballs
    file.remove(gz_files)
    # Remove empty folder
    file.remove(paste0(here::here(), "/", geo_code))
  }
  # TODO: Transform into a MultiAssayExperiment-object prior to returning object (MAE_Sun)
  # Return numeric matrix
  as.matrix(gex)
}

generate_cbioportal <- function(
  genes, # List of gene symbols to iterate over
  geneticProfiles, # for cgdsr calls, platform and dataset specific string
  caseList, # for cgdsr calls, platform and dataset specific string
  ...
){

  mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
  # Use the wrapper function to iterative calls for split gene lists
  dat <- getProfileDataWrapper(
    x = mycgds, # cgdsr object
    genes = tcga_gene_names, # All unique gene symbols
    geneticProfiles = "prad_tcga_pub_rna_seq_v2_mrna", # Omics profile
    caseList = "prad_tcga_pub_sequenced" # Case list
  )
  # TODO: Filter out low quality samples based on list provided by Travis
  # Return omics matrix
  dat
}

