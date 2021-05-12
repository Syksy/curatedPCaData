#' Download gene expression from GEO using study specific id and process it
#'
#' @param geo_code character string indicating name of GEO dataset
#' @param file_directory character string indicating path for downloading raw 
#' GEO data
#' @param cleanup logical value to remove intermediate files 
#' @param collapse_fun function to collapse probe(s) or select a probe, 
#' e.g. mean, median, or function that picks a probe with high variance
#' @param ... additional arguments
generate_gex_geo <- function(
  geo_code = c("GSE21032", # Taylor et al. TODO: Alternative more specific accession code "GSE21034" for GEX
               "GSE25136", # Sun et al.
               "GSE8218",  #Wang et al.
               "GSE6919",  # Chandran et al., Yu et al. from three platforms combined
               "GSE18655", # Barwick et al.
               "GSE2109",  #IGC
               "GSE119616"
               ), 
  file_directory, 
  cleanup = TRUE, 
  collapse_fun = function(z) {apply(z, MARGIN = 2, FUN = stats::median)},
  ...
){
  if(!missing(file_directory)) here::set_here(file_directory)
  # Supplementary files include the raw CEL files
  supfiles <- GEOquery::getGEOSuppFiles(geo_code)

  # Sun et al. -----  
  if(geo_code == "GSE25136"){
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

	# Extracting .CEL and packaging names from the GEO-compatible sample names
	colnames(gex) <- gsub(".CEL.gz", "", colnames(affy::exprs(gex)))

	# Find gene annotations
	keys <- AnnotationDbi::mappedkeys(hgu133a.db::hgu133aGENENAME)
	nam <- names(as.character(hgu133a.db::hgu133aALIAS2PROBE)[match(rownames(gex),
								    as.character(hgu133a.db::hgu133aALIAS2PROBE))])
	nam[is.na(nam)] <- "NA"
	# Collapse probes
	gex <- do.call("rbind", by(as.matrix(affy::exprs(gex)), INDICES = nam, 
	                           FUN = collapse_fun))
	
    # Sun et al does not use Hugo names so the below code makes the names for sun 
	  # et al comparable with those in tcga/taylor
    compare_names <- data.frame(original = row.names(gex),
                               current = limma::alias2SymbolTable(row.names(gex),
                                                                  species="Hs"))
    duplicated_hugo_symbols <- compare_names[duplicated(compare_names$current),]$current

    compare_names <- compare_names %>%
     dplyr::mutate(new_names = dplyr::case_when(
       current %in% duplicated_hugo_symbols ~ original,
       TRUE ~ current
     ))

    row.names(gex) <- compare_names$new_names
    
	# Sort genes to alphabetic order for consistency
	gex <- gex[order(rownames(gex)),]
  }

  #Wang et al.  
  else if(geo_code == "GSE8218"){
    # Make sure to function in a working directory where the are no other tarballs present
  gz_files <- list.files()
  gz_files <- gz_files[grep(".gz", gz_files)]
  
    # Read Affymetrix MA
    wang <- affy::ReadAffy()
    colnames(affy::exprs(wang)) <- gsub(".gz|.CEL", "", colnames(wang))
    
    # Careful not to mask 'rma' from 'affy' by the 'rma' from 'oligo'
    gex <- affy::rma(wang)
    
    # Extracting .CEL and packaging names from the GEO-compatible sample names
    colnames(gex) <- gsub(".CEL.gz", "", colnames(affy::exprs(gex)))
    
    # Find gene annotations
    keys <- AnnotationDbi::mappedkeys(hgu133a.db::hgu133aGENENAME)
    nam <- names(as.character(hgu133a.db::hgu133aALIAS2PROBE)[match(rownames(gex),
                                                                    as.character(hgu133a.db::hgu133aALIAS2PROBE))])
    nam[is.na(nam)] <- "NA"
    # Collapse probes
    gex <- do.call("rbind", by(as.matrix(affy::exprs(gex)), INDICES = nam, 
                               FUN = collapse_fun))
    
    # Sun et al does not use Hugo names so the below code makes the names for sun 
    # et al comparable with those in tcga/taylor
    compare_names <- data.frame(original = row.names(gex),
                                current = limma::alias2SymbolTable(row.names(gex),
                                                                   species="Hs"))
    duplicated_hugo_symbols <- compare_names[duplicated(compare_names$current),]$current
    
    compare_names <- compare_names %>%
      dplyr::mutate(new_names = dplyr::case_when(
        current %in% duplicated_hugo_symbols ~ original,
        TRUE ~ current
      ))
    
    row.names(gex) <- compare_names$new_names
    
    # Sort genes to alphabetic order for consistency
    gex <- gex[order(rownames(gex)),]
  }
  
  # Kim et al
  else if(geo_code == "GSE119616"){
    # Open the tarball(s)
    supfiles <- supfiles[-2,]
    utils::untar(tarfile = rownames(supfiles))
    
    # Use editcelfile package to convert decipher annotations
    gz_files <- list.files()
    gz_files <- gz_files[grep(".gz", gz_files)]
    for (x in gz_files){GEOquery::gunzip(x)}
    celfiles <- list.files()
    celfiles <- celfiles[grep(".CEL", celfiles)]

    # Include custom function that was originally small part of editcelfile (private package available upon request only), to avoid unnecessary non-publicly available dependencies
    # Original code by Jonathan Lehrer & Seagle Liu (R package presented without license, modified code credited here)
    celfileHeaderToHuex <- function(celFilePath, header="HuEx-1_0-st-v2", verbose=TRUE){
	  # NEVER change this.
	  writeBytePos=430
	  overwriteBytes=28 # 2 x 14 characters in "HuEx-1_0-st-v2"

	  if (verbose) {
	    print("ORIGINAL:")
	    print(affyio::read.celfile.header(celFilePath)[1])
	  }

	  f <- base::file(celFilePath, "r+b")

	  base::seek(f, writeBytePos, "start", "write")
	  replicate(overwriteBytes, base::writeBin(as.raw(0), f, useBytes = T)) #clear the old name (14 characters)

	  base::seek(f, writeBytePos, "start", "write")
	  encoded <- stringi::stri_enc_toutf32(header)[[1]]
	  base::writeBin(encoded, f, size=2)

	  base::close(f)

	  cfh <- affyio::read.celfile.header(celFilePath)

	  if (verbose) {
	    print("REVISED:")
	    print(cfh[1])
	  }

	  return(invisible(cfh$cdfName ==  header))
    }

    sapply(celfiles, celfileHeaderToHuex, verbose=TRUE)
    
    # Read in the CEL files 
    CELs <- oligo::read.celfiles(celfiles)
    
    # Perform RMA normalization
    RMAs <- oligo::rma(CELs)
    
    # Obtain gene and sample information
    Biobase::featureData(RMAs) <- oligo::getNetAffx(RMAs, "transcript")
    # GSM######-type names from GEO
    nam0 <- unlist(lapply(strsplit(affy::list.celfiles(), "_"), 
                          FUN = function(z) z[[1]])) 
    # Two naming conventions for the files; picking the PCA###-style 
    nam1 <- unlist(lapply(strsplit(affy::list.celfiles(), "_"), 
                          FUN = function(z) z[[3]])) 
    nam2 <- gsub(".CEL.gz", "", unlist(lapply(strsplit(affy::list.celfiles(), "_"),
                                              FUN = function(z) z[[4]])))
    # Some samples were suffixed with HuEx, while others had Exonl prefix
    nam <- paste(nam0, "_", ifelse(nam1 == "Exon1", nam2, nam1), sep="")
    
    # Extract gene names
    genenames <- unlist(lapply(Biobase::fData(RMAs)[,"geneassignment"], 
                               FUN = function(z) { strsplit(z, " // ")[[1]][2] }))
    
    # Transform into a matrix and remove empty gene names
    gex <- as.matrix(Biobase::exprs(RMAs))
    gex <- gex[-which(is.na(genenames)),]
    rownames(gex) <- genenames[-which(is.na(genenames))]
    # Map the coventionally used Taylor sample names instead of GEO codes 
    # Compatible with e.g. cBioPortal sample names
    # Give unique names with GSM#####.{PCA,PAN}##### combination for uniqueness
    ##colnames(gex) <- nam
    # Use the unique GSM###-names
    colnames(gex) <- unlist(lapply(nam, FUN=function(z){ strsplit(z, "_")[[1]][1] }))
    
    # Sort genes to alphabetic order for consistency
    gex <- gex[order(rownames(gex)),]}
  
  # Taylor  et al.-----
  else if(geo_code == "GSE21032"){ # TODO: Alternative more specific accession code "GSE21034"
	  # Open the tarball(s)
    	utils::untar(tarfile = rownames(supfiles))
	# Read in the CEL files - note: requires a substantial amount of RAM for all 370 samples
	 CELs <- oligo::read.celfiles(affy::list.celfiles())
	
	# Perform RMA normalization
	  RMAs <- oligo::rma(CELs)
	
	# Obtain gene and sample information
	  Biobase::featureData(RMAs) <- oligo::getNetAffx(RMAs, "transcript")
	# GSM######-type names from GEO
	  nam0 <- unlist(lapply(strsplit(affy::list.celfiles(), "_"), 
	                        FUN = function(z) z[[1]])) 
	# Two naming conventions for the files; picking the PCA###-style 
	  nam1 <- unlist(lapply(strsplit(affy::list.celfiles(), "_"), 
	                      FUN = function(z) z[[3]])) 
	  nam2 <- gsub(".CEL.gz", "", unlist(lapply(strsplit(affy::list.celfiles(), "_"),
	                                          FUN = function(z) z[[4]])))
	# Some samples were suffixed with HuEx, while others had Exonl prefix
	  nam <- paste(nam0, "_", ifelse(nam1 == "Exon1", nam2, nam1), sep="")

	# Extract gene names
	  genenames <- unlist(lapply(Biobase::fData(RMAs)[,"geneassignment"], 
	                           FUN = function(z) { strsplit(z, " // ")[[1]][2] }))

	# Transform into a matrix and remove empty gene names
  	gex <- as.matrix(Biobase::exprs(RMAs))
  	gex <- gex[-which(is.na(genenames)),]
  	rownames(gex) <- genenames[-which(is.na(genenames))]
  	# Map the coventionally used Taylor sample names instead of GEO codes 
  	# Compatible with e.g. cBioPortal sample names
  	# Give unique names with GSM#####.{PCA,PAN}##### combination for uniqueness
  	##colnames(gex) <- nam
  	# Use the unique GSM###-names
  	colnames(gex) <- unlist(lapply(nam, FUN=function(z){ strsplit(z, "_")[[1]][1] }))
  	
  	# Sort genes to alphabetic order for consistency
  	gex <- gex[order(rownames(gex)),]
  }
	

  # Chandran et al.-----  
  else if(geo_code == "GSE6919"){
	  # Open the tarball(s)
	  utils::untar(tarfile = rownames(supfiles))

	# Three different platforms were used; need to read them separately with ReadAffy

	gse <- GEOquery::getGEO("GSE6919", GSEMatrix = TRUE)
	# ...
	# GSM152839 - GSM152855 : U95C
	# GSM152856 - GSM152880 : U95Av2
	# GSM152881 - GSM152905 : U95B
	# GSM152906 - GSM152930 : U95C
	# GSM152931 - GSM152991 : U95Av2
	# GSM152992 - GSM153053 : U95B
	# GSM153054 - GSM153114	: U95C
	# ...
	
	# Read Affymetrix MA

	# U95Av2
	GPL_Av2 <- rownames(Biobase::pData(gse[[1]]))
	Chandran_Av2 <- affy::ReadAffy(filenames=paste0(GPL_Av2, ".CEL.gz"))

	# In U95B non-valid files, i.e.
	# Error: the following are not valid files:
    	#   GSM152822.CEL.gz

	# U95B
	broken <- c("GSM152822")
	GPL_U95B <- rownames(Biobase::pData(gse[[2]]))
	GPL_U95B <- GPL_U95B[-which(GPL_U95B %in% broken)]
	Chandran_U95B <- affy::ReadAffy(filenames=paste0(GPL_U95B, ".CEL.gz"))

	# U95C
	Chandran_U95C <- affy::ReadAffy(filenames=paste0(GPL_U95C, ".CEL.gz"))
	GPL_U95C <- rownames(Biobase::pData(gse[[3]]))

	# Careful not to mask 'rma' from 'affy' by the 'rma' from 'oligo'
	gex_Av2 <- affy::rma(Chandran_Av2)
	gex_U95B <- affy::rma(Chandran_U95B)
	gex_U95C <- affy::rma(Chandran_U95C)
	
	# Annotation dbs mapping manufacturer ids to gene symbols
	map_Av2 <- as.list(hgu95av2.db::hgu95av2SYMBOL[AnnotationDbi::mappedkeys(hgu95av2.db::hgu95av2SYMBOL)])
	map_U95B <- as.list(hgu95b.db::hgu95bSYMBOL[AnnotationDbi::mappedkeys(hgu95b.db::hgu95bSYMBOL)])
	map_U95C <- as.list(hgu95c.db::hgu95cSYMBOL[AnnotationDbi::mappedkeys(hgu95c.db::hgu95cSYMBOL)])
	# Extract mapped gene symbols per row (not all have a gene annotation)
  	mapped_Av2 <- c(map_Av2[match(rownames(gex_Av2), names(map_Av2))])
  	mapped_Av2 <- unlist(lapply(mapped_Av2, FUN=function(x) ifelse(is.null(x[1]), NA, x[1])))
  	mapped_U95B <- c(map_U95B[match(rownames(gex_U95B), names(map_U95B))])
  	mapped_U95B <- unlist(lapply(mapped_U95B, FUN=function(x) ifelse(is.null(x[1]), NA, x[1])))
  	mapped_U95C <- c(map_U95C[match(rownames(gex_U95C), names(map_U95C))])
  	mapped_U95C <- unlist(lapply(mapped_U95C, FUN=function(x) ifelse(is.null(x[1]), NA, x[1])))

	# Collapse per gene symbols over probes	
	gex_Av2_mapped <- do.call("rbind", (by(as.matrix(gex_Av2), INDICES=mapped_Av2, FUN=function(x){
		apply(x, MARGIN=2, FUN=collapse_fun)
	})))
	gex_U95B_mapped <- do.call("rbind", (by(as.matrix(gex_U95B), INDICES=mapped_U95B, FUN=function(x){
		apply(x, MARGIN=2, FUN=collapse_fun)
	})))
	gex_U95C_mapped <- do.call("rbind", (by(as.matrix(gex_U95C), INDICES=mapped_U95C, FUN=function(x){
		apply(x, MARGIN=2, FUN=collapse_fun)
	})))
	
	# Merge according to common gene symbols
	tmp <- merge(gex_Av2_mapped, gex_U95B_mapped, by="row.names", all.x=TRUE, all.y=TRUE)
	rownames(tmp) <- tmp[,1]
	tmp <- tmp[,-1]
	tmp <- merge(tmp, gex_U95C_mapped, by="row.names", all.x=TRUE, all.y=TRUE)
	rownames(tmp) <- tmp[,1]
	tmp <- tmp[,-1]
	colnames(tmp) <- gsub(".CEL.gz", "", colnames(tmp))
	tmp <- tmp[order(rownames(tmp)),]

	gex <- tmp
  }

  # GSE2109 - IGC
  else if(geo_code == "GSE2109"){
    # Open the tarball(s)
    utils::untar(tarfile = rownames(supfiles[17,]))
    clinical <- rio::import("data-raw/clinical_igc.RData")
    # Make sure to function in a working directory where the are no other tarballs present
    rownames(clinical) <- paste0(rownames(clinical),'.CEL.gz')
    
    gz_files <- rownames(clinical)
    gz_files <- gz_files[grep(".gz", gz_files)]
    gz_files <- as.data.frame(gz_files)
    
    # Read Affymetrix MA
    igc <- affy::ReadAffy(filenames=gz_files$gz_files)
    colnames(affy::exprs(igc)) <- gsub(".gz|.CEL", "", colnames(igc))
    
    # Careful not to mask 'rma' from 'affy' by the 'rma' from 'oligo'
    gex <- affy::rma(igc)
    
    # Extracting .CEL and packaging names from the GEO-compatible sample names
    colnames(gex) <- gsub(".CEL.gz", "", colnames(affy::exprs(gex)))
    
    # Find gene annotations
    keys <- AnnotationDbi::mappedkeys(hgu133a.db::hgu133aGENENAME)
    nam <- names(as.character(hgu133a.db::hgu133aALIAS2PROBE)[match(rownames(gex),
                                                                    as.character(hgu133a.db::hgu133aALIAS2PROBE))])
    nam[is.na(nam)] <- "NA"
    # Collapse probes
    gex <- do.call("rbind", by(as.matrix(affy::exprs(gex)), INDICES = nam, 
                               FUN = collapse_fun))
    
    # Sun et al does not use Hugo names so the below code makes the names for sun 
    # et al comparable with those in tcga/taylor
    compare_names <- data.frame(original = row.names(gex),
                                current = limma::alias2SymbolTable(row.names(gex),
                                                                   species="Hs"))
    duplicated_hugo_symbols <- compare_names[duplicated(compare_names$current),]$current
    
    compare_names <- compare_names %>%
      dplyr::mutate(new_names = dplyr::case_when(
        current %in% duplicated_hugo_symbols ~ original,
        TRUE ~ current
      ))
    
    row.names(gex) <- compare_names$new_names
    
    # Sort genes to alphabetic order for consistency
    gex <- gex[order(rownames(gex)),]
  }
  # Barwick et al.-----  
  
  else if(geo_code == "GSE18655"){
	# Custom DASL
	# .gz
	GEOquery::gunzip(rownames(supfiles), overwrite=TRUE)
	tmp <- read.table("./GSE18655/GSE18655_HCP_Toronto_raw.txt", header=TRUE, row.names=1, skip=4)
	# Contains
	#                 X1_rep1  X1_rep2      X10 X100_rep1 X100_rep2
	#GI_10092618-S-3 30329.56 28241.03 28005.69  30165.83  28913.85
	#GI_10092618-S-1 25222.37 22141.40 25871.17  26419.28  23814.04
	#GI_10092618-S-2 30202.68 30705.60 32448.47  30405.46  29000.44
	
	# Download GPL annotations for the custom platform
	gpl <- GEOquery::getGEO("GSE18655", GSEMatrix = FALSE, getGPL = TRUE)
	map_barwick <- gpl@gpls[[1]]@dataTable@table
	# Contains
	#	     ID SequenceSource      GB_ACC
	#1 GI_10092618-S         RefSeq NM_020529.1
	#2 GI_10337586-S         RefSeq NM_020996.1
	#3 GI_10834981-S         RefSeq NM_000599.1
	
	# Drop the third '-' split suffix from tmp rownames
	map_tmp <- unlist(lapply(rownames(tmp), FUN=function(x) { paste(strsplit(x, "-")[[1]][1:2], collapse="-") } ))
	# Collapse over multiple hits to same RefSeq within a sample
	gex <- do.call("rbind", by(tmp, INDICES=map_barwick[match(map_tmp, map_barwick$ID),"GB_ACC"], FUN=collapse_fun))
	# Omit RefSeq versions from mapping to gene symbols
	rownames(gex) <- unlist(lapply(rownames(gex), FUN=function(x) { strsplit(x, ".", fixed=TRUE)[[1]][1] }))
	
	# Hugo symbols
	symbols <- curatedPCaData:::curatedPCaData_genes[match(rownames(gex), curatedPCaData:::curatedPCaData_genes[,"refseq_mrna"]),"hgnc_symbol"]
	gex <- gex[!is.na(symbols),]
	symbols <- symbols[!is.na(symbols)]
	rownames(gex) <- symbols
	gex <- gex[order(rownames(gex)),]	
  }
  
  # Unknown GEO id (throw an error) -----
  else{
  	stop("Unknown GEO id, see allowed parameter values for geo_code")
  }

  ## TODO: Move cleanup inside each GEO as file types and custom files are very cohort specific
  # Remove downloaded files
  if(cleanup){
    # First GEO download
    file.remove(rownames(supfiles))
    # Tarballs
    #file.remove(gz_files)
    # Remove empty folder
    file.remove(paste0(here::here(), "/", geo_code))
  }

  gex <- as.matrix(gex)
  gex <- gex %>% janitor::remove_empty(which = c("rows", "cols"))
  gex
}

#' Download copy number variant data from GEO using study specific id and process it
#' 
#' @param geo_code character string indicating name of GEO dataset. Default is "GSE21035"
#' (Taylor et al)
#' @param file_directory character string indicating path for downloading raw 
#' GEO data
#' @param cleanup logical value to remove intermediate files 
#' @param ... additional arguments
generate_cna_geo <- function(
  geo_code = c("GSE21035", # Taylor et al.
               "GSE54691" # Hieronymus et al.
               ),
  file_directory, 
  cleanup = TRUE, 
  ...
){
  if(!missing(file_directory)) here::set_here(file_directory)
  # Supplementary files include the raw CEL files
  #supfiles <- GEOquery::getGEOSuppFiles(geo_code, filter_regex='tar')
  supfiles <- GEOquery::getGEOSuppFiles(geo_code)
  

  # Hieronymus et al. requires going to a subfolder
  # you should NEVER use setwd() in a package - causes massive downstream problems 
  # need to remove 
  # if(geo_code == "GSE54691"){
  #   setwd("GSE54691") # Contains both CNA input as well as clinical data txt.gz
  # }
  
  # Open the tarball(s)
  #rownames(supfiles) <- shQuote(rownames(supfiles))
  utils::untar(tarfile = rownames(supfiles))
  
  ##
  # same rCGH pipeline applied to datasets:
  # Taylor et al.
  # Hieronymus et al.
  ##
  if(geo_code %in% c("GSE21035", "GSE54691")){
  	# For now, the package 'rCGH' has to be available in the workspace,
    # otherwise below functions will fail on e.g. rCGH::adjustSignal and when trying to find 'hg18'
  	# Read in Agilent 2-color data
  	cna <- lapply(list.files(geo_code), FUN = function(z) { 
  		try({
  			cat("\n\nProcessing: ",z,"\n\n") 
  			# Taylor et al. (and & Hieronymus et al.)
  			#if(geo_code %in% c("GSE21035", "GSE54691")){
			rCGH::readAgilent(z, genome = "hg38", sampleName = gsub(".txt.gz", "", z)) 
  			# Hieronymus et al. - any reason to use hg19 since the locations are known?
  			#}else if(geo_code == "GSE54691"){
  			#	rCGH::readAgilent(z, genome = "hg19", sampleName = gsub(".txt.gz", "", z)) 
  			#}
  		})
  	})
	#> list.files()[which(unlist(lapply(cna, FUN=class))=="try-error")]
	#[1] "GSM525755.txt" "GSM525763.txt"
	# Some files appear broken in Taylor et al; missing columns?
	
	# Not all files can always be successfully processed
	tryerr <- which(lapply(cna, FUN = class) == "try-error")
	if(length(tryerr)>0){
		# Warn of try-errors
		warning(paste("Error while processing files: ", 
			# Collapse file names
			paste(list.files()[which(unlist(lapply(cna, FUN=class))=="try-error")], collapse=", ")
		))
		# Omit data that could not be succcessfully read
		cna <- cna[-tryerr]
	}
	
	# Signal adjustments
	cna <- lapply(cna, FUN = function(z){
		try({
			rCGH::adjustSignal(z) 
		})
	})
	# Segmentation
	cna <- lapply(cna, FUN = function(z){
		try({
			rCGH::segmentCGH(z) 
		})
	})
	# EM-algorithm normalization
	cna <- lapply(cna, FUN = function(z){
		try({
			rCGH::EMnormalize(z) 
		})
	})
	# Remove additional suffixes from sample names
	cna <- lapply(cna, FUN = function(z){ 
		try({
			if(!"rCGH-Agilent" %in% class(z)){
				stop("This function is intended for Agilent aCGH analyzed with rCGH R Package (class \'rCGH-Agilent\')")		
			}
			# e.g. transform "GSM525575.txt|.gz" -> "GSM525575"
			z@info["sampleName"] <- gsub(pattern = ".gz|.txt", replacement = "", z@info["fileName"])
			z
		}) 
	})
	# Save sample names separately (of 'length(cna)')
	samplenames <- unlist(lapply(cna, FUN = function(z) { z@info["sampleName"] }))
	# Reformat samplenames in Hieronymus
	if(geo_code =="GSE54691"){
		# Pick GSM-part in GSM###_PCA###
		samplenames <- unlist(lapply(samplenames, FUN=function(z) { strsplit(z, "_")[[1]][[1]]}))
	}	
	# Get segmentation table
	cna <- lapply(cna, FUN = function(z){
		try({
			rCGH::getSegTable(z)
		})
	})
	# Get per-gene table
	cna <- lapply(cna, FUN = function(z){
		try({
			rCGH::byGeneTable(z)
		})
	})
	# Extract all unique gene symbols present over all samples
	genenames <- unique(unlist(lapply(cna, FUN = function(z) { z$symbol })))
	# Bind genes to rows, name samples afterwards
	cna <- do.call("cbind", lapply(cna, FUN = function(z){
		# Return CNAs as Log2Ratios
		z[match(genenames, z$symbol), "Log2Ratio"]
	}))
	# Name rows and columns to genes and sample names, respectively
	rownames(cna) <- genenames
	colnames(cna) <- samplenames
	# CNA matrix is ready
	cna <- as.matrix(cna)
  }
  # Other (placeholder) - should probably place this in the beginning so it 
  # doesnt try to run through all the beginning steps 
  else if(geo_code == ""){
    stop("Must supply a GEO id")
  }
  # Unknown
  else{
    stop("Unknown GEO id, see allowed parameter values for geo_code")
  }

  # Remove downloaded files
  # TODO: Make sure appropriate permissions exist for removing files
  if(cleanup){
    # First GEO download
    file.remove(rownames(supfiles))
    # TODO: Tarballs
    #file.remove(gz_files)
    # Remove empty folder
    file.remove(paste0(here::here(), "/", geo_code))
  }
  # Return numeric matrix
  cna <- as.matrix(cna)
  cna <- cna %>% janitor::remove_empty(which = c("rows", "cols"))
  cna

}


#' Download generic 'omics data from cBioPortal using dataset specific query
#' 
#' @param genes character vector of genes to query
#' @param geneticProfiles charatcer string of cbioportal genetic profiles
#' @param caseList charcter string of patient IDs for that genetic profile
#' @param delay numberic value for delay time between querying gene sets
#' @param splitsize number of genes in each query
#' @param verb logical value for displaying progress bar 
generate_cbioportal <- function(
  genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)), 
  geneticProfiles = c("prad_tcga_pub_rna_seq_v2_mrna", #TCGA GEX 
                      "prad_tcga_pub_gistic", # TCGA CNA (GISTIC)
                      "prad_tcga_pub_linear_CNA", # TCGA CNA (Capped relative linear copy-number values)
                      "prad_tcga_pub_mutations", # TCGA mutation data
                      # prad_tcga_pub_fusion, # TCGA fusions
                      "prad_mskcc_mrna_median_Zscores", # Taylor et al. GEX (z-score normalized)
                      "prad_mskcc_cna", # Taylor et al. CNA
                      "prad_broad_mrna", # PRAD Broad GEX
                      "prad_broad_cna", # PRAD Broad CNA
                      "prad_eururol_2017_rna_seq_mrna", # PRAD Eururol GEX
                      "prad_eururol_2017_cna", # PRAD Eururol CNA
                      "prad_su2c_2019_mrna_seq_fpkm_capture_all_sample_Zscores", # Abida et al. FPKM values capture assay
                      "prad_su2c_2019_mrna_seq_fpkm_polya_all_sample_Zscores", # Abida et al. FPKM values PolyA
                      "prad_su2c_2019_gistic" #Abida et al. CNA
                      ), # for cgdsr calls, platform and dataset specific string
  caseList = c("prad_tcga_pub_sequenced", # TCGA
               "prad_mskcc_sequenced", # Taylor et al.
               "prad_broad_sequenced", # PRAD Broad
               "prad_eururol_2017_sequenced", # PRAD Eururol
               "prad_su2c_2019_sequenced" #Abida et al.
               ), # for cgdsr calls, platform and dataset specific string
  delay = 0.05, 
  splitsize = 100, 
  verb = TRUE
){
  # If given genes is a list (with slots for various annotation types), try to extract hugo gene symbols
  if(class(genes)=="list"){
    genes <- genes$hgnc_symbol
  }
  # Establisigh connection to cBioPortal
  mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
  # Split gene name vector into suitable lengths
  genesplit <- rep(1:ceiling(length(genes)/splitsize), 
                   each = splitsize)[1:length(genes)]
  splitgenes <- split(genes, f = genesplit)
  # Fetch split gene name lists as separate calls
  pb <- progress::progress_bar$new(total = length(splitgenes))
  # Bind the API calls as per columns
  ret <- as.matrix(do.call("cbind", lapply(1:length(splitgenes), FUN = function(z){
      if(verb == TRUE) pb$tick()
      # Sleep if necessary to avoid API call overflow
      Sys.sleep(delay)
      # Fetch each split gene name list from the URL-based API, essentially a wrapper for cgdsr's own function
      cgdsr::getProfileData(mycgds, genes = splitgenes[[z]], 
                            geneticProfiles = geneticProfiles, caseList = caseList)
    })))  
  
  ret <- t(ret)
  # Remove fully empty rows/columns (redundancy)
  ret <- ret %>% janitor::remove_empty(which = c("rows", "cols"))
  
}

#' Download mutation data from cBioPortal and format them into an oncoprint-friendly matrix
#'
#' @param study_id TODO
#' @param genes TODO
#' @param delay TODO
#' @param splitsize TODO
#' @param verb TODO
#' @param oncoprintify TODO
#' @examples
#' oncop_tcga <- curatedPCaData:::generate_cbioportal_oncoprint(study_id="tcga", oncoprintify=TRUE)
#' oncop_taylor <- curatedPCaData:::generate_cbioportal_oncoprint(study_id="taylor", oncoprintify=TRUE)
#' oncop_barbieri <- curatedPCaData:::generate_cbioportal_oncoprint(study_id="barbieri", oncoprintify=TRUE)
#' oncop_ren <- curatedPCaData:::generate_cbioportal_oncoprint(study_id="ren", oncoprintify=TRUE)
#' # TODO
#' # oncop_abida <- curatedPCaData:::generate_cbioportal_oncoprint(study_id="abida", oncoprintify=TRUE)
generate_cbioportal_oncoprint <- function(
	study_id, # tcga, taylor, barbieri, ren, or abida
	genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)),
	delay = 0.05,
	splitsize = 100,
	verb = TRUE,
	oncoprintify = TRUE
){
	study_id <- tolower(study_id)
	# Remember cBioPortal genetic profile names and query according to study id
	# Mutations
	geneticProfileMutations = c(
		tcga = "prad_tcga_pub_mutations", # TCGA mutations
		taylor = "prad_mskcc_mutations", # Taylor et al. mutations
		barbieri = "prad_broad_mutations", # Barbieri et al. mutations
		ren = "prad_eururol_2017_mutations", # Ren et al. mutations
		abida = "prad_su2c_2019_mutations" # Abida et al. mutations
	)
	# CNA listings for GISTIC
	geneticProfileCNA = c( # Gistic-level mutation info, {-2,-1,0,1,2}
		tcga = "prad_tcga_pub_gistic", # TCGA (GISTIC instead of capped relative linear copy numbers)
		taylor = "prad_mskcc_cna", # Taylor et al., GISTIC
		barbieri = "prad_broad_cna", # Barbieri et al., GISTIC
		ren = "prad_eururol_2017_cna", # Ren et al. (GISTIC instead of capped relative linear copy numbers)
		abida = "prad_su2c_2019_gistic" # Abida et al., GISTIC
	)
	# Sample ID list
	caseList = c(
		tcga = "prad_tcga_pub_sequenced", # TCGA samples
		taylor = "prad_mskcc_sequenced", # Taylor et al. samples
		barbieri = "prad_broad_sequenced", # Barbieri et al. samples
		ren = "prad_eururol_2017_sequenced", # Ren et al. samples
		abida = "prad_su2c_2019_sequenced" # Abida et al. samples 		
	)

	# If given genes is a list (with slots for various annotation types), try to extract hugo gene symbols
	if(class(genes)=="list"){
		genes <- genes$hgnc_symbol
	}
	# Establisigh connection to cBioPortal
	mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
	# Split gene name vector into suitable lengths
	genesplit <- rep(1:ceiling(length(genes)/splitsize), each = splitsize)[1:length(genes)]
	splitgenes <- split(genes, f = genesplit)
	# Fetch split gene name lists as separate calls
	pb <- progress::progress_bar$new(total = length(splitgenes))
	# Bind the API calls as per columns
	if(verb) print("Downloading mutation data...")
	mut <- do.call("rbind", lapply(1:length(splitgenes), FUN = function(z){
		# Sleep if necessary to avoid API call overflow
		Sys.sleep(delay)
		# Fetch each split gene name list from the URL-based API, essentially a wrapper for cgdsr's own function
		cgdsr::getMutationData(mycgds, 
			genes = splitgenes[[z]], # 
			geneticProfile = geneticProfileMutations[study_id], 
			caseList = caseList[study_id]
		)[,c("case_id", "gene_symbol", "mutation_type")]
	}))
	if(verb) print("Downloading copy number alteration (GISTIC) data...")
	cna <- as.matrix(do.call("cbind", lapply(1:length(splitgenes), FUN = function(z){
		# Sleep if necessary to avoid API call overflow
		Sys.sleep(delay)
		# Fetch each split gene name list from the URL-based API, essentially a wrapper for cgdsr's own function
		cgdsr::getProfileData(mycgds, genes = splitgenes[[z]], 
				geneticProfile = geneticProfileCNA[study_id], 
				caseList = caseList[study_id]
			)
		}
	)))  
	cna <- t(cna)	
	cna[which(cna=="NaN")] <- NA

	if(!oncoprintify){
		# Return raw mutation and cna calls
		list(mut, cna)
	}else{
		# Return an oncoprint-friendly, textified matrix with ;-separators
		# Create base text matrix from GISTIC CNAs
		if(verb) print("Appending CNAs to the oncoprint")
		oncop <- t(apply(cna, MARGIN=1, FUN=function(z){
			c("DEEP_DEL", "SHALLOW_DEL", "", "LOW_GAIN", "HIGH_AMP")[z+3]
		}))
		dimnames(oncop) <- dimnames(cna)
		# Add samples from mutation list that may have been lacking from CNA matrix
		if(!all(mut$case_id %in% colnames(oncop))){
			samps <- matrix(NA, nrow=nrow(oncop), ncol=sum(!unique(mut$case_id) %in% colnames(oncop)))
			colnames(samps) <- unique(mut$case_id)[!unique(mut$case_id) %in% colnames(oncop)]
			oncop <- cbind(oncop, samps)
		}
		if(verb) print("Appending mutations to the oncoprint")
		# Add gene names from mutation list that may have been lacking from CNA matrix
		if(!all(mut$gene_symbol %in% rownames(oncop))){
			genes <- matrix(NA, nrow=sum(!unique(mut$gene_symbol) %in% rownames(oncop)), ncol=ncol(oncop))
			rownames(genes) <- unique(mut$gene_symbol)[!unique(mut$gene_symbol) %in% rownames(oncop)]
			oncop <- rbind(oncop, genes)
		}
		# Sanitize '-' symbols to '.' in mutations for genes and sample names similar to how CNA matrix has them
		mut[,"gene_symbol"] <- gsub("-", ".", mut[,"gene_symbol"])
		mut[,"case_id"] <- gsub("-", ".", mut[,"case_id"])
		# Append mutations
		for(i in 1:nrow(mut)){
			if(is.na(oncop[mut[i,"gene_symbol"], mut[i,"case_id"]])){
				oncop[mut[i,"gene_symbol"], mut[i,"case_id"]] <- mut[i,"mutation_type"]
			}else{
				oncop[mut[i,"gene_symbol"], mut[i,"case_id"]] <- paste(oncop[mut[i,"gene_symbol"], mut[i,"case_id"]], mut[i,"mutation_type"], sep=";")
			}
		}
		#apply(mut, MARGIN=1, FUN=function(z){
		#	if(z["gene_symbol"] %in% rownames(oncop) & z["case_id"] %in% colnames(oncop)){
		#		oncop[z["gene_symbol"],z["case_id"]] <- ifelse(is.na(oncop[z["gene_symbol"],z["case_id"]]), oncop[z["gene_symbol"],z["case_id"]], paste(oncop[z["gene_symbol"],z["case_id"]], z["mutation_type"], sep=";"))
		#	}else{
		#		warning(paste("Gene symbol", z["gene_symbol"], "of case_id", z["case_id"], "not found from oncoprint matrix"))
		#	}
		#})
		if(verb) print("Subsetting samples to match those present in MAE")
		# Sample names from the corresponding MAE object
		cols <- switch(study_id,
			"tcga" = MultiAssayExperiment::colData(curatedPCaData::mae_tcga)$sample_name,
			"taylor" = MultiAssayExperiment::colData(curatedPCaData::mae_taylor)$patient_id,
			"barbieri" = MultiAssayExperiment::colData(curatedPCaData::mae_barbieri)$sample_name,
			"ren" = MultiAssayExperiment::colData(curatedPCaData::mae_ren)$sample_name,
			"abida" = MultiAssayExperiment::colData(curatedPCaData::mae_abida)$sample_name
		)
		#oncop <- oncop %>% janitor::remove_empty(which = c("rows", "cols"))
		# NA columns if samples are not present
		oncop <- oncop[,match(cols, colnames(oncop))]
		colnames(oncop) <- cols
		oncop
	}
}

#' Download and generate omics from the ICGC
#'
#' @param icgc_id character string indicating name of ICGC dataset
#' @param file_directory character string indicating path for downloading the ICGC files
#' @param omic character string indicating which omic to download and generate for the specified icgc_id
#'
#' ICGC Publication Policy for embargoes etc: http://www.icgc.org/icgc/goals-structure-policies-guidelines/e3-publication-policy
#' ICGC Publication guidelines: http://docs.icgc.org/portal/publication/#current-moratorium-status-for-icgc-projects . 
#' (e.g. 'All data shall become free of a publication moratorium when either the data is published by the ICGC member project or 
#'  one year after a specified quantity of data (e.g. genome dataset from 100 tumours per project) has been released via the ICGC 
#'  database or other public databases. In all cases data shall be free of a publication moratorium two years after its initial release.')
#' "PRAD-CA: No Embargo. Data available without limitations"
#' "PRAD-UK: No Embargo. Data available without limitations"
#' "PRAD-FR: No Embargo. Data available without limitations"
#' "PRAD-CN: Publication limitations in place until 2020-04-30" (expired) (Only simple mutations, excluded)
#' "EOPC-DE: No Embargo. Data available without limitations" (Need to verify biological applicability) 
#'
#' NOTE: Sometimes the downloads seem to fail randomly; perhaps a fixed amount of retries ought to be allowed?
generate_icgc <- function(
	icgc_id = "PRAD_CA", # Study which ought to be downloaded; Canadian Prostate Adenocarcima study as default; note ICGC uses format 'PRAD-CA' but '_' is used for R-friendliness
	set = "gex", # Which dataset (patient or sample data / omics platform) to try to extract from the data; valid values: 'clinical', 'gex', 'cna', ...
	file_directory, # Temporary download location
	collapseFUN = mean, # Function to use to average/median etc over multiple genes/probes mapped to same hugo symbol
	verb = 0 # Level of verbosity; 0 = minimal, 1 = informative, 2 = debugging
){
	# Use upper case similar to ICGC's naming conventions, also map '-' into '_' as it is not a mathematical operator and is safer
	icgc_id <- base::gsub('-', '_', base::toupper(icgc_id))
	# Currently the studies from Canada, UK and France have enough samples & omics to fit to the package
	if(!icgc_id %in% c("PRAD_CA", "PRAD_FR", "PRAD_UK")){
		stop("The queried ICGC study ought to be one of: 'PRAD_CA', 'PRAD_UK', or 'PRAD_FR'")		
	}
	# If provided, set to a custom temporary download directory
	if(!missing(file_directory)) here::set_here(file_directory)
	# At the time of writing, the latest public release is 28
	# The downloadable files depend on study and fixed URLs are listed here-in
	icgc_sets <- list(
		# Canada (Release 28 fixed reference, instead of 'current', for reproducibility)
		"PRAD_CA" = c(
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-CA/copy_number_somatic_mutation.PRAD-CA.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-CA/donor.PRAD-CA.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-CA/donor_exposure.PRAD-CA.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-CA/donor_family.PRAD-CA.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-CA/donor_therapy.PRAD-CA.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-CA/exp_array.PRAD-CA.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-CA/meth_array.PRAD-CA.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-CA/sample.PRAD-CA.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-CA/simple_somatic_mutation.open.PRAD-CA.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-CA/specimen.PRAD-CA.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-CA/structural_somatic_mutation.PRAD-CA.tsv.gz"			
		),
		# France (Release 28 fixed reference, instead of 'current', for reproducibility)
		"PRAD_FR" = c(
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-FR/copy_number_somatic_mutation.PRAD-FR.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-FR/donor.PRAD-FR.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-FR/donor_family.PRAD-FR.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-FR/donor_surgery.PRAD-FR.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-FR/exp_array.PRAD-FR.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-FR/exp_seq.PRAD-FR.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-FR/sample.PRAD-FR.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-FR/simple_somatic_mutation.open.PRAD-FR.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-FR/specimen.PRAD-FR.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-FR/structural_somatic_mutation.PRAD-FR.tsv.gz"
		),
		# United Kingdom (Release 28 fixed, instead of 'current', for reproducibility)
		"PRAD_UK" = c(
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-UK/copy_number_somatic_mutation.PRAD-UK.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-UK/donor.PRAD-UK.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-UK/donor_exposure.PRAD-UK.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-UK/donor_family.PRAD-UK.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-UK/donor_therapy.PRAD-UK.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-UK/meth_array.PRAD-UK.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-UK/sample.PRAD-UK.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-UK/simple_somatic_mutation.open.PRAD-UK.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-UK/specimen.PRAD-UK.tsv.gz",
			"https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/PRAD-UK/structural_somatic_mutation.PRAD-UK.tsv.gz"
		)
		# Possibly? (Early Onset Prostate Cancer, Germany)
		# https://dcc.icgc.org/releases/release_28/Projects/EOPC-DE
	)
	# Small internal function to assist with the downloads
	.icgcDownload <- function(url){
		# Pick the filename from the end of the URL
		filename <- strsplit(url, "/")
		filename <- filename[[1]][[length(filename[[1]])]]
		# Download file into parsed *.tsv.gz 
		utils::download.file(url=url, destfile=filename)
		# gunzip the files open
		GEOquery::gunzip(filename, overwrite=TRUE)
		filename
	}
	# Loop over the list of urls, download and gunzip them
	files <- icgc_sets[[icgc_id]]
	# Select subset of files to download
	if(set == "clinical"){ # Clinical data
		files <- grep("sample|donor|specimen", files, value=TRUE)
	}else if(set == "gex"){ # Gene expression (array or sequencing)
		files <- grep("exp_array|exp_seq", files, value=TRUE)
	}else if(set == "cna"){ # Copy number alterations
		files <- grep("copy_number_somatic_mutation", files, value=TRUE)
	}else if(set == "mut"){ # Point mutations
		files <- grep("simple_somatic_mutation", files, value=TRUE)
	}else if(set == "met"){ # Methylation
		files <- grep("meth_array", files, value=TRUE)
	}else if(set == "str"){ # Larger structural aberrations
		files <- grep("structural_somatic_mutation", files, value=TRUE)
	}else{
		stop("Invalid parameter 'set'; should be one of: 'clinical, 'gex', 'cna', 'mut', 'met', or 'str'")
	}
	# list apply the .icgcDownload-function to all eligible files
	files <- unlist(lapply(files, FUN=.icgcDownload))
	# Verbosity
	if(verb>=1) print(paste("Files downloaded & extracted:", paste(files, collapse=", ")))
	# Gunzip has extracted the compressed files, removing the suffix accordingly
	files <- gsub(".gz", "", files)
	# Pick an 'omics based on the requested 'omic' and process into a suitable data.frame or matrix
	if(set == "clinical"){ # Clinical data
		# Main id file, should contain only simple file with prefix 'sample.*'
		ret <- read.table(grep("sample", files, value=TRUE), sep="\t", header=TRUE, stringsAsFactors=FALSE)
		# Other files to append:
		# donor.* contains multiple useful clinical features
		# donor_family.* possibly interesting family history related to cancer
		# Append useful donor.* -file contained information to the clinical data
		don <- read.table(grep("donor.", files, value=TRUE, fixed=TRUE), sep="\t", header=TRUE, stringsAsFactors=FALSE)
		ret <- cbind(ret, don[match(ret$icgc_donor_id, don$icgc_donor_id),])
		# Append useful specimen.* -file containing information e.g. for gleason grade
		spe <- read.table(grep("specimen.", files, value=TRUE, fixed=TRUE), sep="\t", header=TRUE, stringsAsFactors=FALSE)
		ret <- cbind(ret, spe[match(ret$icgc_donor_id, spe$icgc_donor_id),])
	}else if(set == "gex"){ # Gene expression (array or sequencing)
		# In RNA-seq, use raw read counts and normalize them
		if(any(grepl("exp_seq", files))){
			ret <- read.table(grep("exp_seq", files, value=TRUE), sep="\t", header=TRUE, stringsAsFactors=FALSE)[,c("icgc_sample_id", "gene_id", "raw_read_count")]
		# In microarrays, need to use the prenormalized expression values
		}else if(any(grepl("exp_array", files))){
			ret <- read.table(grep("exp_array", files, value=TRUE), sep="\t", header=TRUE, stringsAsFactors=FALSE)[,c("icgc_sample_id", "gene_id", "normalized_expression_value")]
		}else{
			stop("Unknown gene expression file type")
		}
		# Mapping of genes to hugo symbols is dependent on the dataset
		if(icgc_id == "PRAD_CA"){
			# Normalized expression values based on refseq gene ids, add hugo symbols
			ret$hgnc_symbol <- curatedPCaData:::curatedPCaData_genes[match(ret$gene_id, curatedPCaData:::curatedPCaData_genes$refseq_mrna),"hgnc_symbol"]
			# Remove entries for which gene_id or hgnc_symbol is missing
			ret <- ret[-which(is.na(ret$gene_id) | is.na(ret$hgnc_symbol)),]
			# Create a gex matrix of samples as columns and genes as rows
			genes <- sort(unique(ret$hgnc_symbol))
			genes <- genes[-which(genes=="")]
			samples <- unique(ret$icgc_sample_id)
			# Loop over samples, inner loop for mapping genes and then collapsing if multiple measurements exist for same hgnc symbol
			ret <- do.call("cbind", by(ret, 
				INDICES=ret$icgc_sample_id, 
				FUN=function(x) { 
					unlist(lapply(genes, FUN=function(z){
						collapseFUN(x[match(z, x$hgnc_symbol),"normalized_expression_value"])
					}))
				}
			))
			rownames(ret) <- genes
		}else if(icgc_id == "PRAD_FR"){
			# Normalizing raw counts
			
		}
	}else if(set == "cna"){ # Copy number alterations	
		ret <- read.table(grep("copy_number_somatic", files, value=TRUE), sep="\t", header=TRUE, stringsAsFactors=FALSE)[,c("icgc_sample_id", "gene_affected", "mutation_type", "copy_number", "chromosome", "chromosome_start", "chromosome_end", "assembly_version")]
	}else if(set == "mut"){ # Point mutations
		ret <- read.table(grep("simple_somatic_mutation", files, value=TRUE), sep="\t", header=TRUE, stringsAsFactors=FALSE)[,c("icgc_sample_id", "gene_affected", "mutation_type", "chromosome", "chromosome_start", "chromosome_end", "assembly_version", "chromosome_strand", "reference_genome_allele", "mutated_from_allele", "mutated_to_allele", "quality_score", "consequence_type", "aa_mutation", "cds_mutation", "transcript_affected")]
	}else if(set == "met"){ # Methylation
		ret <- read.table(grep("meth_array", files,, value=TRUE), sep="\t", header=TRUE, stringsAsFactors=FALSE)[,c("icgc_sample_id")]
	}else if(set == "str"){ # Larger structural aberrations
		ret <- read.table(grep("structural_somatic_mutation", files, value=TRUE), sep="\t", header=TRUE, stringsAsFactors=FALSE)[,c("icgc_sample_id", "variant_type", "sv_id", "placement", "annotation", "interpreted_annotation", "chr_from", "chr_from_strand", "chr_from_range", "chr_from_flanking_seq", "chr_to", "chr_to_bkpt", "chr_to_strand", "chr_to_range", "chr_to_flanking_seq", "assembly_version", "gene_affected_by_bkpt_to", "gene_affected_by_bkpt_to", "transcript_affected_by_bkpt_from", "transcript_affected_by_bkpt_to", "bkpt_from_context", "bkpt_to_context", "gene_build_version", "platform")]
	}
	# Return constructed object (TODO: or raw data to be further processed/debugged)
	ret 
}

		