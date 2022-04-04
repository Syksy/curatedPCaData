#' Download gene expression from GEO using study specific id and process it
#'
#' @param geo_code character string indicating name of GEO dataset
#' @param file_directory character string indicating path for downloading raw 
#' GEO data
#' @param cleanup logical value to remove intermediate files 
#' @param collapse_fun function to collapse probe(s) or select a probe, 
#' e.g. mean, median, or function that picks a probe with high variance
#' @param ... additional arguments
#' @noRd
#' @keywords internal
generate_gex_geo <- function(
	# Unique GEO identifier
	geo_code = c(
		"GSE18655",  # Barwick et al.
		"GSE6919",   # Chandran et al., Yu et al. from three platforms combined
		"GSE134051", # Friedrich et al.
		"GSE119616", # Kim et al.
		"GSE14206",  # Kunderfranco et al.
		"GSE2109",   # IGC
		"GSE25136",  # Sun et al.
		"GSE21032",  # Taylor et al. Alternative more specific accession code "GSE21034" for GEX
		"GSE5132",   # True et al.
		"GSE8218",   # Wang et al.
		"GSE6956",   # Wallace et al.
		"GSE157548"  # Weiner et al.
		), 
	# Indicate whether the 'oligo' or the 'affy' package should be the primary means for processing the CEL-data		
	pckg = "oligo", 
	# Working directory
	file_directory, 
	# Clean up (i.e. delete) files post-download/-processing
	cleanup = TRUE, 
	# Regex filter for GEOquery download
	filter_regex,
	# Function for collapsing rows
	collapse_fun = function(z) {apply(z, MARGIN = 2, FUN = stats::median)},
	...
){
	if(!missing(file_directory)) here::set_here(file_directory)
	# Supplementary files include the raw CEL files, possibility to use regular expression filter
	if(missing(filter_regex)){
		supfiles <- GEOquery::getGEOSuppFiles(geo_code)
	}else{
		supfiles <- GEOquery::getGEOSuppFiles(geo_code, filter_regex = filter_regex)
	}

	if(!pckg %in% c("oligo", "affy", "limma", "other")) stop(paste0("Invalid processing method parameter pckg (should be either 'oligo', 'affy', 'limma', or 'other'):", pckg))


  
	# Barwick et al.
	if(geo_code == "GSE18655"){
		# Following lumi-package variance stabilizing transformation normalization
		# https://www.bioconductor.org/packages//2.13/bioc/vignettes/lumi/inst/doc/lumi.pdf
	
		# Custom DASL, just .gz (gunzipped) raw file
		GEOquery::gunzip(rownames(supfiles), overwrite=TRUE)
		tmp <- read.table("./GSE18655/GSE18655_HCP_Toronto_raw.txt", header=TRUE, row.names=1, skip=4)
		# Contains
		#                 X1_rep1  X1_rep2      X10 X100_rep1 X100_rep2
		#GI_10092618-S-3 30329.56 28241.03 28005.69  30165.83  28913.85
		#GI_10092618-S-1 25222.37 22141.40 25871.17  26419.28  23814.04
		#GI_10092618-S-2 30202.68 30705.60 32448.47  30405.46  29000.44

		# No probe level variance information available; log-transforming after quantile normalization		
		gex <- log2(preprocessCore::normalize.quantiles(as.matrix(tmp)))
		dimnames(gex) <- dimnames(tmp)
		tmp <- gex
		
		# Download GPL annotations for the custom platform
		gpl <- GEOquery::getGEO(geo_code, GSEMatrix = FALSE, getGPL = TRUE)
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
		
		# Collapse multiple instances of same gene symbol
		gex <- do.call("rbind", by(as.matrix(gex), INDICES=rownames(gex), FUN=function(x){
			collapse_fun(x)
		}))
	}
  
  	# Chandran et al.
  	else if(geo_code == "GSE6919"){
  		# Open the tarball(s)
  		utils::untar(tarfile = rownames(supfiles))
		gz_files <- list.files()
		gz_files <- gz_files[grep(".gz", gz_files)]	
  	
  		if(pckg == "oligo"){
  			# Download the GSE that contains indicators which samples were run on the most modern chip type Av2
			gse <- GEOquery::getGEO(geo_code, GSEMatrix = TRUE)
			gpl_av2 <- rownames(Biobase::pData(gse[[1]]))
			#gex <- affy::ReadAffy(filenames=paste0(gpl_av2, ".CEL.gz"))
			gex <- oligo::read.celfiles(paste0(gpl_av2, ".CEL.gz"))
			# Normalize background convolution of noise and signal using RMA (median-polish)
			gex <- oligo::rma(gex)
			# Extract expression matrix with probe ids
			gex <- oligo::exprs(gex)
			# Sanitize column names
			colnames(gex) <- gsub(".CEL.gz", "", colnames(gex))
			# > Biobase::featureData(gex) <- oligo::getNetAffx(gex, "transcript")
			# Error in oligo::getNetAffx(gex, "transcript") : 
			#  NetAffx Annotation not available in 'pd.hg.u95av2'. Consider using 'biomaRt'.
			# Using hgu95av2.db directly:
			map <- AnnotationDbi::select(hgu95av2.db::hgu95av2.db, rownames(gex), c("SYMBOL"))
			map <- map[match(rownames(gex),map$PROBEID),]
			# Collapsing gene expression for a given symbol:
			gex <- do.call("rbind", by(gex, INDICES=map$SYMBOL, FUN=collapse_fun))
  		}
  		else if(pckg == "affy"){
  			## DEPRECATED ANNOTATION OVER CHIP TYPES
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
  	}
  	
	# Friedrich et al. 
	else if(geo_code == "GSE134051"){
		# TDL: Re-implement
  		# Open the tarball(s)
  		utils::untar(tarfile = rownames(supfiles))
		gz_files <- list.files()
		gz_files <- gz_files[grep(".gz", gz_files)]	
		for (x in gz_files){GEOquery::gunzip(x)}
		files <- gsub(".gz", "", gz_files)
		# Single color Agilent array 
		if(pckg == "limma"){
			# Download GPL for probe identifier annotations
			gpl <- GEOquery::getGEO(geo_code, GSEMatrix = FALSE, getGPL = TRUE)		
			# Extract annotations in ENSG####
			gpl <- GEOquery::Table(slot(gpl, "gpls")[[1]])
			
			# Process using limma pipeline for Agilent single color custom arrays
			# Read in raw data
			gex <- limma::read.maimages(dir(".", "txt"), "agilent.median", green.only = TRUE)
			# Perform background correction (and apparently log-transformation of intensities)
			# 'method="normexp" a convolution of normal and exponential distributions is fitted to the foreground intensities using the 
			# background intensities as a covariate, and the expected signal given the observed foreground becomes the corrected intensity.'
			gex <- limma::backgroundCorrect(gex, method="normexp", offset=1)
			# ^ normexp.method="rma" is a possibility
			# Normalize between arrays using quantile normalization
			gex <- limma::normalizeBetweenArrays(gex, method="quantile")
			# Extract probe average expressions and use Agilent probe identifiers
			gex <- limma::avereps(gex, ID=gex$genes$ProbeName)
			# Transform into a matrix and truncate column names
			gex <- as.matrix(gex)
			colnames(gex) <- lapply(colnames(gex), FUN=function(x) { strsplit(x, "_")[[1]][1] })
			
			# Map rownames to ENSG#### 
			ensg <- gpl[match(rownames(gex), gpl$ID),"SPOT_ID"]
			ensg[which(ensg == "NoEntry")] <- NA
			genes <- curatedPCaData:::curatedPCaData_genes$hgnc_symbol[match(ensg, curatedPCaData:::curatedPCaData_genes$ensembl_gene_id)]
			genes[which(genes == "")] <- NA
			# Collapse using hugo gene symbols
			gex <- do.call("rbind", by(gex, INDICES=genes, FUN=collapse_fun))
			# Include non-redundant names (non-NAs)
			gex <- gex[!rownames(gex) %in% c(""),]
		}else if(pckg == "oligo"){
			stop("Unavailable")
		}else if(pckg == "affy"){
			stop("Old code not processing raw data")
			# TDL: Original code by FC moved from download-data.R for concordance with other raw data processing and fixed

			# load series and platform data from GEO
			fr_gset <- GEOquery::getGEO("GSE134051", GSEMatrix =TRUE, getGPL=TRUE)

			labels = Biobase::fData(fr_gset[[1]])
			gtab = curatedPCaData:::curatedPCaData_genes

			if (length(fr_gset) > 1) idx <- grep("GPL26898", attr(gset, "names")) else idx <- 1
			fr_ex <- Biobase::exprs(fr_gset[[idx]])

			# replacing row names with gene ids
			##############################################
			labels$ensb = substr(labels$SPOT_ID, 1, 15)
			rownames(fr_ex) = labels$ensb
			fr_ex <- fr_ex[rownames(fr_ex) != 'NoEntry', ]
			fr_ex <- fr_ex[substr(rownames(fr_ex), 1, 4) != 'XLOC', ]
			fr_ex <- fr_ex[is.element(rownames(fr_ex), gtab[,1]), ]
			gtab2 <- gtab[match(rownames(fr_ex), gtab[,1]), ]

			gtab2[which(gtab2[,3] == ''), 3] <- gtab2[which(gtab2[,3] == ''), 1]

			rownames(fr_ex) <- gtab2[,3]
			## Typo ?
			#gex <- aggregate(fr_ex, by = list(rownames(ex)), mean)
			gex <- aggregate(fr_ex, by = list(rownames(fr_ex)), mean)
			rownames(gex) <- gex[,1]
			gex <- gex[, -1]
			gex <- gex[order(rownames(gex)),]
		}
	}
  
	# Kim et al.
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
    
	}

	# Kunderfranco et al.
	# GPL887	Agilent-012097 Human 1A Microarray (V2) G4110B (Feature Number version)
	else if(geo_code == "GSE14206"){
		# Two colour Agilent array
		utils::untar(tarfile = rownames(supfiles))
		gz_files <- list.files()
		gz_files <- gz_files[grep(".gz", gz_files)]	
		for (x in gz_files){GEOquery::gunzip(x)}
		files <- gsub(".gz", "", gz_files)
		# Files inlcude 'GPL887_old_annotation.txt', grep this away
		files <- files[-grep("GPL887", files)]
		if(pckg=="limma"){
			# New code RAW data processed in limma
			# Download GPL for probe identifier annotations
			gpl <- GEOquery::getGEO(geo_code, GSEMatrix = FALSE, getGPL = TRUE)		
			# Extract annotations in ENSG####
			gpl <- GEOquery::Table(slot(gpl, "gpls")[[1]])
			# Read in raw data
			gex <- limma::read.maimages(files, "agilent", green.only = FALSE)
			# Background correction
			gex <- limma::backgroundCorrect(gex, method="normexp", offset=1)
			# Normalize with global loess
			gex <- limma::normalizeWithinArrays(gex, method="loess")
			# Average expression
			gex <- limma::avereps(gex, ID=gex$genes$ProbeName)
			# Extract sample names (GSM#####)
			samples <- unlist(lapply(strsplit(colnames(gex$M)[seq(from=1, to=ncol(gex$M), by=2)], "_"), FUN=function(x) { x[[1]][1] }))
			# For diagnostics (mean should be at M=0): > limma::plotMA(gex, array=1)
			# Checking dye-swapped correlations (flipped log ration)
			#> quantile(unlist(lapply(seq(from=1, to=ncol(gex$M), by=2), FUN=function(x) { cor(gex$M[,x], -gex$M[,x+1], method="spearman") })), probs=seq(0,1,by=.1))
			#        0%        10%        20%        30%        40%        50%        60%        70%        80%        90%       100% 
			#-0.2760865  0.2179309  0.2766869  0.3413689  0.3967651  0.4219280  0.4615115  0.5176102  0.5468144  0.5776107  0.6099629
			# ...
			# Visual diagnostics 
			# 
			
			# For now take the average of the normal and flipped dye swapped log ratios:
			gex$M <- do.call("cbind", lapply(seq(from=1, to=ncol(gex$M), by=2), FUN=function(x) { apply(cbind(-gex$M[,x], gex$M[,x+1]), MARGIN=1, FUN=mean) }))
			rownames(gex$M) <- gex$genes$ProbeName
			colnames(gex$M) <- samples
			gex <- do.call("rbind", by(gex$M, INDICES=gpl[match(rownames(gex$M),gpl$SPOT_ID),"GENE_SYMBOL"], FUN=collapse_fun))
			# Keep only non-NA & non-"" rows
			gex <- gex[!(rownames(gex) == '' | is.na(rownames(gex))), ]
		}else{
			# Old code by FC - does not process RAW data
			gpl <- GEOquery::getGEO(geo_code, GSEMatrix = TRUE, getGPL = TRUE)
			labels <- Biobase::fData(gpl[[1]])
			gex <- Biobase::exprs(gpl[[1]])
			rownames(gex) <- labels$GENE_SYMBOL
			gex <- gex[-which(rownames(gex) == ''), ]
		}
	}

	# IGC
	else if(geo_code == "GSE2109"){
		## Note: IGC has >2k samples, many of which are not prostate cancer
		
		# Open the tarball(s)
		utils::untar(tarfile = rownames(supfiles))

		# Make sure to function in a working directory where the are no other tarballs present
		gz_files <- list.files()
		gz_files <- gz_files[grep(".gz", gz_files)]

		# Need phenodata as samples as from mixed cancers
		gset <- GEOquery::getGEO(geo_code, GSEMatrix = TRUE)
		# Subset to just prostate cancer samples	
		pca_gsm <- rownames(Biobase::pData(gset[[1]]))[grep("Prostate", Biobase::pData(gset[[1]])[,"title"])]
		gz_files_pca <- grep(paste0(pca_gsm, collapse="|"), gz_files, value=TRUE)

		if(pckg == "oligo"){
			## Subset to just prostate cancer as IGC covers in this GEO submission other cancer types also	
			# Read CEL
			gex <- oligo::read.celfiles(gz_files_pca)
			# Normalize background convolution of noise and signal using RMA (median-polish)
			gex <- oligo::rma(gex)
			# Extract expression matrix with probe ids
			gex <- oligo::exprs(gex)
			# Sanitize column names
			colnames(gex) <- gsub(".CEL.gz", "", colnames(gex))
			# Map the probes to gene symbols stored in curatedPCaData:::curatedPCaData_genes for hgu133a
			genes <- curatedPCaData:::curatedPCaData_genes_affy_hg_u133a_2[match(rownames(gex), curatedPCaData:::curatedPCaData_genes_affy_hg_u133a_2$affy_hg_u133a_2),"hgnc_symbol"]
			# Collapse probes that target the same gene
			gex <- do.call("rbind", by(gex, INDICES=genes, FUN=collapse_fun))
		}else if(pckg == "affy"){
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
		}
	}

  	# Sun et al. 
	else if(geo_code == "GSE25136"){
		# Open the tarball(s)
		utils::untar(tarfile = rownames(supfiles))

		# Make sure to function in a working directory where the are no other tarballs present
		gz_files <- list.files()
		gz_files <- gz_files[grep(".gz", gz_files)]

		if(pckg == "oligo"){
			# Read CEL
			gex <- oligo::read.celfiles(gz_files)
			# Normalize background convolution of noise and signal using RMA (median-polish)
			gex <- oligo::rma(gex)
			# Extract expression matrix with probe ids
			gex <- oligo::exprs(gex)
			# Sanitize column names
			colnames(gex) <- gsub(".CEL.gz", "", colnames(gex))
			# Map the probes to gene symbols stored in curatedPCaData:::curatedPCaData_genes for hgu133a
			genes <- curatedPCaData:::curatedPCaData_genes_affy_hg_u133a[match(rownames(gex), curatedPCaData:::curatedPCaData_genes_affy_hg_u133a$affy_hg_u133a),"hgnc_symbol"]
			# Collapse probes that target the same gene
			gex <- do.call("rbind", by(gex, INDICES=genes, FUN=collapse_fun))
		}else if(pckg == "affy"){	
			# Read Affymetrix MA
			Sun <- affy::ReadAffy()
			colnames(affy::exprs(Sun)) <- gsub(".gz|.CEL", "", colnames(Sun))

			# Careful not to mask 'rma' from 'affy' by the 'rma' from 'oligo'
			gex <- affy::rma(Sun)

			# Extracting .CEL and packaging names from the GEO-compatible sample names
			colnames(gex) <- gsub(".CEL.gz", "", colnames(affy::exprs(gex)))

			# Find gene annotations
			keys <- AnnotationDbi::mappedkeys(hgu133a.db::hgu133aGENENAME)
			nam <- names(as.character(hgu133a.db::hgu133aALIAS2PROBE)[match(rownames(gex), as.character(hgu133a.db::hgu133aALIAS2PROBE))])
			nam[is.na(nam)] <- "NA"
			# Collapse probes
			gex <- do.call("rbind", by(as.matrix(affy::exprs(gex)), INDICES = nam, FUN = collapse_fun))

			# Gene mapping and uniqueness via updateAnno-functionality
			gex <- curatedPCaData:::updateAnno(x=gex, main="hgnc_symbol", type="Aliases", collapse_fun=collapse_fun)
		}
	}

	# Taylor et al.
	# [HuEx-1_0-st] Affymetrix Human Exon 1.0 ST Array [probe set (exon) version]
	# Affymetrix Human Exon 1.0 ST Array [CDF: HuEx_1_0_st_v2_main_A20071112_EP.cdf]
	else if(geo_code == "GSE21032"){ # TODO: Alternative more specific accession code "GSE21034"
		# Open the tarball(s)
		utils::untar(tarfile = rownames(supfiles))
		# Make sure to function in a working directory where the are no other tarballs present
		gz_files <- list.files()
		gz_files <- gz_files[grep(".gz", gz_files)]
		# Subset to prostate cancer or normal patient samples
		gz_files <- grep("PCA|PAN", gz_files, value=TRUE)
		
		if(pckg == "oligo"){
			# Read CEL
			gex <- oligo::read.celfiles(gz_files)
			# Normalize background convolution of noise and signal using RMA (median-polish)
			gex <- oligo::rma(gex)
			# Extract annotated gene information 
			Biobase::featureData(gex) <- oligo::getNetAffx(gex, "transcript")
			genes <- unlist(lapply(Biobase::fData(gex)[,"geneassignment"], FUN = function(z) { strsplit(z, " // ")[[1]][2] }))
			# Extract expression matrix with probe ids
			gex <- oligo::exprs(gex)
			# Sanitize column names
			colnames(gex) <- unlist(lapply(colnames(gex), FUN=function(x) { grep("PCA|PAN", gsub(".CEL.gz", "", strsplit(x, "_")[[1]]), value=TRUE)} ))
			# Collapse the mapped probes under the same gene name
			gex <- do.call("rbind", by(gex, INDICES=genes, FUN=collapse_fun))
			# Omit duplicated entries that are replicated columns
			gex <- gex[,!duplicated(colnames(gex))]
		}
		else if(pckg == "affy"){
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
		}
	}

	# True et al.
	# GPL3834	FHCRC Human Prostate PEDB cDNA Array
	# GPL3836	FHCRC Human Prostate PEDB cDNA Array v3 (-> single sample only! 11th, GSM115769)
	else if(geo_code == "GSE5132"){
		# Open the tarball(s)
		utils::untar(tarfile = rownames(supfiles))
		# Make sure to function in a working directory where the are no other tarballs present
		gz_files <- list.files()
		gz_files <- gz_files[grep(".gz", gz_files)]

		if(pckg == "limma"){
			# Extract latest GPL annotations from GEO
			gpl <- GEOquery::getGEO(geo_code, GSEMatrix = FALSE, getGPL = TRUE)
			# Add meta Block column for mapping between GPL and genes read in with limma::read.maimages; GEO's GPL-files do not share identifiers with raw .gpr files
			gpl@gpls[[1]]@dataTable@table$Block <- rep(1:32, each=484)
			#gpl@gpls[[1]]@dataTable@table$Block <- rep(1:32, times=484)
			gpl@gpls[[2]]@dataTable@table$Block <- rep(1:32, each=420)			
			#gpl@gpls[[2]]@dataTable@table$Block <- rep(1:32, times=420)
			# Construct ID with xx_yy_zz where xx = Block, yy = Column, zz = Row (includes leading zeroes)
			gpl@gpls[[1]]@dataTable@table$ID <- paste(sprintf("%02d", gpl@gpls[[1]]@dataTable@table$Block), sprintf("%02d", gpl@gpls[[1]]@dataTable@table$Column), sprintf("%02d", gpl@gpls[[1]]@dataTable@table$Row), sep="_")
			gpl@gpls[[2]]@dataTable@table$ID <- paste(sprintf("%02d", gpl@gpls[[2]]@dataTable@table$Block), sprintf("%02d", gpl@gpls[[2]]@dataTable@table$Column), sprintf("%02d", gpl@gpls[[2]]@dataTable@table$Row), sep="_")
			# Gunzip .gpr files and read them via limma
			lapply(gz_files, FUN=GEOquery::gunzip)
			# FHCRC Human Prostate PEDB cDNA Array v3 vs v4
			#gpl_1 <- paste0(c("GSM115759", paste0("GSM11576", 0:8)), ".gpr")
			#gpl_2 <- grep(".gpr", list.files(), value=TRUE)
			#gpl_2 <- gpl_2[which(!gpl_2 %in% gpl_1)]
			gpl_1 <- grep(".gpr", list.files(), value=TRUE)
			gpl_1 <- gpl_1[!gpl_1 == "GSM115769.gpr"]
			gpl_2 <- "GSM115769.gpr"
			#gex <- limma::read.maimages(grep(".gpr", list.files(), value=TRUE), source="genepix")
			gex1 <- limma::read.maimages(gpl_1, source="genepix")
			gex2 <- limma::read.maimages(gpl_2, source="genepix")
			# Platform 2 has IDs in different format
			#gex2$genes$ID2 <- gex2$genes$ID
			#gex2$genes$ID <- paste(sprintf("%02d", gex2$genes$Block), sprintf("%02d", gex2$genes$Column), sprintf("%02d", gex2$genes$Row), sep="_")
			# Sort GPL data tables to be in same order as the read MA images
			#gpl@gpls[[1]]@dataTable@table <- gpl@gpls[[1]]@dataTable@table[match(gpl@gpls[[1]]@dataTable@table$ID, gex1$genes$ID),]
			#gpl@gpls[[2]]@dataTable@table <- gpl@gpls[[2]]@dataTable@table[match(gpl@gpls[[2]]@dataTable@table$ID, gex2$genes$ID),]
			# Assign gene names to the metadata
			gex1$genes$genename <- gpl@gpls[[1]]@dataTable@table[,"Related Gene Symbol"]
			gex2$genes$genename <- gpl@gpls[[2]]@dataTable@table[,"Hugo"]
			# Background correction separately for both platforms
			gex1 <- limma::backgroundCorrect(gex1, method="normexp", offset=1)
			gex2 <- limma::backgroundCorrect(gex2, method="normexp", offset=1)
			# Normalize with global loess
			gex1 <- limma::normalizeWithinArrays(gex1, method="loess")
			gex2 <- limma::normalizeWithinArrays(gex2, method="loess")
			# Average expression
			gex1 <- limma::avereps(gex1)
			gex2 <- limma::avereps(gex2)
			# Collapse for gene names
			gex1 <- do.call("rbind", by(as.matrix(gex1), INDICES=gex1$genes$genename, FUN=collapse_fun))
			#gex2 <- do.call("rbind", by(as.matrix(gex2), INDICES=gex2$genes$genename, FUN=collapse_fun))
			name <- colnames(gex2) # Store single sample name which would be dropped otherwise
			gex2 <- as.matrix(by(as.matrix(gex2), INDICES=gex2$genes$genename, FUN=collapse_fun))
			colnames(gex2) <- name
			# Omit unannotated collapsed probes
			gex1 <- gex1[which(!rownames(gex1) %in% c("", NA)),,drop=FALSE]
			gex2 <- gex2[which(!rownames(gex2) %in% c("", NA)),,drop=FALSE]
			# Median impute the genes not measured in v3 in order to retain full dimension of 3615 of v4 array instead of downtoning to 2727 present in v3
			# This imputation happens only for a single paired sample:
			gex2 <- gex2[match(rownames(gex1),rownames(gex2)),,drop=FALSE]
			rownames(gex2) <- rownames(gex1)
			for(r in which(is.na(gex2))){
				gex2[r,] <- median(gex1[r,])
			}
			# Merge data based on intersecting genes
			# 11th sample is a special case for v3 array and must be placed separately
			gex <- cbind(gex1[,1:10], GSM115769 = gex2, gex1[,11:31])			
			
			# Channels were swapped in samples #3,4, 7 and 8, 12, 14, 15, 17, 18, 20, 22, 23, 24, 26, 30, 31
			# In these tumor vs normal log ratio is inverted
			for(i in c(3,4,7,8,12,14,15,17,18,20,22,23,24,26,30,31)){
				gex[,i] <- -gex[,i]
			}
			# 11th sample (GSM115769) was oddly processed with v3 array, while samples 1-10, 12-32 were processed with v4 array
		}
		else if(pckg == "other"){
			# As mentioned in the clinical section this data has been split in two datasest, one of 31 samples and one of just 1 sample
			gset <- GEOquery::getGEO(geo_code, GSEMatrix = TRUE, getGPL = TRUE)
			labels1 = Biobase::fData(gset[[1]])
			labels2 = Biobase::fData(gset[[2]])

			# First part of the data
			gex1 <- Biobase::exprs(gset[[1]])
			rownames(gex1) = labels1$"Related Gene Symbol"
			gex1 = gex1[-which(rownames(gex1) == ''), ]
			gex1 = aggregate(gex1, by = list(rownames(gex1)), mean, na.rm = T)
			rownames(gex1) = gex1[, 1]
			gex1 = gex1[, -1]
			# Second part of the data
			gex2 <- Biobase::exprs(gset[[2]])
			rownames(gex2) = labels2$Hugo
			gex2 = gex2[-which(rownames(gex2) == ''), , drop = F]
			gex2 = cbind(gex2, 1)
			gex2 = aggregate(gex2, by = list(rownames(gex2)), mean, na.rm = T)
			rownames(gex2) = gex2[, 1]
			gex2 = gex2[, -1]
			gex2 = gex2 [, -2, drop = F]

			# Intersect to common genes
			common_genes = intersect(rownames(gex1), rownames(gex2))
			gex1 = gex1[is.element(rownames(gex1), common_genes), ,drop = F]
			gex2 = gex2[is.element(rownames(gex2), common_genes), ,drop = F]

			# the two datasets are merged respecting the order of the GEO sample IDs
			if(identical(rownames(gex1), rownames(gex2))) gex = cbind(gex1[,1:10], gex2[,1], gex1[,11:31])
			# the appropriate name is used for the new column
			colnames(gex)[11] = colnames(gex2)
		}
	}
  
	# Wallace et al.	
	else if(geo_code == "GSE6956"){
  
		unmatched_healty_tissue <- c('GSM160418', 'GSM160419', 'GSM160420', 'GSM160421', 'GSM160422', 'GSM160430')
		# Open the tarball(s)
		utils::untar(tarfile = rownames(supfiles))
		# Make sure to function in a working directory where the are no other tarballs present
		gz_files <- list.files()
		gz_files <- gz_files[grep(".gz", gz_files)]
	
		if(pckg == "oligo"){
			# Read CEL
			gex <- oligo::read.celfiles(gz_files)
			# Normalize background convolution of noise and signal using RMA (median-polish)
			gex <- oligo::rma(gex)
			# Extract expression matrix with probe ids
			gex <- oligo::exprs(gex)
			# Sanitize column names
			colnames(gex) <- gsub(".CEL.gz", "", colnames(gex))
			# Map the probes to gene symbols stored in curatedPCaData:::curatedPCaData_genes for hgu133a
			genes <- curatedPCaData:::curatedPCaData_genes_affy_hg_u133a_2[match(rownames(gex), curatedPCaData:::curatedPCaData_genes_affy_hg_u133a_2$affy_hg_u133a_2),"hgnc_symbol"]
			# Collapse probes that target the same gene
			gex <- do.call("rbind", by(gex, INDICES=genes, FUN=collapse_fun))
		}
		else if(pckg == "affy"){
			# Read CEL files and name columns accordingly
			wallace <- affy::ReadAffy()
			colnames(affy::exprs(wallace)) <- gsub(".gz|.CEL", "", colnames(wallace))
			# RMA normalization
			gex <- affy::rma(wallace)
			# Removing .CEL and packaging names from the GEO-compatible sample names
			colnames(gex) <- gsub(".CEL.gz", "", colnames(affy::exprs(gex)))
			keys <- AnnotationDbi::mappedkeys(hgu133a.db::hgu133aGENENAME)
			nam <- names(as.character(hgu133a.db::hgu133aALIAS2PROBE)[match(rownames(gex), as.character(hgu133a.db::hgu133aALIAS2PROBE))])
			nam[is.na(nam)] <- "NA"
			# Collapse probes
			gex <- do.call("rbind", by(as.matrix(affy::exprs(gex)), INDICES=nam, FUN=collapse_fun))
			# Return numeric matrix

			## TDL: Store (unmatched) healthy tissues, can be useful for e.g. assessing tumor purity assessment where normal samples are required
			#gex <- gex[, !is.element(colnames(gex), unmatched_healty_tissue)]
		}
	}
  
	# Wang et al.  
	else if(geo_code == "GSE8218"){
		# Open the tarball(s)
		utils::untar(tarfile = rownames(supfiles))

		# Make sure to function in a working directory where the are no other tarballs present
		gz_files <- list.files()
		gz_files <- gz_files[grep(".gz", gz_files)]

		if(pckg == "oligo"){
			# Read CEL
			gex <- oligo::read.celfiles(gz_files)
			# Normalize background convolution of noise and signal using RMA (median-polish)
			gex <- oligo::rma(gex)
			# Extract expression matrix with probe ids
			gex <- oligo::exprs(gex)
			# Sanitize column names
			colnames(gex) <- gsub(".CEL.gz", "", colnames(gex))
			# Map the probes to gene symbols stored in curatedPCaData:::curatedPCaData_genes for hgu133a
			genes <- curatedPCaData:::curatedPCaData_genes_affy_hg_u133a[match(rownames(gex), curatedPCaData:::curatedPCaData_genes_affy_hg_u133a$affy_hg_u133a),"hgnc_symbol"]
			# Collapse probes that target the same gene
			gex <- do.call("rbind", by(gex, INDICES=genes, FUN=collapse_fun))
		}else if(pckg == "affy"){
			# Read Affymetrix MA
			wang <- affy::ReadAffy()
			colnames(affy::exprs(wang)) <- gsub(".gz|.CEL", "", colnames(wang))

			# Careful not to mask 'rma' from 'affy' by the 'rma' from 'oligo'
			gex <- affy::rma(wang)

			# Extracting .CEL and packaging names from the GEO-compatible sample names
			colnames(gex) <- gsub(".CEL.gz", "", colnames(affy::exprs(gex)))
			
			## !!
			## TDL:
			## NOTE! These probe names are wrong, the correct array is v2 not ordinary hg u133a
			## !!
			
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
		}
	}
  	
  	# Weiner et al.
  	else if(geo_code == "GSE157548"){
		# Open the tarball(s)
		utils::untar(tarfile = rownames(supfiles))

		# Make sure to function in a working directory where the are no other tarballs present
		gz_files <- list.files()
		gz_files <- gz_files[grep(".gz", gz_files)]
  		if(pckg == "oligo"){
			# Read CEL
			#gex <- oligo::read.celfiles(gz_files)
			# Too large for a high-end desktop: 
			#> Error: cannot allocate vector of size 40.9 Gb 		
			#gex <- oligo::read.celfiles(gz_files[1:400])
			gex <- oligo::read.celfiles(gz_files)
			# Normalize background convolution of noise and signal using RMA (median-polish)
			gex <- oligo::rma(gex)
			# Extract annotated gene information 
			Biobase::featureData(gex) <- oligo::getNetAffx(gex, "transcript")
			genes <- unlist(lapply(Biobase::fData(gex)[,"geneassignment"], FUN = function(z) { strsplit(z, " // ")[[1]][2] }))
			# Extract expression matrix with probe ids
			gex <- oligo::exprs(gex)
			# Sanitize column names
			colnames(gex) <- gsub(".CEL.gz", "", colnames(gex))
			# Collapse the mapped probes under the same gene name
			gex <- do.call("rbind", by(gex, INDICES=genes, FUN=collapse_fun))
			# Sanitize sample names to only contain GSM####
			colnames(gex) <- unlist(lapply(colnames(gex), FUN=function(x){ strsplit(x, "_")[[1]][1] }))
			# Store intensities up to 6th decimal point to conserve space
			gex <- round(gex, 6)
  		}  		
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

	# Remove rows with no names or redundancy
	if(any(is.na(rownames(gex))) | any(rownames(gex) %in% c(""))){
		gex <- gex[-which(is.na(rownames(gex)) | rownames(gex) %in% c("")),]
	}
	# Sort genes to alphabetic order for consistency
	gex <- gex[order(rownames(gex)),]
	# Cast to matrix type, remove empty rows/columns and return gene expression matrix
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
#' @noRd
#' @keywords internal
generate_cna_geo <- function(
	geo_code = c(
		"GSE54691",	# Hieronymus et al.
		"GSE21035"	# Taylor et al. (aCGH only GSE-subset)
	),
	file_directory, 
	filter_regex,
	cleanup = TRUE, 
	...
){
	if(!missing(file_directory)) here::set_here(file_directory)
	# Supplementary files include the raw CEL files, possibility to use regular expression filter
	if(missing(filter_regex)){
		supfiles <- GEOquery::getGEOSuppFiles(geo_code)
	}else{
		supfiles <- GEOquery::getGEOSuppFiles(geo_code, filter_regex = filter_regex)
	}
  
	# Open the tarball(s)
	utils::untar(tarfile = rownames(supfiles))
  
	##
	# same rCGH pipeline applied to datasets:
	# Taylor et al. : GPL4091	Agilent-014693 Human Genome CGH Microarray 244A (Feature number version)
	# Hieronymus et al. : GPL8737	Agilent-021529 Human CGH Whole Genome Microarray 1x1M (G4447A) (Probe Name version)
	##
	if(geo_code %in% c(
		"GSE21035", # Taylor
		"GSE54691" # Hieronymus
	)){
		# Extract sample names in Taylor et al. and format to PCA####, for this need to access GSM-level metadata
		if(geo_code == "GSE21035"){
			# Access metadata based on GSM####
			samplenames <- lapply(gsub(".txt.gz", "", grep("GSM", list.files(), value=TRUE)), FUN=function(x){ GEOquery::getGEO(x, GSEMatrix=FALSE, getGPL=FALSE)@header$title })
			# Parse to just PCA#### from the title, splitting on whitespace character ' '
			samplenames <- unlist(lapply(samplenames, FUN=function(x) { strsplit(x, " ")[[1]][3] }))
			# File names as they are
			filenames <- grep("GSM", list.files(), value=TRUE)
			# Omit cell lines
			filenames <- filenames[grep("PCA", samplenames)]
			samplenames <- samplenames[grep("PCA", samplenames)]
		}
		# Hieronymus et al names based on GSM####
		else{
			# File names end in suffix .txt.gz, sample names have this string replaced
			samplenames <- gsub(".txt.gz", "", grep("GSM", list.files(), value=TRUE))
			filenames <- grep("GSM", list.files(), value=TRUE)
		}
		# For now, the package 'rCGH' has to be available in the workspace,
		require(rCGH)
		# otherwise below functions will fail on e.g. rCGH::adjustSignal and when trying to find 'hg18'
		# Read in Agilent 2-color data
		cna <- lapply(1:length(filenames), FUN = function(i) { 
			try({
				cat("\nProcessing: ",filenames[i],"\n") 
				rCGH::readAgilent(filenames[i], genome = "hg38", sampleName = samplenames[i]) 
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
		## Samplenames picked prior to this
		# Remove additional suffixes from sample names
		#cna <- lapply(cna, FUN = function(z){ 
		#	try({
		#		if(!"rCGH-Agilent" %in% class(z)){
		#			stop("This function is intended for Agilent aCGH analyzed with rCGH R Package (class \'rCGH-Agilent\')")		
		#		}
		#		# e.g. transform "GSM525575.txt|.gz" -> "GSM525575"
		#		z@info["sampleName"] <- gsub(pattern = ".gz|.txt", replacement = "", z@info["fileName"])
		#		z
		#	}) 
		#})
		# Save sample names separately (of 'length(cna)')
		#samplenames <- unlist(lapply(cna, FUN = function(z) { z@info["sampleName"] }))
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
		cna <- cna[!is.na(rownames(cna)),]
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
	# Sort genes to alphabetic order for consistency
	cna <- cna[order(rownames(cna)),]
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
#' @noRd
#' @keywords internal
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
                      "prad_su2c_2019_gistic", #Abida et al. CNA
                      "prad_broad_2013_cna" #BACA
                      ), # for cgdsr calls, platform and dataset specific string
  caseList = c("prad_tcga_pub_sequenced", # TCGA
               "prad_mskcc_sequenced", # Taylor et al.
               "prad_broad_sequenced", # PRAD Broad
               "prad_eururol_2017_sequenced", # PRAD Eururol
               "prad_su2c_2019_sequenced",#Abida et al.
               "prad_broad_2013_sequenced"#BACA
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
#' @noRd
#' @keywords internal
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
			c("HOMDEL", "HETLOSS", "", "LOW_GAIN", "HIGH_AMP")[z+3]
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
			if(is.na(oncop[mut[i,"gene_symbol"], mut[i,"case_id"]]) | oncop[mut[i,"gene_symbol"], mut[i,"case_id"]] == ""){
				oncop[mut[i,"gene_symbol"], mut[i,"case_id"]] <- mut[i,"mutation_type"]
			}else{
				oncop[mut[i,"gene_symbol"], mut[i,"case_id"]] <- paste(oncop[mut[i,"gene_symbol"], mut[i,"case_id"]], mut[i,"mutation_type"], sep=";")
			}
		}
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

#' Download mutation data from cBioPortal using cgsdr package
#'
#' @param study_id TODO
#' @param genes TODO
#' @param delay TODO
#' @param splitsize TODO
#' @param verb TODO
#' @examples
#' ren_mut=generate_cgdsr_mut("ren",genes=curatedPCaData::curatedPCaData_genes$hgnc_symbol)
#' barbieri_mut=generate_cgdsr_mut("barbieri",genes=curatedPCaData::curatedPCaData_genes$hgnc_symbol)
#' @noRd
#' @keywords internal
generate_cgdsr_mut <- function(
  study_id, # tcga, taylor, barbieri, ren, or abida
  genes = sort(unique(curatedPCaData:::curatedPCaData_genes$hgnc_symbol)),
  delay = 0.05,
  splitsize = 100,
  verb = TRUE
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
    )
  }))
  
  list(mut)
  #}
}


#' Download generic 'omics data from cBioPortal using the cbioportaldata package
#' 
#' @param profile character string of variable
#' @param caseList character string of patient IDs for that genetic profile
#' @examples
#' cna.gistic_abida <- curatedPCaData:::generate_cbioportaldata("prad_su2c_2019","cna")
#' gex.relz_abida <- curatedPCaData:::generate_cbioportaldata("prad_su2c_2019","gex")
#' @noRd
#' @keywords internal
generate_cbioportaldata <- function(caselist,profile){
  prof=cBioPortalData::cBioDataPack(caselist,ask = FALSE)
  
  if(profile=="cna"){
    if (caselist == "prad_broad"){
      res=prof[["cna"]]
      a=rowData(res)[match(rownames(res),rownames(rowData(res))),"Hugo_Symbol"]
      res2=assay(res)
      res2<-res2[!is.na(a),]
      a=a[!is.na(a)]
      rownames(res2) <- a
      symbols <- curatedPCaData:::curatedPCaData_genes[match(rownames(res2), curatedPCaData:::curatedPCaData_genes[,"hgnc_symbol"]),"hgnc_symbol"]
      res2 <- res2[!is.na(symbols),]
      symbols <- symbols[!is.na(symbols)]
      rownames(res2) <- symbols
      colnames(res2)<-gsub("-",".",colnames(res2))
      res2<-res2[rowSums(is.na(res2)) != ncol(res2), ]
      return(res2)}
    
    else if(caselist=="prad_eururol_2017"){
      res=prof[["cna"]]
      res2=assay(res)
      symbols <- curatedPCaData:::curatedPCaData_genes[match(rownames(res2), curatedPCaData:::curatedPCaData_genes[,"hgnc_symbol"]),"hgnc_symbol"]
      res2 <- res2[!is.na(symbols),]
      symbols <- symbols[!is.na(symbols)]
      rownames(res2) <- symbols
      res2<-res2[rowSums(is.na(res2)) != ncol(res2), ]
      return(res2)
      
    }
    
    else if(caselist=="prad_su2c_2019"){
      res=metadata(prof)$cna
      res=as.data.frame(res)
      res2=res[!(duplicated(res$Hugo_Symbol)|duplicated(res$Hugo_Symbol, fromLast=TRUE)),, drop=FALSE]
      res2=res2[!grepl("-AS1", res2$Hugo_Symbol),]
      rownames(res2)<-res2$Hugo_Symbol
      res2<-res2[,-1]
      symbols <- curatedPCaData:::curatedPCaData_genes[match(rownames(res2), curatedPCaData:::curatedPCaData_genes[,"hgnc_symbol"]),"hgnc_symbol"]
      
      res2 <- res2[!is.na(symbols),]
      symbols <- symbols[!is.na(symbols)]
      rownames(res2) <- symbols
      colnames(res2)<-gsub("-",".",colnames(res2))
      res2<-res2[rowSums(is.na(res2)) != ncol(res2), ]
      return(res2)
    }
    
    else if (caselist=="prad_mskcc"){
      study=downloadStudy("prad_mskcc")
      ut <- untarStudy(study[[1]])
      res2=rio::import(paste0(ut,"/prad_mskcc/","data_cna.txt"))
      res2<-res2[,-2]
      rownames(res2)<-make.names(res2$Hugo_Symbol,unique = TRUE)
      symbols <- curatedPCaData:::curatedPCaData_genes[match(rownames(res2), curatedPCaData:::curatedPCaData_genes[,"hgnc_symbol"]),"hgnc_symbol"]
      res2 <- res2[!is.na(symbols),]
      symbols <- symbols[!is.na(symbols)]
      rownames(res2) <- symbols
      res2<-res2[,-1]
      res2<-res2[rowSums(is.na(res2)) != ncol(res2), ]
      
      same_barcode=colnames(res2)[grepl("PCA", colnames(res2))]
      ind=which(colnames(res2) %in% same_barcode=="TRUE")
      res3<-res2[, c(ind)]
      
      return(res3)}
    
    else if (caselist == "prad_broad_2013"){
      res=prof[["cna"]]
      res2=assay(res)
      symbols <- curatedPCaData:::curatedPCaData_genes[match(rownames(res2), curatedPCaData:::curatedPCaData_genes[,"hgnc_symbol"]),"hgnc_symbol"]
      res2 <- res2[!is.na(symbols),]
      symbols <- symbols[!is.na(symbols)]
      rownames(res2) <- symbols
      colnames(res2)<-gsub("-",".",colnames(res2))
      res2<-res2[rowSums(is.na(res2)) != ncol(res2), ]
      return(res2)
      
    }
  }
  
  if(profile=="mut"){
    ragexp=prof[["mutations"]]
    if(caselist=="prad_eururol_2017"){
      symbols <- curatedPCaData:::curatedPCaData_genes[match(rownames(ragexp), curatedPCaData:::curatedPCaData_genes[,"hgnc_symbol"]),"hgnc_symbol"]
      ragexp2 <- ragexp[!is.na(symbols),]
      symbols <- symbols[!is.na(symbols)]
      rownames(ragexp2) <- symbols
      return(ragexp2)
      
    }else if (caselist == "prad_broad" || caselist=="prad_su2c_2019"){
      colnames(ragexp)<-gsub("-",".",colnames(ragexp))
      symbols <- curatedPCaData:::curatedPCaData_genes[match(rownames(ragexp), curatedPCaData:::curatedPCaData_genes[,"hgnc_symbol"]),"hgnc_symbol"]
      ragexp2 <- ragexp[!is.na(symbols),]
      symbols <- symbols[!is.na(symbols)]
      rownames(ragexp2) <- symbols
      colnames(ragexp2)<-gsub("-",".",colnames(ragexp2))
      return(ragexp2)
      
    }else if (caselist=="prad_mskcc"){
      same_barcode=colnames(ragexp)[grepl("PCA", colnames(ragexp))]
      ind=which(colnames(ragexp) %in% same_barcode=="TRUE")
      ragexp2<-ragexp[, c(ind)]
      symbols <- curatedPCaData:::curatedPCaData_genes[match(rownames(ragexp2), curatedPCaData:::curatedPCaData_genes[,"hgnc_symbol"]),"hgnc_symbol"]
      ragexp3 <- ragexp2[!is.na(symbols),]
      symbols <- symbols[!is.na(symbols)]
      rownames(ragexp3) <- symbols
      
      return(ragexp3)
      
    }else if (caselist=="prad_broad_2013"){
      symbols <- curatedPCaData:::curatedPCaData_genes[match(rownames(ragexp), curatedPCaData:::curatedPCaData_genes[,"hgnc_symbol"]),"hgnc_symbol"]
      ragexp2 <- ragexp[!is.na(symbols),]
      symbols <- symbols[!is.na(symbols)]
      rownames(ragexp2) <- symbols
      colnames(ragexp2)<-gsub("-",".",colnames(ragexp2))
      return(ragexp2)
      
    }
  }
  if(profile=="gex"){
    if(caselist=="prad_eururol_2017"){
      gex=prof[["mrna_seq_rpkm_zscores_ref_all_samples"]]
      gex2=assay(gex)
      symbols <- curatedPCaData:::curatedPCaData_genes[match(rownames(gex2), curatedPCaData:::curatedPCaData_genes[,"hgnc_symbol"]),"hgnc_symbol"]
      gex2 <- gex2[!is.na(symbols),]
      symbols <- symbols[!is.na(symbols)]
      rownames(gex2) <- symbols
      gex2<-gex2[rowSums(is.na(gex2)) != ncol(gex2), ]
      return(gex2)
    }else if(caselist=="prad_broad"){
      gex=prof[["mrna_agilent_microarray_zscores_ref_all_samples"]]
      gex2=assay(gex)
      gex2=gex2[rowSums(is.na(gex2)) != ncol(gex2), ]
      symbols <- curatedPCaData:::curatedPCaData_genes[match(rownames(gex2), curatedPCaData:::curatedPCaData_genes[,"hgnc_symbol"]),"hgnc_symbol"]
      gex2 <- gex2[!is.na(symbols),]
      symbols <- symbols[!is.na(symbols)]
      rownames(gex2) <- symbols
      colnames(gex2)<-gsub("-",".",colnames(gex2))
      gex2<-gex2[rowSums(is.na(gex2)) != ncol(gex2), ]
      return(gex2)
    }else if(caselist=="prad_su2c_2019"){
      gex=metadata(prof)$mrna_seq_fpkm_polya_zscores_ref_all_samples
      gex=as.data.frame(gex)
      gex2=gex[!(duplicated(gex$Hugo_Symbol)|duplicated(gex$Hugo_Symbol, fromLast=TRUE)),, drop=FALSE]
      gex2=gex2[!grepl("-AS1", gex2$Hugo_Symbol),]
      rownames(gex2)<-gex2$Hugo_Symbol
      gex2<-gex2[,-1]
      symbols <- curatedPCaData:::curatedPCaData_genes[match(rownames(gex2), curatedPCaData:::curatedPCaData_genes[,"hgnc_symbol"]),"hgnc_symbol"]
      
      gex2 <- gex2[!is.na(symbols),]
      symbols <- symbols[!is.na(symbols)]
      rownames(gex2) <- symbols
      colnames(gex2)<-gsub("-",".",colnames(gex2))
      gex2<-gex2[rowSums(is.na(gex2)) != ncol(gex2), ]
      return(gex2)
    }
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
#' @noRd
#' @keywords internal
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

#' Download GDC processed data via xenabrowser.net Santa-Cruz interface
#'
#' @noRd
#' @keywords internal
generate_xenabrowser <- function(
	id = "TCGA-PRAD", # Study ID (by expectation TCGA's Prostate Adenocarcinoma
	type = c("gex", "cna", "mut", "clinical"), # First instance of vector is used to determine what is extracted
	# Function for collapsing rows for identical gene symbols
	collapse_fun = function(z) {apply(z, MARGIN = 2, FUN = stats::median)},	
	# If Sample IDs should be truncated down to Patient ID level (leave out last segment of the '-' or '.' separators)
	truncate = TRUE,
	# Number of digits to store for the data object; for large matrices this may be required to stay beneath 100 MB, or to get rid of insignificant digits
	digits,
	# If intermediate files ought to be removed
	cleanup = TRUE,
	...
){
	# Small internal function to assist with the downloads from xenabrowser.net
	.xenabrowserDownload <- function(url, gz=TRUE){
		# Pick the filename from the end of the URL
		filename <- strsplit(url, "/")
		filename <- filename[[1]][[length(filename[[1]])]]
		# Download file into parsed *.tsv.gz 
		utils::download.file(url=url, destfile=filename)
		# if .gz, gunzip the files open
		if(gz) GEOquery::gunzip(filename, overwrite=TRUE)
		gsub(".gz", "", filename)
	}
	# Cases (typically TCGA-PRAD)
	if(id == "TCGA-PRAD"){
		# Release Mid 2019ish
		urls = list(
			"htseq.fpkm" = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-PRAD.htseq_fpkm-uq.tsv.gz",
			"gistic" = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-PRAD.gistic.tsv.gz",
			"mutect2" = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-PRAD.mutect2_snv.tsv.gz",
			"phenotype" = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-PRAD.GDC_phenotype.tsv.gz",
			"os" = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-PRAD.survival.tsv",
			"genemap" = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v22.annotation.gene.probeMap"
		)
		# Gene expression
		if(type == "gex"){
			# Download the FPKM-UQ values processed via HTSeq
			file <- .xenabrowserDownload(urls[["htseq.fpkm"]])
			dat <- read.table(file, sep="\t", header=TRUE, row.names=1)
			if(cleanup) file.remove(file)
			# Download the gene mapping file
			file <- .xenabrowserDownload(urls[["genemap"]], gz=FALSE)
			map <- read.table(file, sep="\t", header=TRUE, row.names=1)			
			if(cleanup) file.remove(file)
			# Rearrange to match the ordering of the data matrix
			map <- map[match(rownames(dat), rownames(map)),]
			#> all(rownames(map) == rownames(dat))
			#[1] TRUE
			# Reassign gene names from ENSEMBL gene ids to gene symbols collapsing identical names and arranging to alphabetic order
			dat <- do.call("rbind", by(dat, INDICES=map[,"gene"], FUN=collapse_fun))
			dat <- dat[order(rownames(dat)),]
			#> dim(dat)
			#[1] 58387   551
			#
			# If column names should truncate last segment (sample -> patient id level truncation)
			if(truncate>=1){
				print("Truncating to 3 dot-separated names...")
				# Remove fourth element separated by '.'
				colnames(dat) <- unlist(lapply(colnames(dat), FUN=function(x) { paste(strsplit(x, ".", fixed=TRUE)[[1]][1:3], collapse=".") }))
			# Lesser truncation with suffix '.01A' -> '.01' as in cBio
			}else if(truncate>=0){
				print("Substituting '.01A|.01B' with '.01' ...")
				# Sub '.01A" with the default ".01"
				colnames(dat) <- gsub(".01A|.01B", ".01", colnames(dat))
			}
		# Copy number alterations (discretized by GISTIC)
		}else if(type == "cna"){
			# Download the copy number alteration values processed via GISTIC
			file <- .xenabrowserDownload(urls[["gistic"]])
			dat <- read.table(file, sep="\t", header=TRUE, row.names=1)
			if(cleanup) file.remove(file)
			# Download the gene mapping file
			file <- .xenabrowserDownload(urls[["genemap"]], gz=FALSE)
			map <- read.table(file, sep="\t", header=TRUE, row.names=1)			
			if(cleanup) file.remove(file)
			# Rearrange to match the ordering of the data matrix
			map <- map[match(rownames(dat), rownames(map)),]
			# > all(rownames(map) == rownames(dat))
			#[1] TRUE
			#> dim(dat)
			#[1] 19729   502
			#
			## Smaller dimension than for gex as expected
			dat <- do.call("rbind", by(dat, INDICES=map[,"gene"], FUN=collapse_fun))
			dat <- dat[order(rownames(dat)),]
			# If column names should truncate last segment (sample -> patient id level truncation)
			if(truncate>=1){
				print("Truncating to 3 dot-separated names...")
				# Remove fourth element separated by '.'
				colnames(dat) <- unlist(lapply(colnames(dat), FUN=function(x) { paste(strsplit(x, ".", fixed=TRUE)[[1]][1:3], collapse=".") }))
			# Lesser truncation with suffix '.01A' -> '.01' as in cBio
			}else if(truncate>=0){
				print("Substituting '.01A|.01B' with '.01' ...")
				# Sub '.01A" with the default ".01"
				colnames(dat) <- gsub(".01A|.01B", ".01", colnames(dat))
			}
		# Small mutations (SNV / INDELs called by Mutect2)
		}else if(type == "mut"){
			# Download the MuTect2 somatic mutation calls
			file <- .xenabrowserDownload(urls[["mutect2"]])
			# Suitable for RaggedExperiment style data storage
			dat <- read.table(file, sep="\t", header=TRUE)
			if(cleanup) file.remove(file)
			tcga_mut<-dat[,c(3:5,1,2,6:11)]
			colnames(tcga_mut)[1:3]=c("seqnames","start","end")
			tcga_mut$Sample_ID<-gsub("-",".",tcga_mut$Sample_ID)
			tcga_mut$Sample_ID<-gsub("01A","01",tcga_mut$Sample_ID)
			names(tcga_mut)[names(tcga_mut) == 'effect'] <- "Variant_Classification"
			#a=subset(tcga_mut, Sample_ID %in% colnames(mae_tcga[["gex.fpkm"]]))
			GRL <- makeGRangesListFromDataFrame(tcga_mut, split.field = "Sample_ID",
			                                    names.field = "gene",keep.extra.columns = TRUE)
			ragexp_tcga=RaggedExperiment::RaggedExperiment(GRL)
			return(ragexp_tcga)
		# Clinical data matrix construction
		}else if(type == "clinical"){
			# Generic phenotype information
			file <- .xenabrowserDownload(urls[["phenotype"]])
			phenotype <- read.table(file, sep="\t", header=TRUE, row.names=1, quote="#")
			if(cleanup) file.remove(file)
			# Overall Survival
			file <- .xenabrowserDownload(urls[["os"]], gz=FALSE)
			os <- read.table(file, sep="\t", header=TRUE, row.names=1, quote="#")
			if(cleanup) file.remove(file)
			# Combine the two			
			dat <- cbind(phenotype, os[match(rownames(phenotype), rownames(os)),])
			# .11A are healthy samples (GEX)
			#rownames(dat) <- gsub(".01A|.11A|01B", "", gsub("-", ".", rownames(dat)))
		# Unknown data type
		}else{
			stop(paste("Invalid query type for xenabrowser:", type))
		}
	}
	# Round to certain digits if requested
	if(!missing(digits)) dat <- round(dat, digits)
	# Return the processed dat
	dat
}

