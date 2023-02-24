###
#
# Variable functions helping exporting of data/information from curatedPCaData
#
###

#' Export ExperimentHub-friendly metadata csv prototype
#'
#' This function greps over MAE-objects available via curatedPCaData and starts constructing ExperimentHub-compatible metadata.csv, first as a data.frame
#'
#' @noRd
#' @keywords internal
export_metadata <- function(
	timestamp, # Time stamp to append in file names
	path = getwd(), # Location on drive to save the files into
	...
){
	# If user does not provide a time stamp, construct one from the system time ('%Y%m%d'-format)
	if(missing(timestamp)){
		timestamp <- gsub("-", "", substr(Sys.time(), 0, 10))
	}

	# List MAEs available in curatedPCaData
	# Get available maes in curatedPCaData package and specify a path to the folder where to save
	maes <- grep("mae", utils::data(package="curatedPCaData")$result, value=TRUE)
	
	# Metadata fields that are pre-known
	biocv <- "3.16" #BioConductor version
	species <- "Homo sapiens" # Human samples only

	# Listing for description fields
	descriptions <- c(
		# CNA copy number alteration fields with differing levels of information
		cna.gistic = "Copy number alteration GISTIC",
		cna.logr = "Copy number alteration log-ratios", 
		# GEX - Gene expression fields with different processing pipelines
		gex.relz = "gene expression", 
		gex.logq = "gene expression",
		gex.rma = "gene expression",
		gex.logr = "gene expression",
		gex.rsem.log = "gene expression", 
		# MUT - Mutation data
		mut = "mutation", 
		cibersort = "cibersort deconvolution", 
		xcell = "xcell deconvolution",
		epic = "epic deconvolution", 
		quantiseq = "quantiseq deconvolution", 
		mcp = "mcp deconvolution",
		estimate = "estimate devonvolution",
		scores= "Various gene expression risk and marker scores"
	)
	# Format an empty df with required fields
	metadata <- data.frame(
		Title = NA_character_,	
		Description = NA_character_,	
		BiocVersion = NA_character_,	
		Genome = NA_character_,	
		SourceType = NA_character_,	
		SourceUrl = NA_character_,
		SourceVersion = NA_character_,	
		Species = NA_character_,
		TaxonomyId = NA_character_,
		Coordinate_1_based = NA_character_,
		DataProvider = NA_character_,
		Maintainer = NA_character_,
		RDataClass = NA_character_,
		DispatchClass = NA_character_,
		ResourceName = NA_character_,
		RDataPath = NA_character_,
		Tags = NA_character_
	)
	for(i in 1:length(maes)){  #length(maes)
		mae <- maes[i]  # Get specific mae & separate the study name for future naming
		study <- gsub("mae_","",mae)
		ex <- names(experiments(get(mae))) # Get the names of object ("experiments") in mae

		# Go through all "experiments" and save the object into its own .RDs file
		#lapply(ex,FUN=function(x){

		for(x in ex){
			oname <- paste0(study,"_",x,"_",timestamp)  # object name
			fname <- paste0(oname,".Rds")  # file name

			assign(oname,get(mae)[[x]]) # Assign the object to the specific name
			save(list=c(oname), file=paste0(path,"/",fname))

			# Add info to metadata
			descrp <- paste(oname,descriptions[x],"data of",study,"cohort in curatedPCaData package",sep=" ")

		# Add with <<- to modify global variable from lapply
			metadata <- rbind(metadata, 
				c(Title = oname,
				Description = descrp,
				BiocVersion = biocv,
				Species = species,
				RDataClass = class(get(oname))[1],
				DispatchClass = "Rds",
				ResourceName = fname,
				RDataPath = paste0("curatedPCaData/", fname),
				Tags = x
			)
		}

		# Save clinical info
		oname <- paste0(study,"_colData_",vtag)  # object name
		fname <- paste0(oname,".Rds")  # file name

		assign(oname,MultiAssayExperiment::colData(get(mae))) # Assign the object to the spesific name
		save(list=c(oname), file=paste0(path,"/",fname))

		# Add info to metadata
		descrp <- paste(oname,"Clinical metadata (colData-slot) of",study,"cohort in curatedPCaData package",sep=" ")

		metadata_raw <- rbind(metadata_raw,
			c(Title=oname,
			Description=descrp,
			BiocVersion=biocv,
			Species=species,
			RDataClass=class(get(oname))[1],
			DispatchClass="Rds",
			ResourceName=fname,
			RDataPath=paste0("curatedPCaData/",fname),
			Tags="clinical")
		)

		# Save sampleMap
		oname <- paste0(study,"_sampleMap_",vtag)  # object name
		fname <- paste0(oname,".Rds")  # file name

		assign(oname,sampleMap(get(mae))) # Assign the object to the spesific name
		save(list=c(oname), file=paste0(path,"/",fname))

		# Add info to metadata
		descrp <- paste(oname,"MAE-object sampleMap of",study,"cohort in curatedPCaData package",sep=" ")

		metadata_raw <- rbind(metadata_raw,
			c(Title=oname,
			Description=descrp,
			BiocVersion=biocv,
			Species=species,
			RDataClass=class(get(oname))[1],
			DispatchClass="Rds",
			ResourceName=fname,
			RDataPath=paste0("curatedPCaData/",fname),Tags="sampleMap")
		)

	}
}