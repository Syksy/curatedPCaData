###
#
# Variable functions helping exporting of data/information from curatedPCaData
#
###

#' Export ExperimentHub-friendly metadata csv prototype
#'
#' This function greps over MAE-objects available via curatedPCaData and starts constructing ExperimentHub-compatible metadata.csv, first as a data.frame
#' 
#' @examples
#' metadat <- export_metadata(timestamp = "20230215", export = FALSE)
#'
#'
#' @noRd
#' @keywords internal
export_metadata <- function(
	timestamp, # Time stamp to append in file names
	path = getwd(), # Location on drive to save the files into
	export = TRUE, # Whether the actual .Rds objects should be exported; if FALSE, the relevant metadata is created as if the Rds hada been created
	...
){
	# If user does not provide a time stamp, construct one from the system time ('%Y%m%d'-format)
	if(missing(timestamp)){
		timestamp <- gsub("-", "", substr(Sys.time(), 0, 10))
	}

	# List MAEs available in curatedPCaData
	# Get available maes in curatedPCaData package and specify a path to the folder where to save
	maes <- grep("mae", utils::data(package="curatedPCaData")$result, value=TRUE)
	# Load MAEs
	data(list = maes, package = "curatedPCaData")
	
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
	metadata <- data.frame()
	#	Title = character(0),	
	#	Description = character(0),	
	#	BiocVersion = character(0),	
	#	Genome = character(0),	
	#	SourceType = character(0),	
	#	SourceUrl = character(0),
	#	SourceVersion = character(0),	
	#	Species = character(0),
	#	TaxonomyId = character(0),
	#	Coordinate_1_based = character(0),
	#	DataProvider = character(0),
	#	Maintainer = character(0),
	#	RDataClass = character(0),
	#	DispatchClass = character(0),
	#	ResourceName = character(0),
	#	RDataPath = character(0),
	#	Tags = character(0)
	#)
	for(i in 1:length(maes)){  #length(maes)
		mae <- maes[i]  # Get specific mae & separate the study name for future naming
		study <- gsub("mae_","",mae)
		ex <- names(MultiAssayExperiment::experiments(get(mae))) # Get the names of object ("experiments") in mae

		# Go through all "experiments" and save the object into its own .RDs file
		for(x in ex){
			oname <- paste0(study,"_",x,"_",timestamp)  # object name
			fname <- paste0(oname,".Rds")  # file name

			assign(oname,get(mae)[[x]]) # Assign the object to the specific name
			if(export) save(list=c(oname), file=paste0(path,"/",fname))

			# Add info to metadata
			descrp <- paste(oname,descriptions[x],"data of",study,"cohort in curatedPCaData package",sep=" ")

			metadata <- rbind(metadata, 
				c(
					Title = oname,
					Description = descrp,
					BiocVersion = biocv,
					Genome = "",
					SourceType = "",
					SourceURL = "",
					SourceVersion = "",
					Species = species,
					TaxonomyId = "",
					Coordinate_1_based = "",
					DataProvider = "",
					Maintainer = "",
					RDataClass = class(get(oname))[1],
					DispatchClass = "Rds",
					ResourceName = fname,
					RDataPath = paste0("curatedPCaData/", fname),
					Tags = x
				)
			)
		}

		# Save clinical info
		oname <- paste0(study,"_colData_",timestamp)  # object name
		fname <- paste0(oname,".Rds")  # file name

		assign(oname,MultiAssayExperiment::colData(get(mae))) # Assign the object to the spesific name
		if(export) save(list=c(oname), file=paste0(path,"/",fname))

		# Add info to metadata
		descrp <- paste(oname,"Clinical metadata (colData-slot) of",study,"cohort in curatedPCaData package",sep=" ")

		metadata <- rbind(metadata,
			c(
				Title=oname,
				Description=descrp,
				BiocVersion=biocv,
				Genome = "",
				SourceType = "",
				SourceURL = "",
				SourceVersion = "",
				Species=species,
				TaxonomyId = "",
				Coordinate_1_based = "",
				DataProvider = "",
				Maintainer = "",
				RDataClass=class(get(oname))[1],
				DispatchClass="Rds",
				ResourceName=fname,
				RDataPath=paste0("curatedPCaData/",fname),
				Tags="clinical"
			)
		)

		# Save sampleMap
		oname <- paste0(study,"_sampleMap_",timestamp)  # object name
		fname <- paste0(oname,".Rds")  # file name

		assign(oname,MultiAssayExperiment::sampleMap(get(mae))) # Assign the object to the spesific name
		if(export) save(list=c(oname), file=paste0(path,"/",fname))

		# Add info to metadata
		descrp <- paste(oname,"MAE-object sampleMap of",study,"cohort in curatedPCaData package",sep=" ")

		metadata <- rbind(metadata,
			c(Title=oname,
			Description=descrp,
			BiocVersion=biocv,
			Genome = "",
			SourceType = "",
			SourceURL = "",
			SourceVersion = "",
			Species=species,
			TaxonomyId = "",
			Coordinate_1_based = "",
			DataProvider = "",
			Maintainer = "",
			RDataClass=class(get(oname))[1],
			DispatchClass="Rds",
			ResourceName=fname,
			RDataPath=paste0("curatedPCaData/",fname),Tags="sampleMap"
			)
		)
	}
	# Set column names properly
	colnames(metadata) <- c("Title", "Description", "BiocVersion", "Genome", "SourceType", "SourceURL", "SourceVersion", "Species", "TaxonomyId", "Coordinate_1_based", "DataProvider", "Maintainer", "RDataClass", "DispatchClass", "ResourceName", "RDataPath", "Tags")
	
	metadata
}
metadat <- export_metadata(timestamp = "20230215", export = FALSE)
metadat


