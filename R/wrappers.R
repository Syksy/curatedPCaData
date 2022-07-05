#' Wrapper function to summarize mutations from a Raggedexperiment object
#' 
#' @param ragexp Raggedexperiment object from MAE
#' @param field Which column to collapse from inside the RaggedExperiment object
#' @return A summarized matrix with mutation information for different genes across all samples
#'
#' @examples
#' mut_abida<-curatedPCaData:::wrapper_raggedexp(mae_abida[["mut"]])
#' 
#' @noRd
#' @keywords internal
wrapper_raggedexp<-function(ragexp, field="Variant_Classification"){
I=RaggedExperiment::sparseAssay(ragexp,field)
I=as.data.frame(I)
I$gene=rownames(I)
#I=I[,c(44,1:43)]
I$gene=gsub("\\..*","",I$gene)


v=I %>% 
  dplyr::group_by(gene) %>% 
  dplyr::summarise_all(toString)

v=as.data.frame(v)

b <- data.frame(lapply(v, function(x) {
  gsub("NA, NA, NA", "NA", x)
}))
D<-data.frame(lapply(b, function(x) {
  gsub("NA, NA", "NA", x)
}))
E<-data.frame(lapply(D, function(x) {
  gsub("NA,", "", x)
}))
G<-data.frame(lapply(E, function(x) {
  gsub(", NA", "", x)
}))
# h<-data.frame(lapply(G, function(x) {
#   gsub("", "NA",x)
# }))

rownames(G)=G$gene
G <- G[ , ! names(G) %in% "gene"]

j=as.data.frame(apply(G,2,function(x)gsub('\\s+', '',x)))
j=as.data.frame(j)
j$gene<-rownames(j)

h<-data.frame(lapply(j, function(x) {
  gsub(",NA", "", x)
}))

i<-data.frame(lapply(h, function(x) {
  gsub(";", ",", x)
}))

rownames(i)<-i$gene

i<-i[,-321]

return(i)
}

#' Wrapper function to help produce oncoprint-friendly output from MAEs
#'
#' @param mae MultiAssayExperiment-object which should be collapsed
#' @param genes If only a subset of genes should be processed, a list for their names can be provided here
#' @param omics Which omics are collapsed together; eligible: "mut" and "cna.gistic"
#' @param join Type of matrix collapsing; "either" will include genes and samples found in either CNA or MUT, while "both" will only include genes and samples that were found in both
#' @param field Name of the field to extract from RaggedExperiments
#' @param map Mapping of values to names when collapsing; e.g. GISTIC values to words indicating copy number changes. List with element corresponding to 
#'
#' @return Collapsed oncoprint-friendly matrix
#'
#' @noRd
#' @keywords internal
wrapper_oncoprintify <- function(
	mae,
	genes,
	omics = c("mut", "cna.gistic"),
	join = "either",
	field = "Variant_Classification",
	map = list(
		"cna.gistic" = c("-2" = "Deep deletion", "-1" = "Deletion", "0" = "", "1" = "Gain", "2" = "High gain")
	)
){
	# Either a custom set of genes or the whole set of available genes
	if(missing(genes)){
		genes <- rownames(mutmat)
	}
	
	# Mutation portion
	if("mut" %in% omics){
		if(!"mut" %in% names(mae)){
			warning("No mutation data found in the MAE object")
			mutrag <- matrix(NA, nrow=0, ncol=0)
		}else{
			mutrag <- as.data.frame(RaggedExperiment::sparseAssay(mae[["mut"]], field))
			# Loop and bind same genes
			mutrag <- do.call("rbind", 
				# Loop over the uniquefied gene names and collapse them together into concatenated character strings
				by(mutrag, INDICES=gsub("\\.[1-99]", "", rownames(mutrag)), FUN=function(x){ apply(x, MARGIN=2,
					FUN=function(y){
						paste(na.omit(y), collapse=";")
					})
				})
			)
			# Subset to desired genes (if found)
			mutmat <- mutrag[which(rownames(mutrag) %in% genes),]		
		}
	}

	# CNA portion; mapping 	
	if("cna.gistic" %in% omics){
		if(!"cna.gistic" %in% names(mae)){
			warning("No GISTIC CNA data found in the MAE object")
			cnamat <- matrix(NA, nrow=0, ncol=0)
		}else{
			# Pick selected genes and map GISTIC coding to oncoprint-names
			cnamat <- apply(mae[["cna.gistic"]][which(rownames(mae[["cna.gistic"]]) %in% genes),], MARGIN=2, FUN=function(x){
				map$cna.gistic[match(x, names(map$cna.gistic))]
			})
			rownames(cnamat) <- rownames(mae[["cna.gistic"]][which(rownames(mae[["cna.gistic"]]) %in% genes),])
		}
	}
	
	# Populate joint matrices
	if(all(c("cna.gistic", "mut") %in% omics)){	
		if(join == "either"){
			cols <- unique(c(colnames(mutmat), colnames(cnamat)))
			rows <- unique(c(rownames(mutmat), rownames(cnamat)))
			oncomat <- matrix("", nrow=length(rows), ncol=length(cols))
			colnames(oncomat) <- cols
			rownames(oncomat) <- rows
			# Need to handle with care rows and columns which might not be present in both
			for(row in rows){
				for(col in cols){
					#print(paste(cnamat[cnarow,cnacol], mutmat[mutrow, mutcol], sep=";"))
					if(!row %in% rownames(mutmat) | !col %in% colnames(mutmat)){
						mutpart <- NA
					}else{
						mutpart <- mutmat[row,col]
					}
					if(!row %in% rownames(cnamat) | !col %in% colnames(cnamat)){
						cnapart <- NA
					}
					else{
						cnapart <- cnamat[row,col]
					}
					oncomat[row,col] <- gsub("NA", "", gsub("^NA$|^;NA|;NA$|^;$|^;|;$", "", paste(cnapart, mutpart, sep=";")))
				}			
			}
		}else if(join == "both"){
			# Only including intersections for rows and columns, so these can safely be called directly
			cols <- intersect(colnames(mutmat), colnames(cnamat))
			rows <- intersect(rownames(mutmat), rownames(cnamat))
			oncomat <- matrix("", nrow=length(rows), ncol=length(cols))
			colnames(oncomat) <- cols
			rownames(oncomat) <- rows
			for(row in rows){
				for(col in cols){
					oncomat[row,col] <- gsub("^;|;$", "", paste(mutmat[row,col], cnamat[row,col], sep=";"))
				}			
			}
		}else{
			stop(paste("Invalid parameter 'join':", join, "(should be 'either' or 'both')"))
		}
	}else if(all(omics == "cna.gistic")){
		oncomat <- cnamat
	}else if(all(omics == "mut")){
		oncomat <- mutmat
	}
	
	# Return the oncoprint-friendly matrix
	oncomat
}

#' A wrapper function for sweeping over whole curatedPCaData-package over all 'omics and datasets
#'
#' @param gene Hugo gene symbol to query for
#' @param aliases A boolean whether gene's aliases are searched for when querying; defaults to FALSE
#' @param exact Whether the query should be exact gene name and not regular expression; defaults to FALSE
#' @param drop Should matrices be dropped to vectors if only single hit occurs; defaults to FALSE
#'
#' @return A list of lists containing all hits for the queried gene
#'
#' @examples
#' wrapper_genesweep("HLA")
#' wrapper_genesweep("TP53", exact=TRUE)
#'
#' @noRd
#' @keywords internal
wrapper_genesweep <- function(
	gene,
	aliases = FALSE,
	exact = FALSE,
	drop = FALSE
){
	res <- list()
	# List of MAE objects
	maes <- grep("mae_", utils::data(package="curatedPCaData")$result[,"Item"], value=TRUE)
	# Query gene over omics
	for(mae in maes){
		i <- length(res)+1
		res[[i]] <- list()
		eval(parse(text=paste0("mae_obj <- curatedPCaData::", mae)))
		omics <- grep("gex|cna|mut", names(mae_obj), value=TRUE)	
		# Loop over omics and query rows accordingly
		for(j in 1:length(omics)){
			omic <- mae_obj[[omics[j]]]
			if(!exact){
				res[[i]][[length(res[[i]])+1]] <- omic[grep(gene, rownames(omic), value=TRUE),]
			}else{
				res[[i]][[length(res[[i]])+1]] <- omic[which(rownames(omic) %in% gene),,drop=drop]
			}
		}
		# Name slots for various omics
		names(res[[i]]) <- omics
	}
	# Name the outermost nested list according to the datasets
	names(res) <- maes
	# Return the resulting list of lists
	res
}

#' A wrapper function for sweeping a specific colData field over all datasets
#'
#' @param col Column name to query for (see template_prad for possible fields)
#' @param exact Whether the query should be exact column name and not regular expression; defaults to FALSE
#' @param drop Should matrices be dropped to vectors if only single hit occurs; defaults to FALSE
#'
#' @return A list of lists containing all hits for the queried gene
#'
#' @examples
#' wrapper_metasweep("gleason")
#' wrapper_metasweep("survival")
#' wrapper_metasweep("sample_type", exact = TRUE)
#'
#' @noRd
#' @keywords internal
wrapper_metasweep <- function(
	col,
	exact = FALSE,
	drop = FALSE
){
	res <- list()
	# List of MAE objects
	maes <- grep("mae_", utils::data(package="curatedPCaData")$result[,"Item"], value=TRUE)
	# Query gene over omics
	res <- lapply(maes, FUN=function(mae){
		eval(parse(text=paste0("mae_obj <- curatedPCaData::", mae)))
		colDat <- MultiAssayExperiment::colData(mae_obj)
		if(exact){
			colDat[,col,drop=drop]
		}else{
			colDat[,grep(col, colnames(colDat)),drop=drop]
		}
		
	})
	# Name the result lists
	names(res) <- maes	
	# Return the resulting data frame
	res
}

#' Unwraps a concatenated character metadata field from a MAE-object
#'
#' @param mae MultiAssayExperiment object from curatedPCaData package
#' @param col Name of the column to unwrap (e.g. 'other_sample', 'other_patient', ...)
#' @param vals Directly a vector of concatenized character values (for example from an external source); will override use of 'mae' and 'col' if provided
#' @param varsep Separator for different variables to perform strsplit on, for example "var1=123|var2=foo|var3=bar", by default '|'.
#' @param valsep Separator that separates variable name from its corresponding value, for example "myvar=123", where '=' is the default separator.
#' @param casts List of class cast functions in importance order to attempt cast values into; if an non-NA warning is detected, the next cast type is tested
#'
#' @return Return a data.frame where the values with a certain separator have been unwrapped from character strings
#' 
#' @examples
#' library(curatedPCaData)

#' unwrap(vals=colData(mae_abida)[,"other_feature"])
#' 
#' @noRd
#' @keywords internal
unwrap <- function(
	mae,
	col,
	vals,
	varsep = "|",
	valsep = "=",
	casts = c(as.numeric, as.character)
){
	if(missing(vals)){
		vals <- MultiAssayExperiment::colData(mae)[,col]
	}
	do.call("rbind", lapply(lapply(vals, FUN=function(x){
		strsplit(x, split=varsep, fixed=TRUE)[[1]]
		}), FUN=function(z){
			unlist(lapply(strsplit(z, split=valsep, fixed=TRUE), FUN=function(q) { q[2] }))
	}))
}
