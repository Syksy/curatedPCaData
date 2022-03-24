###
#
# Various utility functions
#
###

#' Update gene annotations and resolve ambiguity
#' 
#' A function that helps update gene annotations in ambiguous cases
#'
#' @noRd
#' @keywords internal
updateAnno <- function(
	# Matrix of gene expression, with rownames corresponding to mappable entities
	x,
	# Main column to extract, by default hugo gene symbols
	main = "hgnc_symbol",
	# Legitimate mapping types; Aliaes broken down into lists with ';' delimiter, others have non-unique multirow mappins between each other
	type = c("Aliases", "ensembl_gene_id", "ensembl_transcript_id", "refseq_mrna"),
	# Collapsing functions to merge rows for which the end result is multirow mapping (i.e. duplicated probes etc)
	collapse_fun,
	# Whether empty or NA rows should be omitted
	omitNAempty = TRUE,
	# Additional parameters
	...
){
	genes <- curatedPCaData:::curatedPCaData_genes
	# First argument
	type <- type[1]
	# Row splitting known aliases
	if(type == "Aliases"){
		gs <- lapply(genes[,type], FUN=function(g){
			stringr::str_to_upper(strsplit(g, ";")[[1]])
		})
		names(gs) <- genes[,main]
		# Take each hugo symbol instance only once
		gs <- gs[unique(names(gs))]
		genenames <- lapply(rownames(x), FUN=function(g){
			# TRUE/FALSE whether one or more of the aliases matched the name
			names(gs)[which(unlist(lapply(gs, FUN=function(q){
				any(q %in% stringr::str_to_upper(g))
			})))]
		})
		indices <- unlist(lapply(1:nrow(x), FUN=function(i){
			rep(rownames(x)[i], times=length(genenames[[i]]))
		}))
		# Extract mapped rows
		x <- x[indices,]
		# Replace old names with the new mapped aliases
		rownames(x) <- unlist(genenames)
		
	# Other synonym/annotation systems
	}else if(type=="ensembl_gene_id"){
		
	}else if(type=="ensembl_transcript_id"){
	
	}else if(type=="refseq_mrna"){
	
	}else{
		stop(paste("Unknown 'type':", type))
	}
	
	# If collapse function has been defined (i.e. non-missing) collapse duplicates
	if(!missing(collapse_fun)){
		x <- do.call("rbind", by(as.matrix(x), INDICES = rownames(x), FUN = collapse_fun))
	}
	
	# Omit "" or NA rows
	if(omitNAempty){
		# Include those that are not NA or equal to ""
		x <- x[which(!(is.na(rownames(x)) | rownames(x)=="" | rownames(x)=="NA")),]
	}
	
	x
}

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @rdname pipe
#' @keywords internal
#' @noRd
#' @importFrom magrittr %>%
NULL