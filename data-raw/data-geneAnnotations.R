###
#
# Pre-package functions for updating gene annotations (based one unexported object curatedPCaData:::curatedPCaData_genes with biomaRt)
#
###

#' Update gene annotations and resolve ambiguity
updateAnno <- function(
	# Matrix of gene expression, with rownames corresponding to mappable entities
	x,
	# Main column to extract, by default hugo gene symbols
	main = "hgnc_symbol",
	# Legitimate mapping types; Aliaes broken down into lists with ';' delimiter
	type = c("Aliases", "ensembl_gene_id", "ensembl_transcript_id", "refseq_mrna"),
	# Collapsing functions to merge rows for which the end result is multirow mapping (i.e. duplicated probes etc); by default col-wise median
	collapse_fun = function(z) {apply(z, MARGIN = 2, FUN = stats::median)},
	# Whether empty or NA rows should be omitted
	omitNAempty = TRUE,
	# Additional parameters
	...
){
	genes <- curatedPCaData:::curatedPCaData_genes
	# First argument
	type <- type[1]
	# Row splitting known aliases
	if(type[1] == "Aliases"){
		gs <- lapply(genes[,type], FUN=function(g){
			stringr::str_to_upper(strsplit(g, ";")[[1]])
		})
		names(gs) <- genes[,main]
		genenames <- lapply(rownames(x), FUN=function(g){
			# TRUE/FALSE whether one or more of the aliases matched the name
			names(gs)[which(unlist(lapply(gs, FUN=function(q){
				any(q %in% stringr::str_to_upper(g))
			})))]
		})
	# Other synonym/annotation systems
	}else{
		
	}
	
	# If collapse function has been defined (i.e. non-missing) collapse duplicates
	if(!missing(collapse_fun)){
	
	}
	
	# Omit "" or NA rows
	if(omitNAempty){
	
	}
	
	genenames
}

updateAnno(curatedPCaData::mae_sun[["gex"]])

