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
	field = "Variant_Classification",
	map = list(
		"cna.gistic" = c("-2" = "Deep deletion", "-1" = "Shallow deletion", "0" = "", "1" = "Low gain", "2" = "High gain")
	)
){
	mutmat <- as.data.frame(RaggedExperiment::sparseAssay(mae[["mut"]], field))
	if(missing(genes)){
		genes <- rownames(mutmat)
	}
	
	
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
#' 
wrapper_sweep <- function(
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
			if(exact){
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

test1 <- wrapper_sweep("TP53")
test2 <- wrapper_sweep("TP53", exact = TRUE)

