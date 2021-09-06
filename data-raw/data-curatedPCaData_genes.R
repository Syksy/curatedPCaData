# Bioconductor;
# Warning! Old packages may cause some dependencies to fail, while their updating may fail if they're already in use for the session
# Safest is to update all Bioconductor packages before analyses / running pipelines

# BioConductor 'annotate' package for all sorts of conversions and genetic location info, etc
# Annotation for microarrays
# For Entrez <-> Hugo Gene Symbol mapping
# http://bioconductor.org/packages/release/bioc/html/annotate.html

# Genome wide annotation for Human
# For Entrez <-> Hugo Gene Symbol mapping (database)
# https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html

# See if legacy aliases are supported
# grep("alias", biomaRt::listAttributes(biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl"))[,"name"], value=TRUE)

# Example for extracting the affy specific probes as suggested in https://www.biostars.org/p/402677/
# listAttributes(mart)[grep("affy", listAttributes(mart)[,1]),]
  
# Fetch gene names for various annotations
curatedPCaData_genes <- biomaRt::getBM(
	attributes = 
		c(
			# ENSEMBL
			'ensembl_gene_id', 'ensembl_transcript_id',
			# Hugo
			'hgnc_symbol',
			# RefSeq
			'refseq_mrna',
			# Chromosomal information
			'chromosome_name','start_position','end_position',
			# Description
			'description',
			# Array specific probe identifiers needed by processed studies
			"affy_hg_u133a" 			
		),
	mart = biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
)
# Sort using chromosomes and then bp locations
curatedPCaData_genes <- curatedPCaData_genes[order(curatedPCaData_genes$chromosome_name, curatedPCaData_genes$start_position, curatedPCaData_genes$end_position),]
# Connect to annotation database and extract a list of gene name synonyms, aliases, legacy names
db.con <- org.Hs.eg.db::org.Hs.eg_dbconn()
# Query 
aliases <- DBI::dbGetQuery(db.con, "SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;")
curatedPCaData_genes[,"Aliases"] <- unlist(lapply(curatedPCaData_genes$hgnc_symbol, FUN=function(g){
	w <- which(aliases[,"symbol"] == g)
	if(length(w)>0){
		paste0(aliases[w,"alias_symbol"], collapse=";")
	}else{
		""
	}
}))
# Omit row names (wrong order indices)
rownames(curatedPCaData_genes) <- NULL
# Save gene information extraction date as an attribute 'date'
attr(curatedPCaData_genes, 'date') <- Sys.time()
# Save the data frame of various gene annotations for package's internal use
usethis::use_data(curatedPCaData_genes, internal = TRUE, overwrite = TRUE)

