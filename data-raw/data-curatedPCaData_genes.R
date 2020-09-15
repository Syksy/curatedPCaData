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
			'description'
		),
	mart = biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
)
# Sort using chromosomes and then bp locations
curatedPCaData_genes <- curatedPCaData_genes[order(curatedPCaData_genes$chromosome_name, curatedPCaData_genes$start_position, curatedPCaData_genes$end_position),]
# Omit row names (wrong order indices)
rownames(curatedPCaData_genes) <- NULL
# Save gene information extraction date as an attribute 'date'
attr(curatedPCaData_genes, 'date') <- Sys.time()
# Save the data frame of various gene annotations for package's internal use
usethis::use_data(curatedPCaData_genes, internal = TRUE, overwrite = TRUE)

# Using hgnc gene symbols 
#genenames <- unique(genes$hgnc_symbol)
  
# Alternative approach using 'GenomicFeatures'
# need to have RMariaDB installed but no loaded to work because ¯\_(ツ)_/¯
# also we currently never use anything other than genes from biomart as 
# the function is currently written 
#hg38_genes <- GenomicFeatures::makeTxDbFromUCSC(genome="hg38", table="refGene")
# Current makeTxDbFromUCSC would require to retain connection to UCSC database
#refseq_genes <- GenomicFeatures::genes(hg38_genes)

# Return 
#gene_names <- list(hgnc = genenames, hg38_genes = hg38_genes, refseq_genes = refseq_genes)
#gene_names <- list(hgnc = genenames, refseq_genes = refseq_genes)

# save data for internal use 
#usethis::use_data(gene_names, internal = TRUE)

