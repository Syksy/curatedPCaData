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

  
# Fetch gene names
ensembl <- biomaRt::useEnsembl(biomart = "ensembl", 
                               dataset = "hsapiens_gene_ensembl")
genes <- biomaRt::getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id',
                                       'hgnc_symbol','chromosome_name',
                                       'start_position','end_position'),
                        mart = ensembl)

# Using hgnc gene symbols 
genenames <- unique(genes$hgnc_symbol)
  
# Alternative approach using 'GenomicFeatures'
# need to have RMariaDB installed but no loaded to work because ¯\_(ツ)_/¯
# also we currently never use anything other than genes from biomart as 
# the function is currently written 
hg38_genes <- GenomicFeatures::makeTxDbFromUCSC(genome="hg38", table="refGene")
refseq_genes <- GenomicFeatures::genes(hg38_genes)

# Return 
tcga_gene_names <- list(hgnc = genenames, hg38_genes = hg38_genes, refseq_genes = refseq_genes)

# save data for internal use 
usethis::use_data(tcga_gene_names, internal = TRUE)
