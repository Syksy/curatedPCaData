# curatedPCaData Version 1.0.0 _2023-01-18_

* Updated content and style for first full version 1.0
* Cleaned up deposited code used for reproduciblity
* Latest roxygenized documentation and knitr README.Rmd

# curatedPCaData Version 0.9.44 _2023-01-17_

* Fixed character string being introduced to gleason_grade in Taylor clinical metadata

# curatedPCaData Version 0.9.43 _2023-01-09_

* Clarifications to the descriptions provided in the template_prad
* limma used for DE gene analysis in trends vignette
* Mutation data barcodes added to 'mut' slots

# curatedPCaData Version 0.9.42 _2022-12-29_

* Endothelial vignette revisions

# curatedPCaData Version 0.9.41 _2022-12-19_

* Metadata sweeping vignette
* Endothelial scores trends example in EPIC, MCP-counter and xCell

# curatedPCaData Version 0.9.40 _2022-12-14_

* Endothelial cell score application vignette

# curatedPCaData Version 0.9.39 _2022-12-13_

* Fixed few special cases in Gleasons for Kim and Ren
* Additional table(s) in Overview vignette

# curatedPCaData Version 0.9.38 _2022-12-12_

* Old versioning
* Removal of metadata specific to Bioc
* Updates to low QC gene expression sample filtering in TCGA
* Updates to corresponding downstream variables (such as immunedeconv) in TCGA
* Updated Wallace et al. CIBERSORTx

# curatedPCaData Version 0.99.0 _2022-10-06_

* All content revised to conform to Bioconductor package policies

# curatedPCaData Version 0.9.33 _2022-09-27_
  
* Revised vignettes
* Revised parameters for running immunedeconv (esp. quanTIseq)
* Adding ESTIMATE from immunedeconv >= v2.1.0
* Subset to colData sample_name-column when adding new MAE slots (i.e. new assay/platform types)

# curatedPCaData Version 0.9.32 _2022-09-15_

* Clinical data revised for TCGA (xenabrowser-version) with a larger N subset appended by cBio's additional fields when available
* Brief edits to other datasets' clinical metadata as well
* Revised vignettes and new material

# curatedPCaData Version 0.9.31 _2022-09-08_

* CNA GISTIC for Hieronymus et al. replacing the old log-ratios MAE slot
* Fixes in landscapes.Rmd

# curatedPCaData Version 0.9.30 _2022-09-06_

* CNA landscape plots (original code thanks to Mike)

# curatedPCaData Version 0.9.29 _2022-09-05_

* cBioPortal processing pipeline revised
* MAEs and corresponding raw data revised
* CIBERSORTx rerun

# curatedPCaData Version 0.9.28 _2022-07-25_

* Reprocessed data

# curatedPCaData Version 0.9.27 _2022-07-18_

* Lifted computational burden from a vignette to pre-computed workspace
* Fixed class check comparison WARNINGS from BiocCheck with methods::is 

# curatedPCaData Version 0.9.26 _2022-07-15_

* CNA data matrices filtered down to protein coding (i.e. exclude miRNAs and such)

# curatedPCaData Version 0.9.25 _2022-07-14_

* biocViews / BiocType expanded, various BiocCheck fixes
* Corrected internal function/data calls
* Various misc fixes so that R CMD check is clean (except size NOTE)
* Testing Bioconductor instructions for LazyData set to false, replaced by use of 'loadPCa()'-function by the user

# curatedPCaData Version 0.9.24 _2022-07-13_

* R CMD check feedback fixes
* BiocCheck feedback fixes
* biocViews in DESCRIPTION (ExperimentData)

# curatedPCaData Version 0.9.23 _2022-07-08_

* Median shift bug in risk score calculations fixed
* Fixing R CMD check Notes

# curatedPCaData Version 0.9.21 _2022-06-23_

* Updated NAMESPACE and Rds

# curatedPCaData Version 0.9.20 _2022-06-22_

* LazyDataCompression added
* Overhaul of dependencies

# curatedPCaData Version 0.9.19 _2022-06-21_

* Fixes to vignettes
* Additional allowed values in template_prad's 'grade_group' for cases when info such as 4+3 or 3+4 is not available but instead just 7 was available in original source
* Additional necessary content, such as this NEWS file and inst/CITATION

# curatedPCaData Version 0.9.19 _2022-06-17_

* Moving through checklists for BioconductoR compatible experiment data package structure
