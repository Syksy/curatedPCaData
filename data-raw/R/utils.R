###
#
# Various utility functions
#
###

###
#
# Internal functions
#
###

#' Update gene annotations and resolve ambiguity
#'
#' A function that helps update gene annotations in ambiguous cases
#'
#' @noRd
#' @keywords internal
updateAnno <- function( # Matrix of gene expression, with rownames corresponding to mappable entities
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
                       ...) {
  genes <- curatedPCaData_genes
  # First argument
  type <- type[1]
  # Row splitting known aliases
  if (type == "Aliases") {
    gs <- lapply(genes[, type], FUN = function(g) {
      stringr::str_to_upper(strsplit(g, ";")[[1]])
    })
    names(gs) <- genes[, main]
    # Take each hugo symbol instance only once
    gs <- gs[unique(names(gs))]
    genenames <- lapply(rownames(x), FUN = function(g) {
      # TRUE/FALSE whether one or more of the aliases matched the name
      names(gs)[which(unlist(lapply(gs, FUN = function(q) {
        any(q %in% stringr::str_to_upper(g))
      })))]
    })
    indices <- unlist(lapply(1:nrow(x), FUN = function(i) {
      rep(rownames(x)[i], times = length(genenames[[i]]))
    }))
    # Extract mapped rows
    x <- x[indices, ]
    # Replace old names with the new mapped aliases
    rownames(x) <- unlist(genenames)

    # Other synonym/annotation systems
  } else if (type == "ensembl_gene_id") {

  } else if (type == "ensembl_transcript_id") {

  } else if (type == "refseq_mrna") {

  } else {
    stop(paste("Unknown 'type':", type))
  }

  # If collapse function has been defined (i.e. non-missing) collapse duplicates
  if (!missing(collapse_fun)) {
    x <- do.call("rbind", by(as.matrix(x), INDICES = rownames(x), FUN = collapse_fun))
  }

  # Omit "" or NA rows
  if (omitNAempty) {
    # Include those that are not NA or equal to ""
    x <- x[which(!(is.na(rownames(x)) | rownames(x) == "" | rownames(x) == "NA")), ]
  }

  x
}

#' Function that filters rownames based on curatedPCaData's annotated gene list down to protein coding genes
#'
#' This function will use the list extracted from biomaRt and only leave rows which are reported as protein coding.
#' It will for example remove miRNAs, which may confound CNA-analyses.
#'
#' @param x Input data matrix
#'
#' @noRd
#' @keywords internal
filterToProteinCoding <- function(x) {
  x[which(rownames(x) %in% curatedPCaData_genes[which(curatedPCaData_genes$transcript_biotype == "protein_coding"), "hgnc_symbol"]), ]
}

#' Pre-format clinical metadata matrix based on template_prad; csv-file read in as df version
#'
#' Used in data-raw/download-clinical.R
#'
#' @noRd
#' @keywords internal
initial_curated_df <- function(
    df_rownames,
    template_name) {
  # import template - drafted by Jim & Svitlana
  template <- utils::read.csv(template_name, as.is = TRUE)
  output <- matrix(NA,
    ncol = nrow(template),
    nrow = length(df_rownames)
  )
  colnames(output) <- template$col.name
  rownames(output) <- df_rownames
  output <- data.frame(output)
  for (i in 1:ncol(output)) {
    class(output[, i]) <- template[i, "var.class"]
  }
  output$sample_name <- df_rownames
  return(output)
}

#' Pre-format clinical metadata matrix based on template_prad; package's own exported df version
#'
#' Used in data-raw/download-clinical.R
#'
#' @noRd
#' @keywords internal
initial_curated_internal <- function(df_rownames) {
  template <- curatedPCaData::template_prad
  output <- matrix(NA, ncol = nrow(template), nrow = length(df_rownames))
  colnames(output) <- template$col.name
  rownames(output) <- df_rownames
  output <- data.frame(output)
  for (i in 1:ncol(output)) {
    class(output[, i]) <- template[i, "var.class"]
  }
  output$sample_name <- df_rownames
  output
}
