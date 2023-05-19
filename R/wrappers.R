#' Wrapper function to summarize mutations from a Raggedexperiment object
#'
#' @param ragexp Raggedexperiment object from MAE
#' @param field Which column to collapse from inside the RaggedExperiment object
#' @return A summarized matrix with mutation information for different genes across all samples
#'
#' @examples
#' mae_abida <- getPCa("abida")
#' mut_abida <- curatedPCaData:::wrapperRaggedexp(mae_abida[["mut"]])
#'
#' @noRd
#' @keywords internal
wrapperRaggedexp <- function(ragexp, field = "Variant_Classification") {
  I <- RaggedExperiment::sparseAssay(ragexp, field)
  I[is.na(I)] <- ""
  I <- as.data.frame(I)
  I$gene <- rownames(I)
  I$gene <- gsub("\\..*", "", I$gene)

  # v=I %>%
  #  dplyr::group_by(.data$gene) %>%
  #  dplyr::summarise_all(toString)
  v <- I |>
    dplyr::group_by(.data$gene) |>
    dplyr::summarise_all(toString)

  v <- as.data.frame(v)

  rownames(v) <- v$gene
  v <- v[, -which(names(v) %in% "gene")]

  for (i in seq_len(ncol(v))) {
    v[, i] <- gsub(",", "", v[, i])
  }

  for (i in seq_len(ncol(v))) {
    v[, i] <- stringr::str_squish(v[, i])
  }

  return(v)
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
wrapperOncoprintify <- function(
    mae,
    genes,
    omics = c("mut", "cna.gistic"),
    join = "either",
    field = "Variant_Classification",
    map = list(
      "cna.gistic" = c("-2" = "Deep deletion", "-1" = "Deletion", "0" = "", "1" = "Gain", "2" = "High gain")
    )) {
  # Either a custom set of genes or the whole set of available genes
  if (missing(genes)) {
    genes <- rownames(mutmat)
  }
  # Mutation portion
  if ("mut" %in% omics) {
    if (!"mut" %in% names(mae)) {
      warning("No mutation data found in the MAE object")
      mutrag <- matrix(NA, nrow = 0, ncol = 0)
    } else {
      mutrag <- as.data.frame(RaggedExperiment::sparseAssay(mae[["mut"]], field))
      # Loop and bind same genes
      mutrag <- do.call(
        "rbind",
        # Loop over the uniquefied gene names and collapse them together into concatenated character strings
        by(mutrag, INDICES = gsub("\\.[1-99]", "", rownames(mutrag)), FUN = function(x) {
          apply(x,
            MARGIN = 2,
            FUN = function(y) {
              paste(na.omit(y), collapse = ";")
            }
          )
        })
      )
      # Subset to desired genes (if found)
      mutmat <- mutrag[which(rownames(mutrag) %in% genes), ]
    }
  }
  # CNA portion; mapping
  if ("cna.gistic" %in% omics) {
    if (!"cna.gistic" %in% names(mae)) {
      warning("No GISTIC CNA data found in the MAE object")
      cnamat <- matrix(NA, nrow = 0, ncol = 0)
    } else {
      # Pick selected genes and map GISTIC coding to oncoprint-names
      cnamat <- apply(mae[["cna.gistic"]][which(rownames(mae[["cna.gistic"]]) %in% genes), ], MARGIN = 2, FUN = function(x) {
        map$cna.gistic[match(x, names(map$cna.gistic))]
      })
      rownames(cnamat) <- rownames(mae[["cna.gistic"]][which(rownames(mae[["cna.gistic"]]) %in% genes), ])
    }
  }

  # Populate joint matrices
  if (all(c("cna.gistic", "mut") %in% omics)) {
    if (join == "either") {
      cols <- unique(c(colnames(mutmat), colnames(cnamat)))
      rows <- unique(c(rownames(mutmat), rownames(cnamat)))
      oncomat <- matrix("", nrow = length(rows), ncol = length(cols))
      colnames(oncomat) <- cols
      rownames(oncomat) <- rows
      # Need to handle with care rows and columns which might not be present in both
      for (row in rows) {
        for (col in cols) {
          # print(paste(cnamat[cnarow,cnacol], mutmat[mutrow, mutcol], sep=";"))
          if (!row %in% rownames(mutmat) | !col %in% colnames(mutmat)) {
            mutpart <- NA
          } else {
            mutpart <- mutmat[row, col]
          }
          if (!row %in% rownames(cnamat) | !col %in% colnames(cnamat)) {
            cnapart <- NA
          } else {
            cnapart <- cnamat[row, col]
          }
          oncomat[row, col] <- gsub("NA", "", gsub("^NA$|^;NA|;NA$|^;$|^;|;$", "", paste(cnapart, mutpart, sep = ";")))
        }
      }
    } else if (join == "both") {
      # Only including intersections for rows and columns, so these can safely be called directly
      cols <- intersect(colnames(mutmat), colnames(cnamat))
      rows <- intersect(rownames(mutmat), rownames(cnamat))
      oncomat <- matrix("", nrow = length(rows), ncol = length(cols))
      colnames(oncomat) <- cols
      rownames(oncomat) <- rows
      for (row in rows) {
        for (col in cols) {
          oncomat[row, col] <- gsub("^;|;$", "", paste(mutmat[row, col], cnamat[row, col], sep = ";"))
        }
      }
    } else {
      stop("Invalid parameter 'join' (should be 'either' or 'both')")
    }
  } else if (all(omics == "cna.gistic")) {
    oncomat <- cnamat
  } else if (all(omics == "mut")) {
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
#' wrapperGenesweep("HLA")
#' wrapperGenesweep("TP53", exact = TRUE)
#'
#' @noRd
#' @keywords internal
wrapperGenesweep <- function(
    gene,
    aliases = FALSE,
    exact = FALSE,
    drop = FALSE) {
  res <- list()
  # List of MAE objects
  maes <- grep("mae_", utils::data(package = "curatedPCaData")$result[, "Item"], value = TRUE)
  # As default, use mae_tcga if eval fails
  mae_obj <- getPCa("tcga")
  # Query gene over omics
  for (mae in maes) {
    i <- length(res) + 1
    res[[i]] <- list()
    eval(parse(text = paste0("mae_obj <- curatedPCaData::", mae)))
    omics <- grep("gex|cna|mut", names(mae_obj), value = TRUE)
    # Loop over omics and query rows accordingly
    for (j in seq_len(length(omics))) {
      omic <- mae_obj[[omics[j]]]
      if (!exact) {
        res[[i]][[length(res[[i]]) + 1]] <- omic[grep(gene, rownames(omic), value = TRUE), ]
      } else {
        res[[i]][[length(res[[i]]) + 1]] <- omic[which(rownames(omic) %in% gene), , drop = drop]
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

#' Wrapper function for sweeping a specific colData field over all datasets
#'
#' @param col Column name to query for (see template_prad for possible fields)
#' @param exact Whether the query should be exact column name and not regular expression; defaults to FALSE
#' @param drop Should matrices be dropped to vectors if only single hit occurs; defaults to FALSE
#'
#' @return A list of lists containing all hits for the queried gene
#'
#' @examples
#' wrapperMetasweep("gleason")
#' wrapperMetasweep("survival")
#' wrapperMetasweep("sample_type", exact = TRUE)
#'
#' @noRd
#' @keywords internal
wrapperMetasweep <- function(
    col,
    exact = FALSE,
    drop = FALSE) {
  # List of MAE objects
  maes <- grep("mae_", utils::data(package = "curatedPCaData")$result[, "Item"], value = TRUE)
  # As default, use mae_tcga if eval fails
  mae_obj <- getPCa("tcga")
  # Query gene over omics
  res <- lapply(maes, FUN = function(mae) {
    eval(parse(text = paste0("mae_obj <- curatedPCaData::", mae)))
    colDat <- MultiAssayExperiment::colData(mae_obj)
    if (exact) {
      colDat[, col, drop = drop]
    } else {
      colDat[, grep(col, colnames(colDat)), drop = drop]
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
#' mae_abida <- getPCa("abida")
#' unwrap(vals = colData(mae_abida)[, "other_feature"])
#'
#' @noRd
#' @keywords internal
unwrap <- function(
    mae,
    col,
    vals,
    varsep = "|",
    valsep = "=",
    casts = c(as.numeric, as.character)) {
  if (missing(vals)) {
    vals <- MultiAssayExperiment::colData(mae)[, col]
  }
  do.call("rbind", lapply(lapply(vals, FUN = function(x) {
    strsplit(x, split = varsep, fixed = TRUE)[[1]]
  }), FUN = function(z) {
    unlist(lapply(strsplit(z, split = valsep, fixed = TRUE), FUN = function(q) {
      q[2]
    }))
  }))
}

#' A wrapper function for sorting alterations in oncoprint
#'
#' @param alt_matrix Alteration matrix of copy-number,mutations and ERG fusions
#'
#' @return A sorted alteration matrix for oncoprints
#'
#' @examples
#' wrapperSortonco(alt_matrix)
#'
#' @noRd
#' @keywords internal
wrapperSortonco <- function(alt_matrix) {
  order <- c(
    "Amplification;Fusion", "Amplification;Missense Mutation;Fusion", "Amplification;Missense Mutation", "Amplification;Frameshift Mutation;Fusion",
    "Amplification;Frameshift Mutation", "Amplification;Splice Mutation;Fusion", "Amplification;Splice Mutation",
    "Amplification;Inframe Mutation;Fusion", "Amplification;Inframe Mutation",
    "Deep deletion;Fusion", "Deep deletion;Missense Mutation;Fusion", "Deep deletion;Missense Mutation", "Deep deletion;Frameshift Mutation;Fusion",
    "Deep deletion;Frameshift Mutation", "Deep deletion;Splice Mutation;Fusion", "Deep deletion;Splice Mutation",
    "Deep deletion;Inframe Mutation;Fusion", "Deep deletion;Inframe Mutation",
    "Gain;Fusion", "Gain;Missense Mutation;Fusion", "Gain;Missense Mutation", "Gain;Frameshift Mutation;Fusion", "Gain;Frameshift Mutation",
    "Gain;Splice Mutation;Fusion", "Gain;Splice Mutation", "Gain;Inframe Mutation;Fusion", "Gain;Inframe Mutation",
    "Shallow deletion;Fusion", "Shallow deletion;Missense Mutation;Fusion", "Shallow deletion;Missense Mutation", "Shallow deletion;Frameshift Mutation;Fusion",
    "Shallow deletion;Frameshift Mutation", "Shallow deletion;Splice Mutation", "Shallow deletion;Splice Mutation;Fusion",
    "Shallow deletion;Inframe Mutation;Fusion", "Shallow deletion;Inframe Mutation", ";Fusion", "Amplification;", "Amplification",
    "Deep deletion;", "Deep deletion", "Gain;", "Gain", "Shallow deletion;", "Shallow deletion",
    ";Missense Mutation;Fusion", ";Missense Mutation", ";Frameshift Mutation;Fusion", ";Frameshift Mutation",
    ";Splice Mutation;Fusion", ";Splice Mutation", ";Inframe Mutation;Fusion", ";Inframe Mutation", ";", ""
  )
  value <- seq_len(55)
  order_df <- data.frame(order = order, value = value)
  alt_matrix_t <- as.data.frame(t(alt_matrix))
  order1 <- alt_matrix_t[order(match(alt_matrix_t[, 1], order_df$order)), , drop = FALSE]
  df <- data.frame(matrix(ncol = ncol(order1), nrow = 0))
  colnames(df) <- colnames(order1)
  if (ncol(alt_matrix_t) >= 2) {
    for (i in seq_len(nrow(order_df))) {
      if (nrow(order1[order1[, 1] == order_df[i, 1], ]) != 0) {
        col1 <- order1[order1[, 1] == order_df[i, 1], ]
        col1 <- col1[order(match(col1[, 2], order_df$order)), ]

        if (ncol(alt_matrix_t) == 2) {
          df <- rbind(df, col1)
        }
        if (ncol(alt_matrix_t) >= 3) {
          for (j in seq_len(nrow(order_df))) {
            if (nrow(col1[col1[, 2] == order_df[j, 1], ]) != 0) {
              col2 <- col1[col1[, 2] == order_df[j, 1], ]
              col2 <- col2[order(match(col2[, 3], order_df$order)), ]
              if (ncol(alt_matrix_t) == 3) {
                df <- rbind(df, col2)
              }
              if (ncol(alt_matrix_t) >= 4) {
                for (k in seq_len(nrow(order_df))) {
                  if (nrow(col2[col2[, 3] == order_df[k, 1], ]) != 0) {
                    col3 <- col2[col2[, 3] == order_df[k, 1], ]
                    col3 <- col3[order(match(col3[, 4], order_df$order)), ]
                    if (ncol(alt_matrix_t) == 4) {
                      df <- rbind(df, col3)
                    }
                    if (ncol(alt_matrix_t) >= 5) {
                      for (l in seq_len(nrow(order_df))) {
                        if (nrow(col3[col3[, 4] == order_df[l, 1], ]) != 0) {
                          col4 <- col3[col3[, 4] == order_df[l, 1], ]
                          col4 <- col4[order(match(col4[, 5], order_df$order)), ]
                          if (ncol(alt_matrix_t) == 5) {
                            df <- rbind(df, col4)
                          }
                          if (ncol(alt_matrix_t) >= 6) {
                            for (m in seq_len(nrow(order_df))) {
                              if (nrow(col4[col4[, 5] == order_df[m, 1], ]) != 0) {
                                col5 <- col4[col4[, 5] == order_df[m, 1], ]
                                col5 <- col5[order(match(col5[, 6], order_df$order)), ]
                                if (ncol(alt_matrix_t) == 6) {
                                  df <- rbind(df, col5)
                                }
                                if (ncol(alt_matrix_t) >= 7) {
                                  for (n in seq_len(nrow(order_df))) {
                                    if (nrow(col5[col5[, 6] == order_df[n, 1], ]) != 0) {
                                      col6 <- col5[col5[, 6] == order_df[n, 1], ]
                                      col6 <- col6[order(match(col6[, 7], order_df$order)), ]
                                      if (ncol(alt_matrix_t) == 7) {
                                        df <- rbind(df, col6)
                                      }
                                      if (ncol(alt_matrix_t) >= 8) {
                                        for (o in seq_len(nrow(order_df))) {
                                          if (nrow(col6[col6[, 7] == order_df[o, 1], ]) != 0) {
                                            col7 <- col6[col6[, 7] == order_df[o, 1], ]
                                            col7 <- col7[order(match(col7[, 8], order_df$order)), ]
                                            if (ncol(alt_matrix_t) == 8) {
                                              df <- rbind(df, col7)
                                            }
                                            if (ncol(alt_matrix_t) >= 9) {
                                              for (p in seq_len(nrow(order_df))) {
                                                if (nrow(col7[col7[, 8] == order_df[p, 1], ]) != 0) {
                                                  col8 <- col7[col7[, 8] == order_df[p, 1], ]
                                                  col8 <- col8[order(match(col8[, 9], order_df$order)), ]
                                                  if (ncol(alt_matrix_t) == 9) {
                                                    df <- rbind(df, col8)
                                                  }
                                                  if (ncol(alt_matrix_t) >= 10) {
                                                    for (q in seq_len(nrow(order_df))) {
                                                      if (nrow(col8[col8[, 9] == order_df[q, 1], ]) != 0) {
                                                        col9 <- col8[col8[, 9] == order_df[q, 1], ]
                                                        col9 <- col9[order(match(col9[, 10], order_df$order)), ]
                                                        if (ncol(alt_matrix_t) == 10) {
                                                          df <- rbind(df, col9)
                                                        }
                                                      }
                                                    }
                                                  }
                                                }
                                              }
                                            }
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  } else {
    return(order1)
  }
  return(df)
}
