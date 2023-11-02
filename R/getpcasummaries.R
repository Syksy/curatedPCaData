### Multiple summarization functions Provides user-friendly formatted
### information on available data and contents there-in

#' Create summary tables for colData metadata in curatedPCaData
#'
#' Create a variable value availability tables; NA-value in addition to provided
#' 'vals' and Other-values for debugging or anomalies
#'
#' @param maes List of MultiAssayExperiment objects to summarize
#' @param var.name Name of the metadata variable to look for in colData
#' @param vals Possible values for the metadata variable to tabulate
#' @param nas Whether NA or values other than those contained in 'vals' should 
#' be appended as extra columns
#'
#' @return A matrix describing the tabulates 'vals' and their proportions along 
#' with missingness
#'
#' @examples
#' mae_taylor <- getPCa('taylor')
#' mae_tcga <- getPCa('tcga')
#' getPCaSummaryTable(maes = list(Taylor = mae_taylor, TCGA = mae_tcga), 
#'     var.name = 'gleason_grade', vals = 5:10)
#'
#' @export getPCaSummaryTable
getPCaSummaryTable <- function(maes, var.name, vals, nas = TRUE) {
    mat <- matrix(NA, nrow = length(maes), ncol = ifelse(nas, length(vals) + 2, 
        length(vals)))
    rownames(mat) <- names(maes)
    for (i in seq_len(length(maes))) {
        # Get observed metadata values
        obs <- MultiAssayExperiment::colData(maes[[i]])[, var.name]
        # Iterate possible variable values
        for (j in seq_len(length(vals))) {
            if (!sum(obs == vals[j], na.rm = TRUE) == 0) {
                mat[i, j] <- paste0(sum(obs == vals[j], na.rm = TRUE), " (", 
                round(100 * sum(obs == vals[j], na.rm = TRUE)/length(obs), 0),
                "%)")
            } else {
                mat[i, j] <- "-"
            }
        }

    }
    # If Other and N/A values should be reported
    if (nas) {
        colnames(mat) <- c(vals, "Other", "N/A")
        # Other values
        vals <- c(vals, NA, NA_real_, NA_character_)
        mat[i, ncol(mat) - 1] <- paste0(sum(!obs %in% vals, na.rm = TRUE), " (",
            round(100 * sum(!obs %in% vals, na.rm = TRUE)/length(obs), 0), "%)")
        # NA-values
        mat[i, ncol(mat)] <- paste0(sum(is.na(obs), na.rm = TRUE), " (", 
            round(100 * sum(is.na(obs), na.rm = TRUE)/length(obs), 0), "%)")
    } else {
        colnames(mat) <- vals
    }
    mat
}

#' Create survival summary for metadata in curatedPCaData
#'
#' Create summary for survival-like columns (0/1 event with follow-up time) for 
#' colData metadata in curatedPCaData
#'
#' @param maes List of MultiAssayExperiment objects to summarize
#' @param event.name Name for the colData column for the end-point event; 
#' typically one of 'days_to_disease_specific_recurrence' or 
#' 'days_to_overall_survival'
#' @param time.name Name for the colData column for the time for the end-point 
#' until event or censoring; typically one of 
#' 'disease_specific_recurrence_status' or 'overall_survival_status'
#'
#' @return A matrix describing the event and follow-up information as character 
#' strings
#'
#' @examples
#' mae_taylor <- getPCa('taylor')
#' mae_tcga <- getPCa('tcga')
#' getPCaSummarySurv(maes = list(Taylor = mae_taylor, TCGA = mae_tcga), 
#'     event.name = 'disease_specific_recurrence_status', time.name = 
#'     'days_to_disease_specific_recurrence')
#'
#' @importFrom stats quantile
#' @export getPCaSummarySurv
getPCaSummarySurv <- function(maes, event.name, time.name) {
    mat <- matrix(NA, nrow = length(maes), ncol = 5)
    rownames(mat) <- names(maes)
    colnames(mat) <- c("0 (no event)", "1 (event)", "N/A (event)", 
        "Time (days, quantiles)", "N/A (time)")
    for (i in seq_len(length(maes))) {
        times <- MultiAssayExperiment::colData(maes[[i]])[, time.name]
        events <- MultiAssayExperiment::colData(maes[[i]])[, event.name]
        ifelse(all(is.na(events)), mat[i, 1] <- "-", mat[i, 1] <- paste0(sum(
            events == 0, na.rm = TRUE), " (", round(100 * sum(events == 0, 
            na.rm = TRUE)/length(events), 0), "%)"))
        ifelse(all(is.na(events)), mat[i, 2] <- "-", mat[i, 2] <- paste0(sum(
            events == 1, na.rm = TRUE), " (", round(100 * sum(events == 1, 
            na.rm = TRUE)/length(events), 0), "%)"))
        mat[i, 3] <- paste0(sum(is.na(events)), " (", round(100 * sum(
            is.na(events))/length(events), 0), "%)")
        qs <- round(quantile(times, na.rm = TRUE), 0)
        ifelse(all(is.na(qs)), mat[i, 4] <- "-", mat[i, 4] <- paste0("[", qs[1],
            ",", qs[2], ",", qs[3], ",", qs[4], ",", qs[5], "]"))
        mat[i, 5] <- paste0(sum(is.na(times)), " (", round(100 * sum(
            is.na(times))/length(times), 0), "%)")
    }
    mat
}

#' Create a summary table for sample counts for different omics across the 
#' datasets in curatedPCaData
#'
#' Create two types of summary matrices; first one takes all unique assay names,
#' and summarizes sample counts available for all provided studies. The second 
#' matrix takes gene expression, copy number alteration, and somatic mutation 
#' data, and displays sample counts for each omic individually, as well as all 
#' combinations of intersections for samples available for multi-omics.
#'
#' @param maes List of MultiAssayExperiment objects to summarize
#'
#' @return A list with two elements: first is a matrix depicting sample counts 
#' in each unique assay name for the studies; 
#' second one is a matrix with sample count overlap between the three main omics
#' (GEX = gene expression; CNA = copy number alteration; MUT = somatic 
#' mutations).
#'
#' @examples
#' mae_taylor <- getPCa('taylor')
#' mae_tcga <- getPCa('tcga')
#' getPCaSummarySamples(maes = list(Taylor = mae_taylor, TCGA = mae_tcga))
#'
#' @export getPCaSummarySamples
getPCaSummarySamples <- function(maes) {
    # Extract up-to-date assaynames over all MAE-objects
    assaynames <- c()
    for (mae in maes) {
        assaynames <- c(assaynames, names(mae))
    }
    # Unique slot names
    assaynames <- unique(assaynames)
    # Create n count or absent (NA) matrix to be printed out
    assaymat <- matrix("", nrow = length(assaynames), ncol = length(maes))
    # Create a list of lists for sample names in different omics
    sampnames <- list()
    rownames(assaymat) <- sort(assaynames)  # Alphabetic ordering
    colnames(assaymat) <- names(maes)
    for (col in seq_len(length(maes))) {
        sampnames[[length(sampnames) + 1]] <- list()
        mae <- maes[[col]]
        for (row in seq_len(nrow(assaymat))) {
            if (rownames(assaymat)[row] %in% names(mae)) {
                assaymat[row, col] <- 
                    length(colnames(mae[[rownames(assaymat)[row]]]))
                sampnames[[length(sampnames)]][[length(sampnames[[
                    length(sampnames)]]) + 1]] <- 
                        colnames(mae[[rownames(assaymat)[row]]])
                names(sampnames[[length(sampnames)]])[length(sampnames[[
                    length(sampnames)]])] <- rownames(assaymat)[row]
            }
        }
    }
    # Give MAE names
    names(sampnames) <- names(maes)
    # Transpose for more user-friendly printing in rmarkdown
    assaymat <- t(assaymat)
    # Which MAEs had 2 or more omics
    multiomics <- names(which(unlist(lapply(sampnames, FUN = function(x) {
        length(names(x)) > 1
    }))))
    # Create an overlap N count matrix for GEX & CNA, CNA & MUT, GEX & CNA &
    # MUT based on sample names
    overmat <- matrix("", ncol = 7, nrow = length(multiomics))
    colnames(overmat) <- c("GEX", "CNA", "MUT", "GEX & CNA", "GEX & MUT", 
        "CNA & MUT", "GEX & CNA & MUT")
    rownames(overmat) <- multiomics
    for (mult in multiomics) {
        samps <- sampnames[[mult]]
        gex <- samps[[grep("gex", names(samps))[1]]]
        cna <- samps[[grep("cna", names(samps))[1]]]
        mut <- samps[[grep("mut", names(samps))[1]]]
        gex_cna <- length(intersect(gex, cna))
        gex_mut <- length(intersect(gex, mut))
        cna_mut <- length(intersect(cna, mut))
        gex_cna_mut <- length(intersect(intersect(gex, cna), mut))
        overmat[mult, ] <- c(length(gex), length(cna), length(mut), gex_cna, 
            gex_mut, cna_mut, gex_cna_mut)
    }

    # Return the sample N count matrix and the overlap sample matrix
    list(Samples = assaymat, Overlap = overmat)
}

#' Create a summary table of key characteristics for each study
#'
#' @param maes A list of MAE objects or a character vector with study short ids
#' 
#' @return
#' A table with study names, references, sample counts, data source, and 
#' other notes.
#'
#' @examples
#' studies <- getPCaStudies(uniqs = TRUE)
#' tab <- getPCaSummaryStudies(studies)
#'
#' @export getPCaSummaryStudies
getPCaSummaryStudies <- function(maes) {
    if (inherits(maes, "list")) {
        # Check that all list objects are of class MultiAssayExperiment
        if (!all(unlist(
            lapply(maes, FUN=\(x){ inherits(x, "MultiAssayExperiment") }
        )))){
            stop("List should consist only of 
                MAE-objects from the curatedPCaData-package")
        }
        # Check that list naming is ok and format it if necessary
        if (is.null(names(maes)) | length(names(maes)) < length(maes)){
            # If naming vector is not provided, extract it from MAE-objects
            names(maes) <- 
                base::tolower(unlist(lapply(maes, FUN=\(x){ 
                    # Extract study names from the first word in study_name
                    strsplit(
                        MultiAssayExperiment::colData(x)$study_name, 
                    ",| ")[[1]][1] 
                })))
        }
    # Character vector of study names, will use 'getPCa' to fetch MAEs
    } else if (inherits(maes, "character")) {
        maenames <- maes
        maes <- lapply(maenames, FUN = \(id) {
            curatedPCaData::getPCa(id)
        })
        names(maes) <- maenames
    } else {
        stop("Parameter 'maes' should either be a study short id character 
            vector or preferably a named list of MAE objects")
    }
    studies <- matrix("", nrow = length(maes), ncol = 6)
    colnames(studies) <- c("Study short id", "Sample types", 
        "GEX/CNA/MUT platform(s)", "Notes", "Data source", "Reference(s)")
    studies[, 1] <- names(maes)
    # Reformat names with special cases for abbreviations vs. names
    for (index in seq_len(length(maes))) {
        mae <- names(maes)[index]
        samptypes <- table(MultiAssayExperiment::colData(
            maes[[index]])$sample_type)
        sampnames <- names(samptypes)
        samplecollapse <- paste(paste(sampnames, samptypes, sep = ": "), 
            collapse = ", ")
        studies[index, "Sample types"] <- samplecollapse
        # Annotate additional useful information and append to correct
        # positions
        if (mae %in% c("abida", "mae_abida")) {
            studies[index, "Data source"] <- "cBioPortal"
            studies[index, "Reference(s)"] <- "Abida et al."
        } else if (mae %in% c("baca", "mae_baca")) {
            studies[index, "Data source"] <- "cBioPortal"
            studies[index, "Reference(s)"] <- "Baca et al."
        } else if (mae %in% c("barbieri", "mae_barbieri")) {
            studies[index, "Data source"] <- "cBioPortal"
            studies[index, "Reference(s)"] <- "Barbieri et al."
        } else if (mae %in% c("barwick", "mae_barwick")) {
            studies[index, "Data source"] <- "GEO"
            studies[index, "GEX/CNA/MUT platform(s)"] <- "Custom DASL"
            studies[index, "Reference(s)"] <- "Barwick et al."
        } else if (mae %in% c("chandran", "mae_chandran", "yu")) {
            studies[index, "Data source"] <- "GEO"
            studies[index, "GEX/CNA/MUT platform(s)"] <- "GPL8300 [HG_U95Av2]"
            studies[index, "Reference(s)"] <- "Chandran et al., Yu et al."
        } else if (mae %in% c("friedrich", "mae_friedrich")) {
            studies[index, "Data source"] <- "GEO"
            studies[index, "GEX/CNA/MUT platform(s)"] <- "Custom Agilent array"
            studies[index, "Reference(s)"] <- "Friedrich et al."
        } else if (mae %in% c("hieronymus", "mae_hieronymus")) {
            studies[index, "Data source"] <- "GEO"
            studies[index, "GEX/CNA/MUT platform(s)"] <- 
                "GPL8737 Agilent-021529 Human CGH"
            studies[index, "Notes"] <- "CNA only"
            studies[index, "Reference(s)"] <- "Hieronymus et al."
        } else if (mae %in% c("icgcca", "mae_icgcca", "icgc-ca", "zhang")) {
            studies[index, "Data source"] <- "ICGC Data Portal (PRAD-CA)"
            studies[index, "Notes"] <- 
                "Canadian data from International Cancer Genome Collaboratory"
            studies[index, "Reference(s)"] <- "PRAD-CA in Zhang et al."
        } else if (mae %in% c("igc", "mae_igc")) {
            studies[index, "Data source"] <- "GEO"
            studies[index, "GEX/CNA/MUT platform(s)"] <- 
                "GPL570 [HG-U133_Plus_2]"
            studies[index, "Reference(s)"] <- "GEO accession code GSE2109"
        } else if (mae %in% c("kim", "mae_kim")) {
            studies[index, "Data source"] <- "GEO"
            studies[index, "GEX/CNA/MUT platform(s)"] <- "GPL5188 [HuEx-1_0-st]"
            studies[index, "Reference(s)"] <- "Kim et al."
        } else if (mae %in% c("kunderfranco", "mae_kunderfranco", 
            "peraldo-neia", "longoni")) {
            studies[index, "Data source"] <- "GEO"
            studies[index, "GEX/CNA/MUT platform(s)"] <- 
                "GPL887 Agilent-012097 Human 1A Microarray (V2)"
            studies[index, "Reference(s)"] <- 
                "Kunderfranco et al., Peraldo-Neia et al., Longoni et al."
        } else if (mae %in% c("ren", "mae_ren")) {
            studies[index, "Data source"] <- "cBioPortal"
            studies[index, "Reference(s)"] <- "Ren et al."
        } else if (mae %in% c("sun", "mae_sun")) {
            studies[index, "Data source"] <- "GEO"
            studies[index, "GEX/CNA/MUT platform(s)"] <- "GPL96 [HG-U133A]"
            studies[index, "Reference(s)"] <- "Sun et al."
        } else if (mae %in% c("taylor", "mae_taylor", "mskcc")) {
            studies[index, "Data source"] <- "GEO"
            studies[index, "Notes"] <- "Also known as MSKCC"
            studies[index, "GEX/CNA/MUT platform(s)"] <- 
                "GEX: GPL5188 [HuEx-1_0-st], CNA: GPL4091 Agilent CGH"
            studies[index, "Reference(s)"] <- "Taylor et al."
        } else if (mae %in% c("tcga", "mae_tcga")) {
            studies[index, "Data source"] <- "Xenabrowser"
            studies[index, "Reference(s)"] <- 
                "Cancer Genome Atlas Research Network, Goldman et al."
        } else if (mae %in% c("true", "mae_true")) {
            studies[index, "Data source"] <- "GEO"
            studies[index, "GEX/CNA/MUT platform(s)"] <- 
                "GPL3834 FHCRC Human Prostate PEDB cDNA v3 / v4"
            studies[index, "Reference(s)"] <- "True et al."
        } else if (mae %in% c("wallace", "mae_wallace")) {
            studies[index, "Data source"] <- "GEO"
            studies[index, "GEX/CNA/MUT platform(s)"] <- "GPL571 [HG-U133A_2]"
            studies[index, "Reference(s)"] <- "Wallace et al."
        } else if (mae %in% c("wang", "mae_wang")) {
            studies[index, "Data source"] <- "GEO"
            studies[index, "GEX/CNA/MUT platform(s)"] <- "GPL96 [HG-U133A]"
            studies[index, "Reference(s)"] <- "Wang et al., Jia et al."
        } else if (mae %in% c("weiner", "mae_weiner")) {
            studies[index, "Data source"] <- "GEO"
            studies[index, "GEX/CNA/MUT platform(s)"] <- "GPL5175 [HuEx-1_0-st]"
            studies[index, "Reference(s)"] <- "Weiner et al."
        }
    }
    studies
}

#' Create a vector of unique study identifiers available in curatedPCaData
#'
#' This function accesses the metadata available together with the 
#' curatedPCaData package and creates a vector of unique study identifier names 
#' available in the ExperimentHub resources. This vector is useful, for example,
#' when applying the 'getPCa' to create MAE from the studies.
#'
#' @param uniqs Should only unique instances be returned; by default TRUE
#'
#' @return A character vector with short study identifiers, e.g. 'abida', 
#' 'baca', ...
#'
#' @details Uses the internal 'metadata.csv' file, which is also used to provide
#' the metadata
#' available in ExperimentHub. 
#'
#' @examples
#' studyids <- getPCaStudies()
#'
#' @export getPCaStudies
getPCaStudies <- function(uniqs = TRUE) {
    assays_file <- system.file("extdata", "metadata.csv",
        package = "curatedPCaData", mustWork = TRUE)
    assay_metadat <- read.csv(assays_file, stringsAsFactors = FALSE)
    titles <- unlist(lapply(assay_metadat$Title, FUN = \(x) {
        strsplit(x, "_")[[1]][1]
    }))
    if (uniqs) {
        unique(titles)
    } else {
        titles
    }
}
