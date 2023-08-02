###
#
# Multiple summarization functions
# Provides user-friendly formatted information on available data and contents there-in
#
###

#' Create summary tables for colData metadata in curatedPCaData
#'
#' Create a variable value availability tables; NA-value in addition to provided 'vals' and Other-values for debugging or anomalies
#'
#' @param maes List of MultiAssayExperiment objects to summarize
#' @param var.name Name of the metadata variable to look for in colData
#' @param vals Possible values for the metadata variable to tabulate
#' @param nas Whether NA or values other than those contained in 'vals' should be appended as extra columns
#'
#' @return A matrix describing the tabulates 'vals' and their proportions along with missingness
#'
#' @examples
#' mae_taylor <- getPCa("taylor")
#' mae_tcga <- getPCa("tcga")
#' getPCaSummaryTable(maes = list(Taylor = mae_taylor, TCGA = mae_tcga), var.name = "gleason_grade", vals = 5:10)
#'
#' @export getPCaSummaryTable
getPCaSummaryTable <- function(maes, var.name, vals, nas = TRUE)
{
	mat <- matrix(NA, nrow = length(maes), ncol = ifelse(nas, length(vals)+2, length(vals)))
	rownames(mat) <- names(maes)
	for (i in seq_len(length(maes))) {
	  # Get observed metadata values
	  obs <- MultiAssayExperiment::colData(maes[[i]])[,var.name]
	  # Iterate possible variable values
	  for (j in seq_len(length(vals))) {
	    if (!sum(obs == vals[j], na.rm = TRUE) == 0) {
	      mat[i, j] <- paste0(sum(obs == vals[j], na.rm = TRUE), " (", round(100 * sum(obs == vals[j], na.rm = TRUE) / length(obs), 0), "%)")
	    } else {
	      mat[i, j] <- "-"
	    }
	  }

	}
      # If Other and N/A values should be reported
      if(nas){
          colnames(mat) <- c(vals, "Other", "N/A")
          # Other values
          vals <- c(vals, NA, NA_real_, NA_character_)
          mat[i, ncol(mat) - 1] <- paste0(sum(!obs %in% vals, na.rm = TRUE), " (", round(100 * sum(!obs %in% vals, na.rm = TRUE) / length(obs), 0), "%)")
          # NA-values
          mat[i, ncol(mat)] <- paste0(sum(is.na(obs), na.rm = TRUE), " (", round(100 * sum(is.na(obs), na.rm = TRUE) / length(obs), 0), "%)")
      }else{
          colnames(mat) <- vals
      }
      mat
}

#' Create survival summary for metadata in curatedPCaData
#'
#' Create summary for survival-like columns (0/1 event with follow-up time) for colData metadata in curatedPCaData
#'
#' @param maes List of MultiAssayExperiment objects to summarize
#' @param event.name Name for the colData column for the end-point event; typically one of 'days_to_disease_specific_recurrence' or 'days_to_overall_survival'
#' @param time.name Name for the colData column for the time for the end-point until event or censoring; typically one of 'disease_specific_recurrence_status' or 'overall_survival_status'
#'
#' @return A matrix describing the event and follow-up information as character strings
#'
#' @examples
#' mae_taylor <- getPCa("taylor")
#' mae_tcga <- getPCa("tcga")
#' getPCaSummarySurv(maes = list(Taylor = mae_taylor, TCGA = mae_tcga), event.name = "disease_specific_recurrence_status", time.name = "days_to_disease_specific_recurrence")
#'
#' @export getPCaSummarySurv
getPCaSummarySurv <- function(maes, event.name, time.name){
	mat <- matrix(NA, nrow = length(maes), ncol = 5)
	rownames(mat) <- names(maes)
	colnames(mat) <- c("0 (no event)", "1 (event)", "N/A (event)", "Time (days, quantiles)", "N/A (time)")
	for (i in seq_len(length(maes))) {
	  times <- MultiAssayExperiment::colData(maes[[i]])[,time.name]
	  events <- MultiAssayExperiment::colData(maes[[i]])[,event.name]
	  ifelse(all(is.na(events)),
	    mat[i, 1] <- "-",
	    mat[i, 1] <- paste0(sum(events == 0, na.rm = TRUE), " (", round(100 * sum(events == 0, na.rm = TRUE) / length(events), 0), "%)")
	  )
	  ifelse(all(is.na(events)),
	    mat[i, 2] <- "-",
	    mat[i, 2] <- paste0(sum(events == 1, na.rm = TRUE), " (", round(100 * sum(events == 1, na.rm = TRUE) / length(events), 0), "%)")
	  )
	  mat[i, 3] <- paste0(sum(is.na(events)), " (", round(100 * sum(is.na(events)) / length(events), 0), "%)")
	  qs <- round(quantile(times, na.rm = TRUE), 0)
	  ifelse(all(is.na(qs)),
	    mat[i, 4] <- "-",
	    mat[i, 4] <- paste0("[", qs[1], ",", qs[2], ",", qs[3], ",", qs[4], ",", qs[5], "]")
	  )
	  mat[i, 5] <- paste0(sum(is.na(times)), " (", round(100 * sum(is.na(times)) / length(times), 0), "%)")
	}
    mat
}
