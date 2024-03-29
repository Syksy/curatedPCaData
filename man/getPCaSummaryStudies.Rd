% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getpcasummaries.R
\name{getPCaSummaryStudies}
\alias{getPCaSummaryStudies}
\title{Create a summary table of key characteristics for each study}
\usage{
getPCaSummaryStudies(maes)
}
\arguments{
\item{maes}{A list of MAE objects or a character vector with study short ids}
}
\value{
A table with study names, references, sample counts, data source, and 
other notes.
}
\description{
Create a summary table of key characteristics for each study
}
\examples{
studies <- getPCaStudies(uniqs = TRUE)
tab <- getPCaSummaryStudies(studies)

}
