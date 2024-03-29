% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{curatedPCaDatasets_barbieri}
\alias{curatedPCaDatasets_barbieri}
\title{Barbieri et al. MAE-object}
\format{
A MAE object spanning prostate adenocarcinomas from Barbieri et. al
\describe{
    \item{cna.gistic}{matrix with 20844 rows and 109 columns, from GISTIC 
        discretized copy number alteration calls}
    \item{gex.relz}{matrix with 17917 rows and 31 columns, from z-score 
        normalized expression in relative to paired normals.}
    \item{mut}{RaggedExperiment with 5737 rows and 112 columns, mutation 
        data from cbioportal}
    \item{cibersort}{matrix with 22 rows and 31 columns, of cibersort based 
        deconvolution data}
    \item{xcell}{matrix with 39 rows and 31 columns, of xcell based 
        deconvolution data}
    \item{epic}{matrix with 8 rows and 31 columns, of epic based 
        deconvolution data}
    \item{quantiseq}{matrix with 11 rows and 31 columns, of quantiseq based 
        deconvolution data}
    \item{estimate}{data.frame with 4 rows and 31 columns, of cell types 
        based on ESTIMATE method}
    \item{scores}{matrix with 4 rows and 31 columns, of risk scores and AR 
        scores}
    \item{mcp}{matrix with 11 rows and 31 columns, of mcp-counter based 
        deconvolution data}
}
}
\source{
\url{https://www.cbioportal.org/study/summary?id=prad_broad}
}
\value{
A MultiAssayExperiment corresponding to the study and its available 
omics.
}
\description{
MultiAssayExperiment object containing gene expression (gex), copy
number alteration (cna), mutations (mut) and immune cell estimates for 
Barbieri et al.
}
\details{
The data comprises of primary localised prostate tumors from two 
cohorts, the Weill Cornell Medical College (WCMC; New York, NY), and the 
Uropath (Perth, Australia), which commercially provides banked urological 
tissues. None of the samples comes from patients who had received prior 
treatment for prostate cancer.
}
\examples{
mae_barbieri <- getPCa('barbieri')
}
\references{
(\href{https://pubmed.ncbi.nlm.nih.gov/22610119/}{PubMed})
}
\keyword{datasets}
