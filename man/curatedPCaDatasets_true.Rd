% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{curatedPCaDatasets_true}
\alias{curatedPCaDatasets_true}
\title{True et al. MAE-object}
\format{
An MAE object spanning 32 men with prostate cancer
\describe{
    \item{gex.logr}{matrix with 3615 rows and 32 columns, of gene expression 
        data}
    \item{cibersort}{matrix with 22 rows and 32 columns, of cibersort based 
        deconvolution data}
    \item{epic}{matrix with 8 rows and 32 columns, of epic based 
        deconvolution data}
    \item{quantiseq}{matrix with 11 rows and 32 columns, of quantiseq based 
        deconvolution data}
    \item{mcp}{matrix with 7 rows and 32 columns, of mcp-counter based 
        deconvolution data}
    \item{scores}{matrix with 2 rows and 32 columns, of risk scores and AR 
        scores}
    \item{estimate}{data.frame with 4 rows and 32 columns, of cell types 
        based on ESTIMATE method}
}
}
\source{
\url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5132}
}
\value{
A MultiAssayExperiment corresponding to the study and its available 
omics.
}
\description{
MultiAssayExperiment object containing gene expression (gex) and immune cell 
estimates from True et al.
}
\details{
The expression data for this dataset was produced by an older two 
colour chip that provides the relative expression levels for between tumor 
and normal.  Due to the chip type, there was a lack of gene overlap with 
some of the downstream methods, thus some downstream methods are missing.
}
\examples{
mae_true <- getPCa('true')
}
\keyword{datasets}
