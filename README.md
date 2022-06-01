
# curatedPCaData <img src="man/figures/hex.png" align="right" height="139" />

<!-- badges: start -->
<!-- badges: end -->

## Overview

`curatedPCaData` is a collection of publically available and annotated
data resources concerning prostate cancer.

## Installation

You can install `curatedPCaData` from GitHub with:


    # install.packages("devtools")
    devtools::install_github("Syksy/curatedPCaData")

Above may fail depending on the connection stability, as the package is
relatively large. In this case it’s best to install via the tarball. To
build the package tarball from a cloned git repo, run the following in
terminal / command prompt while in the root of the project:

    R CMD build curatedPCaData

One install the self-built tarball or download the a premade latest
package tarball from the releases page and install it with:

    R CMD INSTALL curatedPCaData_x.y.z.tar.gz

<!--- add BioConductor once up --->

## Usage

Simple example use of curated datasets and ’omics there-in:


    library(curatedPCaData)

    curatedPCaData::mae_tcga
    ## A MultiAssayExperiment object of 10 listed
    ##  experiments with user-defined names and respective classes.
    ##  Containing an ExperimentList class object of length 10:
    ##  [1] cna.gistic: matrix with 19645 rows and 322 columns
    ##  [2] gex.fpkm: matrix with 58387 rows and 369 columns
    ##  [3] mut: RaggedExperiment with 29286 rows and 320 columns
    ##  [4] cibersort: matrix with 22 rows and 369 columns
    ##  [5] xcell: matrix with 39 rows and 369 columns
    ##  [6] epic: matrix with 8 rows and 369 columns
    ##  [7] quantiseq: matrix with 11 rows and 369 columns
    ##  [8] mcp: matrix with 11 rows and 369 columns
    ##  [9] scores: matrix with 4 rows and 369 columns
    ##  [10] purity: matrix with 1 rows and 320 columns
    ## Functionality:
    ##  experiments() - obtain the ExperimentList instance
    ##  colData() - the primary/phenotype DataFrame
    ##  sampleMap() - the sample coordination DataFrame
    ##  `$`, `[`, `[[` - extract colData columns, subset, or experiment
    ##  *Format() - convert into a long or wide DataFrame
    ##  assays() - convert ExperimentList to a SimpleList of matrices
    ##  exportClass() - save data to flat files

    curatedPCaData::mae_tcga[["gex.fpkm"]][1:4,1:4]
    ##           TCGA.EJ.5505.01 TCGA.HC.A8D0.01 TCGA.VN.A88I.01 TCGA.KK.A8IA.01
    ## 5_8S_rRNA          0.0000          0.0000          0.0000          0.0000
    ## 5S_rRNA            0.0000          0.0000          0.0000          0.0000
    ## 7SK                0.0000          0.0000          0.0000          0.0000
    ## A1BG               9.5252         10.0481         10.9035         10.7319
    curatedPCaData::mae_tcga[["cna.gistic"]][1:4,1:4]
    ##       TCGA.HC.A632.01 TCGA.HC.8213.01 TCGA.HC.8216.01 TCGA.EJ.A65G.01
    ## A1BG                0               0               0               0
    ## A1CF                0               0               0               0
    ## A2M                -1               0               0               0
    ## A2ML1              -1               0               0               0
    MultiAssayExperiment::colData(curatedPCaData::mae_tcga)[1:3,1:5]
    ## DataFrame with 3 rows and 5 columns
    ##               study_name   patient_id     sample_name        alt_sample_name overall_survival_status
    ##              <character>  <character>     <character>            <character>               <numeric>
    ## TCGA.EJ.5502        TCGA TCGA.EJ.5502 TCGA.EJ.5502.01 d36a0b54-6e09-4ee2-9..                       0
    ## TCGA.YJ.A8SW        TCGA TCGA.YJ.A8SW TCGA.YJ.A8SW.01 8A96A7A8-0413-42B9-9..                       0
    ## TCGA.EJ.5525        TCGA TCGA.EJ.5525 TCGA.EJ.5525.01 351af15c-b213-4621-8..                       0


    curatedPCaData::mae_taylor
    ## A MultiAssayExperiment object of 11 listed
    ##  experiments with user-defined names and respective classes.
    ##  Containing an ExperimentList class object of length 11:
    ##  [1] cna.gistic: data.frame with 16715 rows and 194 columns
    ##  [2] cna.logr: matrix with 22419 rows and 218 columns
    ##  [3] gex.rma: matrix with 17410 rows and 179 columns
    ##  [4] mut: RaggedExperiment with 319 rows and 43 columns
    ##  [5] xcell: matrix with 39 rows and 179 columns
    ##  [6] epic: matrix with 8 rows and 179 columns
    ##  [7] quantiseq: matrix with 11 rows and 179 columns
    ##  [8] mcp: matrix with 11 rows and 179 columns
    ##  [9] scores: matrix with 4 rows and 179 columns
    ##  [10] purity: matrix with 1 rows and 150 columns
    ##  [11] cibersort: matrix with 22 rows and 179 columns
    ## Functionality:
    ##  experiments() - obtain the ExperimentList instance
    ##  colData() - the primary/phenotype DataFrame
    ##  sampleMap() - the sample coordination DataFrame
    ##  `$`, `[`, `[[` - extract colData columns, subset, or experiment
    ##  *Format() - convert into a long or wide DataFrame
    ##  assays() - convert ExperimentList to a SimpleList of matrices
    ##  exportClass() - save data to flat files
    curatedPCaData::mae_sun
    ## A MultiAssayExperiment object of 7 listed
    ##  experiments with user-defined names and respective classes.
    ##  Containing an ExperimentList class object of length 7:
    ##  [1] gex.rma: matrix with 12798 rows and 79 columns
    ##  [2] cibersort: matrix with 22 rows and 79 columns
    ##  [3] xcell: matrix with 39 rows and 79 columns
    ##  [4] epic: matrix with 8 rows and 79 columns
    ##  [5] quantiseq: matrix with 11 rows and 79 columns
    ##  [6] mcp: matrix with 11 rows and 79 columns
    ##  [7] scores: matrix with 4 rows and 79 columns
    ## Functionality:
    ##  experiments() - obtain the ExperimentList instance
    ##  colData() - the primary/phenotype DataFrame
    ##  sampleMap() - the sample coordination DataFrame
    ##  `$`, `[`, `[[` - extract colData columns, subset, or experiment
    ##  *Format() - convert into a long or wide DataFrame
    ##  assays() - convert ExperimentList to a SimpleList of matrices
    ##  exportClass() - save data to flat files

## R Shiny

A web interface built with R Shiny for `curatedPCaData` can be launched
via:

    curatedPCaData::shiny()

## Known issues

Large R-packages directly installed using `devtools::install_github`
sometimes result in error:

    Error in utils::download.file( ...
        download from ' ... ' failed.

There are few options to fix this:

### Adjusting download options

For Windows users, it’s known that setting

    options(download.file.method = "wininet")

may fix this GitHub direct download/install issue. Alternatively,
especially in other OSes, the following download method may fix this
issue:

    options(download.file.method = "libcurl")

### Direct cloning, building, or downloading and installing the package tarball

A more hands-on approach is to use Git to clone the package and then
build and install it using R tools. In your terminal or command prompt,
go to a suitable root directory for git repositories:

    git clone https://github.com/Syksy/curatedPCaData.git
    R CMD build curatedPCaData

alternatively, to speed up package tarball building, add parameter
`--no-build-vignettes`. This will produce a file
`curatedPCaData_x.y.z.tar.gz` in your active directory, where `x.y.z`
correspond to the current package version.

Alternatively, if one does not wish to clone the git and build the
tarball, it’s also possible to directly download latest release tarball
from
<a href="https://github.com/Syksy/curatedPCaData/releases" class="uri">https://github.com/Syksy/curatedPCaData/releases</a>

After this, the `curatedPCaData`R-package tarball can be installed
using:

    R CMD INSTALL curatedPCaData_x.y.z.tar.gz

Please note that some dependencies (such as the packages
`MultiAssayExperiment` and `S4Vectors`) may produce an error during
installation if they are not found for R. In this case these
dependencies need to be installed from their respective repositories
such as CRAN or BioConductor.
