
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

It is possible to install the self-built tarball or download a premade
latest package with (see Releases-page on GitHub):

    R CMD INSTALL curatedPCaData_x.y.z.tar.gz

<!--- add BioConductor once up --->

## Usage

Simple example use of curated datasets and ’omics there-in:


    library(curatedPCaData)

    curatedPCaData::mae_tcga
    ## A MultiAssayExperiment object of 10 listed
    ##  experiments with user-defined names and respective classes.
    ##  Containing an ExperimentList class object of length 10:
    ##  [1] cna.gistic: matrix with 23151 rows and 404 columns
    ##  [2] gex.rsem.log: matrix with 19658 rows and 461 columns
    ##  [3] mut: RaggedExperiment with 15389 rows and 371 columns
    ##  [4] cibersort: matrix with 22 rows and 417 columns
    ##  [5] xcell: matrix with 39 rows and 417 columns
    ##  [6] epic: matrix with 8 rows and 417 columns
    ##  [7] quantiseq: matrix with 11 rows and 417 columns
    ##  [8] mcp: matrix with 11 rows and 417 columns
    ##  [9] scores: matrix with 4 rows and 417 columns
    ##  [10] estimate: data.frame with 4 rows and 417 columns
    ## Functionality:
    ##  experiments() - obtain the ExperimentList instance
    ##  colData() - the primary/phenotype DataFrame
    ##  sampleMap() - the sample coordination DataFrame
    ##  `$`, `[`, `[[` - extract colData columns, subset, or experiment
    ##  *Format() - convert into a long or wide DataFrame
    ##  assays() - convert ExperimentList to a SimpleList of matrices
    ##  exportClass() - save data to flat files

    curatedPCaData::mae_tcga[["gex.rsem.log"]][1:4,1:4]
    ##           TCGA.G9.6348.01 TCGA.CH.5766.01 TCGA.EJ.A65G.01 TCGA.EJ.5527.01
    ## ARHGEF10L          8.8729          8.5581          9.2085          8.7699
    ## HIF3A              5.9049          4.9716          6.7795          5.5978
    ## RNF17              0.4008          0.7574          0.0000          2.5554
    ## RNF10             12.3538         12.2950         11.9701         11.7983

    curatedPCaData::mae_tcga[["cna.gistic"]][1:4,1:4]
    ##         TCGA.2A.A8VL.01 TCGA.2A.A8VO.01 TCGA.2A.A8VT.01 TCGA.2A.A8VV.01
    ## ACAP3                 0               0               0               0
    ## ACTRT2                0               0               0               0
    ## AGRN                  0               0               0               0
    ## ANKRD65               0               0               0               0

    MultiAssayExperiment::colData(curatedPCaData::mae_tcga)[1:3,1:5]
    ## DataFrame with 3 rows and 5 columns
    ##               study_name   patient_id     sample_name        alt_sample_name overall_survival_status
    ##              <character>  <character>     <character>            <character>               <integer>
    ## TCGA.2A.A8VL        TCGA TCGA.2A.A8VL TCGA.2A.A8VL.01 F9F392D3-E3C0-4CF2-A..                       0
    ## TCGA.2A.A8VO        TCGA TCGA.2A.A8VO TCGA.2A.A8VO.01 0BD35529-3416-42DD-A..                       0
    ## TCGA.2A.A8VT        TCGA TCGA.2A.A8VT TCGA.2A.A8VT.01 BFECF807-0658-417B-9..                       0

    curatedPCaData::mae_taylor
    ## A MultiAssayExperiment object of 11 listed
    ##  experiments with user-defined names and respective classes.
    ##  Containing an ExperimentList class object of length 11:
    ##  [1] cna.gistic: matrix with 17832 rows and 194 columns
    ##  [2] cna.logr: matrix with 18062 rows and 218 columns
    ##  [3] gex.rma: matrix with 17410 rows and 179 columns
    ##  [4] mut: RaggedExperiment with 90 rows and 43 columns
    ##  [5] cibersort: matrix with 22 rows and 179 columns
    ##  [6] xcell: matrix with 39 rows and 179 columns
    ##  [7] epic: matrix with 8 rows and 179 columns
    ##  [8] quantiseq: matrix with 11 rows and 179 columns
    ##  [9] mcp: matrix with 11 rows and 179 columns
    ##  [10] scores: matrix with 4 rows and 179 columns
    ##  [11] estimate: data.frame with 4 rows and 179 columns
    ## Functionality:
    ##  experiments() - obtain the ExperimentList instance
    ##  colData() - the primary/phenotype DataFrame
    ##  sampleMap() - the sample coordination DataFrame
    ##  `$`, `[`, `[[` - extract colData columns, subset, or experiment
    ##  *Format() - convert into a long or wide DataFrame
    ##  assays() - convert ExperimentList to a SimpleList of matrices
    ##  exportClass() - save data to flat files

    curatedPCaData::mae_sun
    ## A MultiAssayExperiment object of 8 listed
    ##  experiments with user-defined names and respective classes.
    ##  Containing an ExperimentList class object of length 8:
    ##  [1] gex.rma: matrix with 12784 rows and 79 columns
    ##  [2] cibersort: matrix with 22 rows and 79 columns
    ##  [3] xcell: matrix with 39 rows and 79 columns
    ##  [4] epic: matrix with 8 rows and 79 columns
    ##  [5] quantiseq: matrix with 11 rows and 79 columns
    ##  [6] mcp: matrix with 11 rows and 79 columns
    ##  [7] scores: matrix with 4 rows and 79 columns
    ##  [8] estimate: data.frame with 4 rows and 79 columns
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
