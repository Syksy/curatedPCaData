
# curatedPCaData <img src="man/figures/hex.png" align="right" height="139" />

<!-- badges: start -->
<!-- badges: end -->

## Overview

`curatedPCaData` is a collection of publically available and annotated
data resources concerning prostate cancer.

## Citation

If you use `curatedPCaData`, please consider adding the following
citation (bioRxiv preprint, [direct link
here](https://www.biorxiv.org/content/10.1101/2023.01.17.524403v1)):

    @article {Laajala2023.01.17.524403,
        author = {Laajala, Teemu D and Sreekanth, Varsha and Soupir, Alex and Creed, Jordan and Halkola, Anni S and Calboli, Federico CF and Singaravelu, Kalaimathy and Orman, Michael and Colin-Leitzinger, Christelle and Gerke, Travis and Fidley, Brooke L. and Tyekucheva, Svitlana and Costello, James C},
        title = {curatedPCaData: Integration of clinical, genomic, and signature features in a curated and harmonized prostate cancer data resource},
        year = {2023},
        doi = {10.1101/2023.01.17.524403},
        URL = {https://www.biorxiv.org/content/early/2023/01/19/2023.01.17.524403},
        eprint = {https://www.biorxiv.org/content/early/2023/01/19/2023.01.17.524403.full.pdf},
        journal = {bioRxiv}
    }

## Installation

### Bioconductor installation

In order to install the package from Bioconductor, make sure
`BiocManager` is installed and then call the function to install
`curatedPCaData`:

    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    BiocManager::install("curatedPCaData")

### GitHub installation

A download link to the latest pre-built `curatedPCaData` tarball is
available on the right-side in GitHub under
[Releases](https://github.com/Syksy/curatedPCaData/releases).

You can also install `curatedPCaData` from GitHub inside R with:

    # install.packages("devtools")
    devtools::install_github("Syksy/curatedPCaData")

To build the package tarball from a cloned git repo, run the following
in terminal / command prompt while in the root of the project:

    R CMD build curatedPCaData

It is then possible to install the self-built tarball:

    R CMD INSTALL curatedPCaData_x.y.z.tar.gz

Note that building the package locally will require dependencies to be
present for the R installation.

## Usage

### Vignettes

`curatedPCaData` delivers with a number of vignettes displaying the
package’s generic use as well as summaries and application examples
across the prostate cancer datasets provided there-in. The vignette
`overview` is intended for gaining a first-line comprehensive view into
the package’s contents.

The vignettes can be accessed via `vignette(package = "curatedPCaData")`
or via `?curatedPCaData` section ‘*User guides, package vignettes and
other documentation*’.

A list of available vignettes, subset to suitable topics:

    tools::getVignetteInfo("curatedPCaData")[, c("Topic", "Title")]
    ##      Topic      Title                                
    ## [1,] "analyses" "Analysis examples in curatedPCaData"
    ## [2,] "overview" "Overview to curatedPCaData"

### Downloading data

The function `getPCa` with its main parameter with study shortname is
the primary means of extracting data from a cohort. It will
automatically create a `MultiAssayExperiment`-object of the study:

    library(curatedPCaData)

    mae_tcga <- getPCa("tcga")

    class(mae_tcga)
    ## [1] "MultiAssayExperiment"
    ## attr(,"package")
    ## [1] "MultiAssayExperiment"
    names(mae_tcga)
    ##  [1] "cna.gistic"   "gex.rsem.log" "mut"          "cibersort"    "xcell"        "epic"         "quantiseq"    "mcp"          "estimate"     "scores"

### Brief example

Simple example use of curated datasets and ’omics there-in:

    mae_taylor <- getPCa("taylor")
    mae_sun <- getPCa("sun")

    mae_tcga
    ## A MultiAssayExperiment object of 10 listed
    ##  experiments with user-defined names and respective classes.
    ##  Containing an ExperimentList class object of length 10:
    ##  [1] cna.gistic: matrix with 23151 rows and 492 columns
    ##  [2] gex.rsem.log: matrix with 19658 rows and 461 columns
    ##  [3] mut: RaggedExperiment with 30897 rows and 495 columns
    ##  [4] cibersort: matrix with 22 rows and 461 columns
    ##  [5] xcell: matrix with 39 rows and 461 columns
    ##  [6] epic: matrix with 8 rows and 461 columns
    ##  [7] quantiseq: matrix with 11 rows and 461 columns
    ##  [8] mcp: matrix with 11 rows and 461 columns
    ##  [9] estimate: matrix with 4 rows and 461 columns
    ##  [10] scores: matrix with 4 rows and 461 columns
    ## Functionality:
    ##  experiments() - obtain the ExperimentList instance
    ##  colData() - the primary/phenotype DataFrame
    ##  sampleMap() - the sample coordination DataFrame
    ##  `$`, `[`, `[[` - extract colData columns, subset, or experiment
    ##  *Format() - convert into a long or wide DataFrame
    ##  assays() - convert ExperimentList to a SimpleList of matrices
    ##  exportClass() - save data to flat files

    mae_tcga[["gex.rsem.log"]][1:4, 1:4]
    ##          TCGA.G9.6348.01 TCGA.CH.5766.01 TCGA.EJ.A65G.01 TCGA.EJ.5527.01
    ## A1BG              4.3733          6.0244          7.4927          3.7801
    ## A1BG-AS1          4.5576          6.3326          6.7861          4.5912
    ## A1CF              0.4008          0.7574          0.0000          0.0000
    ## A2M              14.3952         12.8331         12.5017         14.2289

    mae_tcga[["cna.gistic"]][1:4, 1:4]
    ##       TCGA.2A.A8VL.01 TCGA.2A.A8VO.01 TCGA.2A.A8VT.01 TCGA.2A.A8VV.01
    ## A1BG                0               0               0               0
    ## A1CF                0               0              -1               0
    ## A2M                 0               0              -1               0
    ## A2ML1               0               0              -1               0

    colData(mae_tcga)[1:3, 1:5]
    ## DataFrame with 3 rows and 5 columns
    ##                  study_name   patient_id     sample_name        alt_sample_name overall_survival_status
    ##                 <character>  <character>     <character>            <character>               <integer>
    ## TCGA.2A.A8VL.01        TCGA TCGA.2A.A8VL TCGA.2A.A8VL.01 F9F392D3-E3C0-4CF2-A..                       0
    ## TCGA.2A.A8VO.01        TCGA TCGA.2A.A8VO TCGA.2A.A8VO.01 0BD35529-3416-42DD-A..                       0
    ## TCGA.2A.A8VT.01        TCGA TCGA.2A.A8VT TCGA.2A.A8VT.01 BFECF807-0658-417B-9..                       0

    mae_taylor
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
    ##  [10] estimate: matrix with 4 rows and 179 columns
    ##  [11] scores: matrix with 4 rows and 179 columns
    ## Functionality:
    ##  experiments() - obtain the ExperimentList instance
    ##  colData() - the primary/phenotype DataFrame
    ##  sampleMap() - the sample coordination DataFrame
    ##  `$`, `[`, `[[` - extract colData columns, subset, or experiment
    ##  *Format() - convert into a long or wide DataFrame
    ##  assays() - convert ExperimentList to a SimpleList of matrices
    ##  exportClass() - save data to flat files

    mae_sun
    ## A MultiAssayExperiment object of 8 listed
    ##  experiments with user-defined names and respective classes.
    ##  Containing an ExperimentList class object of length 8:
    ##  [1] gex.rma: matrix with 12784 rows and 79 columns
    ##  [2] cibersort: matrix with 22 rows and 79 columns
    ##  [3] xcell: matrix with 39 rows and 79 columns
    ##  [4] epic: matrix with 8 rows and 79 columns
    ##  [5] quantiseq: matrix with 11 rows and 79 columns
    ##  [6] estimate: matrix with 4 rows and 79 columns
    ##  [7] scores: matrix with 4 rows and 79 columns
    ##  [8] mcp: matrix with 11 rows and 79 columns
    ## Functionality:
    ##  experiments() - obtain the ExperimentList instance
    ##  colData() - the primary/phenotype DataFrame
    ##  sampleMap() - the sample coordination DataFrame
    ##  `$`, `[`, `[[` - extract colData columns, subset, or experiment
    ##  *Format() - convert into a long or wide DataFrame
    ##  assays() - convert ExperimentList to a SimpleList of matrices
    ##  exportClass() - save data to flat files

For further hands-on examples, please see for example the
`analyses`-vignette.
