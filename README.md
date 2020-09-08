
# curatedPCaData <img src="man/figures/hex.png" align="right" height="139" />

<!-- badges: start -->

<!-- badges: end -->

## Overview

`curatedPCaData` is a collection of publically available and annotated
data resources concerning prostate cancer.

## Installation

You can install `curatedPCaData` from github with:

``` r
# install.packages("devtools")
devtools::install_github("Syksy/curatedPCaData")
```

<!--- add BioConductor once up --->

## Usage

``` r
library(curatedPCaData)
##Loading required package: S4Vectors
##Loading required package: stats4
##Loading required package: BiocGenerics
##Loading required package: parallel
##
##Attaching package: ‘BiocGenerics’
##
##The following objects are masked from ‘package:parallel’:
##
##    clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##    clusterExport, clusterMap, parApply, parCapply, parLapply,
##    parLapplyLB, parRapply, parSapply, parSapplyLB
##
##The following objects are masked from ‘package:stats’:
##
##    IQR, mad, sd, var, xtabs
##
##The following objects are masked from ‘package:base’:
##
##    anyDuplicated, append, as.data.frame, basename, cbind, colnames,
##    dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
##    grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
##    order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##    rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
##    union, unique, unsplit, which, which.max, which.min
##
##
##Attaching package: ‘S4Vectors’
##
##The following object is masked from ‘package:base’:
##
##    expand.grid    
    
curatedPCaData::mae_tcga
##A MultiAssayExperiment object of 2 listed
## experiments with user-defined names and respective classes.
## Containing an ExperimentList class object of length 2:
## [1] GEX: matrix with 19958 rows and 333 columns
## [2] CNA: matrix with 21735 rows and 333 columns
##Features:
## experiments() - obtain the ExperimentList instance
## colData() - the primary/phenotype DFrame
## sampleMap() - the sample availability DFrame
## `$`, `[`, `[[` - extract colData columns, subset, or experiment
## *Format() - convert into a long or wide DFrame
## assays() - convert ExperimentList to a SimpleList of matrices
 