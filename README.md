---
output: github_document
---



# curatedPCaData <img src="man/figures/hex.png" align="right" height="139" />

<!-- badges: start -->
<!-- badges: end -->

## Overview 

`curatedPCaData` is a collection of publically available and annotated data resources
concerning prostate cancer. 

## Installation

You can install `curatedPCaData` from github with: 

```r
# install.packages("devtools")
devtools::install_github("Syksy/curatedPCaData")
```

<!--- add BioConductor once up --->

## Usage


```r
library(curatedPCaData)
## Loading required package: S4Vectors
## Loading required package: stats4
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ, clusterExport, clusterMap, parApply,
##     parCapply, parLapply, parLapplyLB, parRapply, parSapply, parSapplyLB
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind, colnames, dirname, do.call, duplicated, eval,
##     evalq, Filter, Find, get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank, rbind, Reduce, rownames, sapply,
##     setdiff, sort, table, tapply, union, unique, unsplit, which, which.max, which.min
## 
## Attaching package: 'S4Vectors'
## The following object is masked from 'package:base':
## 
##     expand.grid
## Loading required package: MultiAssayExperiment
## Loading required package: SummarizedExperiment
## Loading required package: GenomicRanges
## Loading required package: IRanges
## 
## Attaching package: 'IRanges'
## The following object is masked from 'package:grDevices':
## 
##     windows
## Loading required package: GenomeInfoDb
## Loading required package: Biobase
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with 'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## Loading required package: DelayedArray
## Loading required package: matrixStats
## 
## Attaching package: 'matrixStats'
## The following objects are masked from 'package:Biobase':
## 
##     anyMissing, rowMedians
## 
## Attaching package: 'DelayedArray'
## The following objects are masked from 'package:matrixStats':
## 
##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges
## The following objects are masked from 'package:base':
## 
##     aperm, apply, rowsum

curatedPCaData::mae_tcga
## A MultiAssayExperiment object of 2 listed
##  experiments with user-defined names and respective classes.
##  Containing an ExperimentList class object of length 2:
##  [1] cna: matrix with 21735 rows and 333 columns
##  [2] gex: matrix with 19958 rows and 333 columns
## Features:
##  experiments() - obtain the ExperimentList instance
##  colData() - the primary/phenotype DFrame
##  sampleMap() - the sample availability DFrame
##  `$`, `[`, `[[` - extract colData columns, subset, or experiment
##  *Format() - convert into a long or wide DFrame
##  assays() - convert ExperimentList to a SimpleList of matrices

curatedPCaData::mae_tcga[["gex"]][1:4,1:4]
##          TCGA.EJ.5502.01 TCGA.HC.7209.01 TCGA.HC.7748.01 TCGA.J4.A83N.01
## A1BG             34.3994         15.4859         25.0083         29.4705
## A1BG.AS1         27.8727         14.6187         27.6506         17.3826
## A1CF              0.0000          0.0000          0.2645          0.0000
## A2M           20704.9666       8330.9059      14777.1633       4068.3117
curatedPCaData::mae_tcga[["cna"]][1:4,1:4]
##          TCGA.EJ.5502.01 TCGA.HC.7209.01 TCGA.HC.7748.01 TCGA.J4.A83N.01
## A1BG               0.003          -0.013           0.007          -0.006
## A1BG.AS1           0.003          -0.013           0.007          -0.006
## A1CF               0.000          -0.040           0.012           0.010
## A2M                0.002           0.009          -0.021          -0.001
MultiAssayExperiment::colData(curatedPCaData::mae_tcga)[1:3,1:5]
## DataFrame with 3 rows and 5 columns
##                     study_name   patient_id     sample_name                      alt_sample_name
##                    <character>  <character>     <character>                          <character>
## TCGA.EJ.5502 TCGA, provisional TCGA.EJ.5502 TCGA.EJ.5502.01 d36a0b54-6e09-4ee2-923b-c1516d78bb03
## TCGA.YJ.A8SW TCGA, provisional TCGA.YJ.A8SW TCGA.YJ.A8SW.01 8A96A7A8-0413-42B9-9173-FD63761DD83A
## TCGA.EJ.5525 TCGA, provisional TCGA.EJ.5525 TCGA.EJ.5525.01 351af15c-b213-4621-8bcb-3f4ddcf72553
##              overall_survival_status
##                            <numeric>
## TCGA.EJ.5502                       0
## TCGA.YJ.A8SW                       0
## TCGA.EJ.5525                       0

curatedPCaData::mae_taylor
## A MultiAssayExperiment object of 4 listed
##  experiments with user-defined names and respective classes.
##  Containing an ExperimentList class object of length 4:
##  [1] cna_cbio: matrix with 17536 rows and 166 columns
##  [2] cna_geo: matrix with 22431 rows and 0 columns
##  [3] gex_cbio: matrix with 21312 rows and 133 columns
##  [4] gex_geo: matrix with 17638 rows and 0 columns
## Features:
##  experiments() - obtain the ExperimentList instance
##  colData() - the primary/phenotype DFrame
##  sampleMap() - the sample availability DFrame
##  `$`, `[`, `[[` - extract colData columns, subset, or experiment
##  *Format() - convert into a long or wide DFrame
##  assays() - convert ExperimentList to a SimpleList of matrices

curatedPCaData::mae_sun
## A MultiAssayExperiment object of 1 listed
##  experiment with a user-defined name and respective class.
##  Containing an ExperimentList class object of length 1:
##  [1] gex: matrix with 12057 rows and 79 columns
## Features:
##  experiments() - obtain the ExperimentList instance
##  colData() - the primary/phenotype DFrame
##  sampleMap() - the sample availability DFrame
##  `$`, `[`, `[[` - extract colData columns, subset, or experiment
##  *Format() - convert into a long or wide DFrame
##  assays() - convert ExperimentList to a SimpleList of matrices
```
