
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
## Loading required package: S4Vectors
## Loading required package: stats4
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
##     union, unique, unsplit, which, which.max, which.min
## 
## Attaching package: 'S4Vectors'
## The following object is masked from 'package:base':
## 
##     expand.grid

curatedPCaData::tcga_gex[1:5,1:5]
##                      AATF ABHD17AP4 ACTR3BP5 BCRP5 C20ORF27
## TCGA.EJ.5502.01  926.7640       NaN      NaN   NaN 302.4839
## TCGA.HC.7209.01 1063.1049       NaN      NaN   NaN 304.6845
## TCGA.HC.7748.01  848.2444       NaN      NaN   NaN 251.8019
## TCGA.J4.A83N.01 1092.4076       NaN      NaN   NaN 431.5684
## TCGA.2A.A8VV.01 1223.2027       NaN      NaN   NaN 433.4628
```

Project TODO check-list
* Functional MAE-objects
* Double-check N-counts and sample lists (e.g. TCGA N=333)
* 3 example datasets
* Instead of pregenerated *_pdata-txt files, have original creation scripts embedded
* Double-check relative paths inside package
* Explore the additional curated fields / software / packages:
* > Neoantigen load (would require BAM/VCF in most software?)
* > Cell composition (xCell)
* > Zone of Origin
* > ...