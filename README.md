
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

# insert basic example here 
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