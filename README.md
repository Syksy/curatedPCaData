---
output: github_document
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

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

Above may fail depending on the connection stability, as the package is relatively large. In this case it's best to install via the tarball.
To build the package tarball from a cloned git repo, run the following in terminal / command prompt while in the root of the project:

```
R CMD build curatedPCaData
```

One install the self-built tarball or download the a premade latest package tarball from the releases page and install it with:

```
R CMD INSTALL curatedPCaData_x.y.z.tar.gz
```


<!--- add BioConductor once up --->

## Usage

Simple example use of curated datasets and 'omics there-in:

```{r example_one, warning = FALSE, message = FALSE}

library(curatedPCaData)


curatedPCaData::mae_tcga

curatedPCaData::mae_tcga[["gex"]][1:4,1:4]
curatedPCaData::mae_tcga[["cna"]][1:4,1:4]
MultiAssayExperiment::colData(curatedPCaData::mae_tcga)[1:3,1:5]


curatedPCaData::mae_taylor
curatedPCaData::mae_sun

```

## R Shiny

A web interface built with R Shiny for `curatedPCaData` can be launched via:

```{r example_one, warning = FALSE, message = FALSE}

curatedPCaData::shiny()

```

## Known issues

### 

Large R-packages directly installed using ```devtools::install_github``` sometimes result in error:

```
Error in utils::download.file( ...
    download from ' ... ' failed.
``` 

There are few options to fix this:

#### Adjusting download options

For Windows users, it's known that setting

```
options(download.file.method = "wininet")
```

may fix this GitHub direct download/install issue. Alternatively, especially in other OSes, the following download method may fix this issue:

```
options(download.file.method = "libcurl")
```

#### Direct cloning, building, or downloading and installing the package tarball

A more hands-on approach is to use Git to clone the package and then build and install it using R tools. In your terminal or command prompt, go to a suitable root directory for git repositories:

```
git clone https://github.com/Syksy/curatedPCaData.git
R CMD build curatedPCaData
``` 

alternatively, to speed up package tarball building, add parameter ```--no-build-vignettes```. This will produce a file ```curatedPCaData_x.y.z.tar.gz``` in your active directory, where ```x.y.z``` correspond to the current package version. 

Alternatively, if one does not wish to clone the git and build the tarball, it's also possible to directly download latest release tarball from https://github.com/Syksy/curatedPCaData/releases

After this, the ```curatedPCaData```R-package tarball can be installed using:

```
R CMD INSTALL curatedPCaData_x.y.z.tar.gz
```

Please note that some dependencies (such as the packages ```MultiAssayExperiment``` and ```S4Vectors```) may produce an error during installation if they are not found for R. In this case these dependencies need to be installed from their respective repositories such as CRAN or BioConductor.
