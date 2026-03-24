# Dynamic Enrichment for Evaluation of Protein Networks

This is the accompanying `R` package for the [DEEPN](https://github.com/emptyewer/DEEPN) project.

## Installation

`deepn` has several package dependencies, most notably [JAGS](http://mcmc-jags.sourceforge.net) and [Bioconductor](https://www.bioconductor.org).  JAGS must be setup outside `R` (see link for instructions), while Bioconductor can be installed within `R` via:

```r
source("https://bioconductor.org/biocLite.R")
biocLite()
```

Once that is done, you can install `deepn` via `remotes` (which itself can be installed via `install.packages("remotes")`):

```r
remotes::install_github("pbreheny/deepn")
```

## Basic usage

Once installed, `deepn` can be used to analyze the output from `deepn` with

```r
require(deepn)
analyzeDeepn("deepn.config")
```

where `"deepn.config"` is the name of the configuration file that DEEPN outputs containing path names to all the `.csv` files needed in the analysis.  This will create a file, `stat.csv`, containing a table of summary statistics for each gene.

The above analysis may take a while to run, depending on the threshold.  To run a "quick" version of the analysis, you can use

```r
analyzeDeepn("deepn.config", debug=TRUE)
```

which applies a very high PPM threshold of 200 prior to analysis, so that it should only take a minute or two to complete.
