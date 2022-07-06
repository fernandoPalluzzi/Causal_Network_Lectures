# Causal Network Lectures
Collection of Causal Network Analysis lectures.
&nbsp;

## Main libraries installation

This lecture is for the R environment version >= 4.0. Reference repo for SEMgraph [HERE](https://github.com/fernandoPalluzzi/SEMgraph)

```r
# Library devtools is required to download the development versions
# from GitHub.

install.packages("devtools")

library(devtools)

# SEMgraph and SEMdata installation

devtools::install_github("fernandoPalluzzi/SEMgraph")

devtools::install_github("fernandoPalluzzi/SEMdata")

# By loading SEMgraph, libraries igraph (network analysis) and 
# lavaan (structural equation modeling) should be automatically loaded.

# Load required libraries
library(SEMgraph)

# Load required datasets
library(SEMdata)
```
&nbsp;

## Additional R packages

```r
# Homo sapiens annotation data (Bioconductor package)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

library(org.Hs.eg.db)


# Applied regression functions: q-q plot (CRAN package)

install.packages("car")

library(car)


# Nonparanormal transform for non-Gaussian data (CRAN package)

install.packages("huge")

library(huge)
```
&nbsp;

## Troubleshooting

To install from CRAN, you can use:

```r
install.packages("PACKAGE")
```

Just substitute "PACKAGE" with the actual CRAN package name.
Here is the list of CRAN dependencies:

- igraph (>= 1.2.1)
- lavaan (>= 0.5-23)
- boot (≥ 1.3-25)
- cate (≥ 1.0.4)
- corpcor (≥ 1.6.9)
- dagitty (≥ 0.3-0)
- diffusr (≥ 0.1.4)
- flip (≥ 2.5.0)
- gdata (≥ 2.18.0)
- ggm (≥ 2.3)
- GGMncv (≥ 2.1.0)
- glmnet (≥ 2.0-18)
- graph (≥ 1.56.0)
- leaps (≥ 3.1)
- Matrix (≥ 1.3-0)
- pbapply (≥ 1.4-3)
- protoclust (≥ 1.6.3)
- RcppEigen (≥ 0.3.3.4.0)

To install from Bioconductor, you can use:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("PACKAGE")
```

Just substitute "PACKAGE" with the actual Bioconductor package name.

Here is the list of Bioconductor dependencies:
- Rgraphviz (≥ 2.22.0)
