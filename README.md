
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RiskMap

<!-- badges: start -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/RiskMap)](https://cran.r-project.org/package=RiskMap)
[![rstudio mirror
downloads](https://cranlogs.r-pkg.org/badges/RiksMap)](https://github.com/r-hub/cranlogs.app)
[![](https://cranlogs.r-pkg.org/badges/grand-total/RiskMap)](https://cran.r-project.org/package=stplanr)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
<!-- badges: end -->

## Overview

`RiskMap` provides functions for geo-statistical analysis of both
continuous and count data using maximum likelihood methods. The models
implemented in the package use stationary Gaussian processes with Matern
correlation function to carry out spatial prediction in a geographical
area of interest. The underpinning theory of the methods implemented in
the package are found in [Diggle and Giorgi
(2019)](https://www.routledge.com/Model-based-Geostatistics-for-Global-Public-Health-Methods-and-Applications/Diggle-Giorgi/p/book/9781032093642).

## Installation

To install the stable version, use:

``` r
install.packages("RiskMap")
```

The development version can be installed using **devtools**:

``` r
# install.packages("devtools") # if not already installed
devtools::install_github("claudiofronterre/RiskMap")
```
