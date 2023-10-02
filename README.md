
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RiskMap

<!-- badges: start -->
<!-- badges: end -->

The goal of `RiskMap` is to provide a set of functions for
visualisation, processing and likelihood-based analysis of
geostatistical data.

## Installation

You can install the development version of RiskMap from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("claudiofronterre/RiskMap")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(RiskMap)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.

## Estimation function design

How should the user specify the model?

``` r
# OPTION 1
y ~ rainfall + gp(x, y) + re(id_school) + re(id_region)

# Arguments for gp function
gp(long = NULL, lat = NULL, kappa = (numeric_value, default = 0.5), nugget = c(T = default, F, fixed_numeric_value), ...)

# Arguments for re function
re(numeric or categorical variable, ...) only needs an index in the dataset

glgm <- function(formula,
                 distr_offset = NULL,
                 cov_offset = NULL,
                 data,
                 family,
                 convert_to_crs = NULL,
                 scale_to_km = TRUE,
                 control_MCMC = NULL,
                 S_samples = NULL,
                 save_samples = F,
                 messages = TRUE) 
```
My solution to incorporate "gp" into the formula

``` r
(tf <- terms(y ~ x + x:z + gp(kappa = 0.5, nugget = TRUE)+w, specials = "gp"))
attr(tf, "specials")    # index 's' variable(s)
gp <- rownames(attr(tf, "factors"))[[attr(tf, "specials")$gp]]
gp_list <- eval(parse(text = gsub("gp","list",gp)))
gp_list
```
