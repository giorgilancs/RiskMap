
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
gp(x, y, kappa = (numeric_value, default = 0.5), nugget = c(T = default, F, fixed_numeric_value), ...)

# Arguments for re function
# re(numeric or categorical variable, ...) only needs an index in the dataset

glgm <- function(formula,
                 m_offset = NULL,
                 data,
                 family,
                 convert_to_crs = NULL,
                 scale_to_km = TRUE,
                 control_MCMC = NULL,
                 kappa = 0.5,
                 start_beta = NULL,
                 start_cov_pars = NULL,
                 start_vars.re = NULL,
                 messages = TRUE,
                 plot_correlogram = TRUE,
                 fix_tau2 = NULL) {

  if(family=="binomial" | family=="poisson") {
    if(is.null(control_MCMC)) stop("if family='binomial' or family='poisson'
                                   'control_MCMC' must be provided.")
    if(class(control_MCMC)!="mcmc.RiskMap") stop("'control_MCMC' must be of class 'mcmc.PrevMap'")
  }

  if(class(formula)!="formula") stop("'formula' must be a 'formula'
                                     object indicating the variables of the
                                     model to be fitted")

  if(class(data)[1]!="sf") stop("'data' must be an object of class 'sf'")
  if(is.na(st_crs(data))) stop("the CRS of the data is missing
                               and must be specified; see ?st_crs")

  if(!is.null(hr_re)) {
    if(class(hr_re)!="formula") stop("the object 'he_re' must be a 'formula' object")
  }

  if(family=="Binomial") {
    if(is.null(m_offset)) stop("if family='binomial', the argument 'm_offset'
                               must be provided")
    m_offset_ch <- as.charaster(as.name(subsitute(m_offset)))
  }

  if(kappa < 0) stop("kappa must be positive.")


  if(geo_re != "GP" & geo_re != "GP+NUG") stop("'geo_re' must be 'GP' or
                                               'GP+NUG'")

  if(family != "gaussian" & method != "binomial" &
     family != "poisson") stop("'family' must be either 'gaussian', 'binomial'
                               or 'poisson'")


  mf <- model.frame(formula,data=data)

  # Extract outcome data
  y <- as.numeric(model.response(mf))
  n <- length(y)

  # Extract covariates matrix
  D <- as.matrix(model.matrix(attr(mf,"terms"),data=data))

  if(!is.null(hr_re)) {
    # Define indices of random effects
    re_mf <- model.frame(hr_re,data=data)
    names_re <- colnames(re_mf)
    n_re <- ncol(re_mf)

    ID_re <- matrix(NA, nrow = n, ncol = n_re)
    re_unique <- list()
    for(i in 1:n_re) {
      re_unique[[names_re[i]]] <- unique(re_mf[,i])
      ID_re[, i] <- sapply(1:n,
                           function(j) which(re_mf[j,i]==re_unique[[names_re[i]]]))
    }
  } else {
    n_re <- 0
    re_unique <- NULL
    ID_re <- NULL
  }


  # Extract coordinates
  if(is.null(convert_to_crs)) {
    data <- st_transform(data, crs = propose_utm(data))
  } else {
    if(!is.numeric(convert_to_crs)) stop("'convert_to_utm' must be a numeric object")
    data <- st_transform(data, crs = convert_to_crs)
  }
  if(messages) cat("The CRS used is", as.list(st_crs(data))$input, "\n")

  coords_o <- st_coordinates(data)
  coords <- unique(coords_o)

  m <- nrow(coords_o)
  ID_coords <- sapply(1:m, function(i)
               which(coords_o[i,1]==coords[,1] &
                     coords_o[i,2]==coords[,2]))
  s_unique <- unique(ID_coords)

  if(scale_to_km) {
    coords_o <- coords_o/1000
    coords <- coords/1000
    if(messages) cat("Distances between locations are computed in kilometers \n")
  } else {
    if(messages) cat("Distances between locations are computed in meters \n")
  }

  if(geo_re=="GP+NUG" & family=="gaussian" & nrow(coords)==nrow(coords_o)) {
    stop("When slectig geo_re='GP+NUG' for family='gaussian' and with no repeated
         observations at multiple locations, the model can only be fitted by
         fixing the variance of the nugget effect through the argument 'fix_tau2'")
  }

  if(geo_re=="GP") {
    fix_tau2 <- 0
  }

  if(!is.null(fix_tau2)) {
    if(geo_re=="GP" & fix_tau2!=0) warning("Because geo_re='GP', the argument passed to
                             'fix_tau2' is ignored and 'fix_tau2=0'")
    if(class(fix_tau2)!="numeric") stop("'fix_tau2' must be numeric")
    if(fix_tau2 < 0) stop("'fix_tau2' must be a non-negative number")
    if(!is.null(start_cov_pars)) {
      if(n_re>0) {
        if(length(start_cov_pars)!=3+n_re) stop("...")
      } else {
        if(length(start_cov_pars)!=3) stop("When passing a value 'fix_tau2' withuot additional non-structured random effects,
                                        'start_cov_pars' should contain three values
                                         corresponding to the starting values for
                                         sigma2, phi and the variance of the measurement error")
      }
    } else {
      if(geo_re=="GP") {
        if(n_re>0) {
          start_cov_pars <- c(1, quantile(dist(coords),0.25), 1, rep(1, n_re))
        } else {
          start_cov_pars <- c(1, quantile(dist(coords), 0.25), 1)
        }
      } else {
        if(n_re>0) {
          start_cov_pars <- c(1, quantile(dist(coords),0.25), 1, 1, rep(1, n_re))
        } else {
          start_cov_pars <- c(1, quantile(dist(coords),0.25), 1, 1)
        }
      }
    }

  } else {
    if(!is.null(start_cov_pars)) {
      if(n_re>0) {
        if(length(start_cov_pars)!=4+n_re) stop("...")
      } else {
        if(length(start_cov_pars)!=4) stop("...")
      }
    } else {
      if(n_re>0) {
        start_cov_pars <- c(1, quantile(dist(coords),0.25), 1, 1, rep(1, n_re))
      } else {
        start_cov_pars <- c(1, quantile(dist(coords),0.25), 1, 1)
      }
    }
  }

  if(!is.null(start_beta)) {
    if(length(start_beta)!=ncol(D)) stop("The values passed to 'start_beta' do not match
                                  the covariates passed to the 'formula'.")
  } else {
    start_beta <- as.numeric(solve(t(D)%*%D)%*%t(D)%*%y)
  }


  if(family=="gaussian") {
    glgm_lm(y, D, coords, ID_coords, ID_re, s_unique, re_unique,
            fix_tau2, method, start_beta, start_cov_pars)
  } else if(family=="binomial" | family=="poisson") {

  }
  res$call <- match.call()
  return(res)
}
```
