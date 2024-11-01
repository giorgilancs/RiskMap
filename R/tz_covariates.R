##' @title Covariates Dataset for Malaria Prediction in Tanzania
##'
##' @description This dataset provides covariates over a 10 by 10 km regular grid covering Tanzania. It is intended to be used together with the `tz_malaria` dataset for spatial prediction of malaria prevalence.
##'
##' \describe{
##'   \item{Population}{Population density in the area (in thousands).}
##'   \item{ITN}{Percentage of households with at least one insecticide-treated net (ITN).}
##'   \item{EVI}{Enhanced Vegetation Index, indicating vegetation density.}
##'   \item{Temperature}{Average temperature in degrees Celsius.}
##'   \item{NTL}{Nighttime light intensity, indicating urbanization and infrastructure.}
##'   \item{Precipitation}{Total precipitation in millimeters.}
##'   \item{utm_x}{UTM (Universal Transverse Mercator) x-coordinate of the grid point.}
##'   \item{utm_y}{UTM (Universal Transverse Mercator) y-coordinate of the grid point.}
##' }
##'
##' @name tz_covariates
##' @docType data
##' @usage data(tz_covariates)
##' @keywords datasets
##' @format A data frame with 8740 observations of 8 variables.
##' The CRS of the UTM coordinates is 32736.
##' @source Giorgi E, Fronterrè C, Macharia PM, Alegana VA, Snow RW, Diggle PJ. 2021 Model building and assessment of the impact of covariates for disease prevalence mapping in low-resource settings: to explain and to predict. J. R. Soc. Interface 18: 20210104. \href{https://doi.org/10.1098/rsif.2021.0104}{https://doi.org/10.1098/rsif.2021.0104}
NULL
