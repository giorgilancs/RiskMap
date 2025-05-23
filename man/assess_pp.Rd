% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_prediction.R
\name{assess_pp}
\alias{assess_pp}
\title{Compute scoring rules using spatial cross-validation}
\usage{
assess_pp(
  object,
  keep_par_fixed = TRUE,
  iter = 1,
  fold = NULL,
  n_size = NULL,
  control_sim = set_control_sim(),
  method,
  min_dist = NULL,
  plot_fold = TRUE,
  messages = TRUE,
  which_metric = c("AnPIT", "CRPS", "SCRPS"),
  ...
)
}
\arguments{
\item{object}{A list of `RiskMap` objects, each representing a model fitted with `glgpm`.}

\item{keep_par_fixed}{Logical; if `TRUE`, parameters are kept fixed across folds, otherwise the model is re-estimated for each fold.}

\item{iter}{Integer; number of times to repeat the cross-validation.}

\item{fold}{Integer; number of folds for cross-validation (required if `method = "cluster"`).}

\item{n_size}{Optional; the size of the test set, required if `method = "regularized"`.}

\item{control_sim}{Control settings for simulation, an output from `set_control_sim`.}

\item{method}{Character; either `"cluster"` or `"regularized"` for the cross-validation method. The `"cluster"` method uses
spatial clustering as implemented by the \code{\link{spatial_clustering_cv}} function from the `spatialEco` package, while the `"regularized"` method
selects a subsample of the dataset by imposing a minimum distance, set by the `min_dist` argument, for a randomly selected
subset of locations.}

\item{min_dist}{Optional; minimum distance for regularized subsampling (required if `method = "regularized"`).}

\item{plot_fold}{Logical; if `TRUE`, plots each fold's test set.}

\item{messages}{Logical; if `TRUE`, displays progress messages.}

\item{which_metric}{Character; either `"CRPS"` or `"SCRPS"` to specify the scoring rule.}

\item{...}{Additional arguments passed to clustering or subsampling functions.}
}
\value{
A list of class `RiskMap.spatial.cv`, containing:
  - `test_set`: A list containing all the test sets used for the validation in
  'sf' class.
  - `score`: A list with either `CRPS` or `SCRPS` scores for each fold, depending on `which_metric`.
  - `refit`: A list of re-fitted models for each fold if `keep_par_fixed = FALSE`.
}
\description{
This function calculates the predictive accuracy of a spatial model fitted to a `RiskMap` object using cross-validation.
It allows model scoring based on specified metrics, with options for two cross-validation methods: spatial clustering and regularized subsampling.
Users can choose between continuous ranked probability score (CRPS) and scaled CRPS (SCRPS) as scoring metrics to evaluate predictive quality.
For each data fold, the function can either refit the model or use fixed parameters, enabling flexible model validation and evaluation.
Additionally, it can generate plots of test sets across folds, providing visual insights into the spatial cross-validation structure.
}
\references{
Bolin, D., & Wallin, J. (2023). Local scale invariance and robustness of proper scoring rules. *Statistical Science*, 38(1), 140–159. \doi{10.1214/22-STS864}.
}
\seealso{
\code{\link{spatial_clustering_cv}}, \code{\link{subsample.distance}}
}
\author{
Emanuele Giorgi
}
