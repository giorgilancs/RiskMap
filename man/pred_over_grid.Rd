% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_prediction.R
\name{pred_over_grid}
\alias{pred_over_grid}
\title{Prediction of the random effects components and covariates effects over a spatial grid using a fitted generalized linear Gaussian process model}
\usage{
pred_over_grid(
  object,
  grid_pred,
  predictors = NULL,
  re_predictors = NULL,
  pred_cov_offset = NULL,
  control_sim = set_control_sim(),
  type = "marginal",
  messages = TRUE
)
}
\arguments{
\item{object}{A RiskMap object obtained from the `glgpm` function.}

\item{grid_pred}{An object of class 'sfc', representing the spatial grid over which predictions
are to be made. Must be in the same coordinate reference system (CRS) as the
object passed to 'object'.}

\item{predictors}{Optional. A data frame containing predictor variables used for prediction.}

\item{re_predictors}{Optional. A data frame containing predictors for unstructured random effects,
if applicable.}

\item{pred_cov_offset}{Optional. A numeric vector specifying covariate offsets at prediction locations.}

\item{control_sim}{Control parameters for MCMC sampling. Must be an object of class "mcmc.RiskMap" as returned by \code{\link{set_control_sim}}.}

\item{type}{Type of prediction. "marginal" for marginal predictions, "joint" for joint predictions.}

\item{messages}{Logical. If TRUE, display progress messages. Default is TRUE.}
}
\value{
An object of class 'RiskMap.pred.re' containing predicted values, uncertainty estimates,
        and additional information.
}
\description{
This function computes predictions over a spatial grid using a fitted model
obtained from the \code{\link{glgpm}} function. It provides point predictions and uncertainty
estimates for the specified locations for each component of the model separately: the spatial random effects;
the unstructured random effects (if included); and the covariates effects.
}
\author{
Emanuele Giorgi \email{e.giorgi@lancaster.ac.uk}

Claudio Fronterre \email{c.fronterr@lancaster.ac.uk}
}
