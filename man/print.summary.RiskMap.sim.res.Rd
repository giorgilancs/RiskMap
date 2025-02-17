% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_prediction.R
\name{print.summary.RiskMap.sim.res}
\alias{print.summary.RiskMap.sim.res}
\title{Print Simulation Results}
\usage{
\method{print}{summary.RiskMap.sim.res}(x, ...)
}
\arguments{
\item{x}{An object of class `RiskMap.sim.res`, as returned by `assess_sim`.}

\item{...}{Additional arguments (not used).}
}
\value{
Invisibly returns `x`.
}
\description{
Prints a concise summary of simulation results from a `RiskMap.sim.res` object, including average metrics by category and a summary of overall correct classification (CC).
}
