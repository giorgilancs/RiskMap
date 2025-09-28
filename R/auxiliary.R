#' @importFrom stats as.formula binomial coef complete.cases
#' @importFrom stats glm median model.frame model.matrix
#' @importFrom stats model.response na.fail na.omit nlminb pnorm
#' @importFrom stats poisson printCoefmat qnorm reformulate rnorm
#' @importFrom stats runif sd step terms terms.formula update

##' @title Convex Hull of an sf Object
##'
##' @description Computes the convex hull of an `sf` object, returning the boundaries of the smallest polygon that can enclose all geometries in the input.
##'
##' @param sf_object An `sf` data frame object containing geometries.
##'
##' @return An `sf` object representing the convex hull of the input geometries.
##'
##' @details The convex hull is the smallest convex polygon that encloses all points in the input `sf` object. This function computes the convex hull by first uniting all geometries in the input using `st_union()`, and then applying `st_convex_hull()` to obtain the polygonal boundary. The result is returned as an `sf` object containing the convex hull geometry.
##'
##' @seealso \code{\link[sf]{st_convex_hull}}, \code{\link[sf]{st_union}}
##'
##' @examples
##' library(sf)
##'
##' # Create example sf object
##' points <- st_sfc(st_point(c(0,0)), st_point(c(1,1)), st_point(c(2,2)), st_point(c(0,2)))
##' sf_points <- st_sf(geometry = points)
##'
##' # Calculate the convex hull
##' convex_hull_result <- convex_hull_sf(sf_points)
##'
##' # Plot the result
##' plot(sf_points, col = 'blue', pch = 19)
##' plot(convex_hull_result, add = TRUE, border = 'red')
##' @importFrom sf st_geometry st_convex_hull st_sf st_union
##' @export
convex_hull_sf <- function(sf_object) {
  # Check if the input is an sf object
  if (!inherits(sf_object, "sf")) {
    stop("`sf_object` must be an sf object")
  }

  # Get the geometry from the sf object
  geometry <- st_geometry(sf_object)

  # Compute the convex hull
  convex_hull <- st_convex_hull(st_union(geometry))

  # Return the convex hull as an sf object
  return(st_sf(geometry = convex_hull))
}


##' @title EPSG of the UTM Zone
##' @description Suggests the EPSG code for the UTM zone where the majority of the data falls.
##' @param data An object of class \code{sf} containing the coordinates.
##' @details The function determines the UTM zone and hemisphere where the majority of the data points are located and proposes the corresponding EPSG code.
##' @return An integer indicating the EPSG code of the UTM zone.
##' @author
##' Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##' @importFrom sf st_transform st_coordinates
##' @export
propose_utm <- function (data) {
  if (class(data)[1] != "sf")
    stop("'data' must be an object of class sf")
  if (is.na(st_crs(data)))
    stop("the CRS of the data is missing and must be specified; see ?st_crs")

  # Transform to WGS84 (EPSG:4326) to ensure coordinates are in lon/lat
  data <- st_transform(data, crs = 4326)

  # Calculate UTM Zone
  utm_z <- floor((st_coordinates(data)[, 1] + 180)/6) + 1
  utm_z_u <- unique(utm_z)

  if (length(utm_z_u) > 1) {
    tab_utm <- table(utm_z)
    if (all(diff(tab_utm) == 0))
      warning("An equal amount of locations falls in different UTM zones")
    utm_z_u <- as.numeric(names(which.max(tab_utm)))
  }

  # Determine Hemisphere (fixing the latitude check)
  ns <- sign(st_coordinates(data)[, 2])  # Use latitude, not longitude
  ns_u <- unique(ns)

  if (length(ns_u) > 1) {
    tab_ns <- table(ns_u)
    if (all(diff(tab_ns) == 0))
      warning("An equal amount of locations falls north and south of the Equator")
    ns_u <- as.numeric(names(which.max(tab_ns)))
  }

  # Construct EPSG code for UTM zone
  if (ns_u == 1) {
    out <- as.numeric(paste0(326, utm_z_u))  # Northern Hemisphere
  } else if (ns_u == -1) {
    out <- as.numeric(paste0(327, utm_z_u))  # Southern Hemisphere
  }

  return(out)
}


##' @title Matern Correlation Function
##' @description Computes the Matern correlation function.
##' @param u A vector of distances between pairs of data locations.
##' @param phi The scale parameter \eqn{\phi}.
##' @param kappa The smoothness parameter \eqn{\kappa}.
##' @param return_sym_matrix A logical value indicating whether to return a symmetric correlation matrix. Defaults to \code{FALSE}.
##' @details The Matern correlation function is defined as
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##' \deqn{\rho(u; \phi; \kappa) = (2^{\kappa-1})^{-1}(u/\phi)^\kappa K_{\kappa}(u/\phi)}
##' where \eqn{\phi} and \eqn{\kappa} are the scale and smoothness parameters, and \eqn{K_{\kappa}(\cdot)} denotes the modified Bessel function of the third kind of order \eqn{\kappa}. The parameters \eqn{\phi} and \eqn{\kappa} must be positive.
##' @return A vector of the same length as \code{u} with the values of the Matern correlation function for the given distances, if \code{return_sym_matrix=FALSE}. If \code{return_sym_matrix=TRUE}, a symmetric correlation matrix is returned.
##' @importFrom sf st_transform st_coordinates
##' @export
matern_cor <- function(u, phi, kappa, return_sym_matrix = FALSE) {
  if (is.vector(u))
    names(u) <- NULL
  if (is.matrix(u))
    dimnames(u) <- list(NULL, NULL)
  uphi <- u / phi
  uphi <- ifelse(u > 0, (((2^(-(kappa - 1)))/ifelse(0, Inf,
                                                    gamma(kappa))) * (uphi^kappa) * besselK(x = uphi, nu = kappa)),
                 1)
  uphi[u > 600 * phi] <- 0

  if(return_sym_matrix) {
    n <- (1 + sqrt(1 + 8 * length(uphi))) / 2
    varcov <- matrix(NA, n, n)
    varcov[lower.tri(varcov)] <- uphi
    varcov <- t(varcov)
    varcov[lower.tri(varcov)] <- uphi
    diag(varcov) <- 1
    out <- varcov
  } else {
    out <- uphi
  }
  return(out)
}

##' @title First Derivative with Respect to \eqn{\phi}
##' @description Computes the first derivative of the Matern correlation function with respect to \eqn{\phi}.
##' @param U A vector of distances between pairs of data locations.
##' @param phi The scale parameter \eqn{\phi}.
##' @param kappa The smoothness parameter \eqn{\kappa}.
##' @return A matrix with the values of the first derivative of the Matern function with respect to \eqn{\phi} for the given distances.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##' @export
matern.grad.phi <- function(U, phi, kappa) {
  der.phi <- function(u, phi, kappa) {
    u <- u + 10e-16
    if(kappa == 0.5) {
      out <- (u * exp(-u / phi)) / phi^2
    } else {
      out <- ((besselK(u / phi, kappa + 1) + besselK(u / phi, kappa - 1)) *
                phi^(-kappa - 2) * u^(kappa + 1)) / (2^kappa * gamma(kappa)) -
        (kappa * 2^(1 - kappa) * besselK(u / phi, kappa) * phi^(-kappa - 1) *
           u^kappa) / gamma(kappa)
    }
    out
  }

  n <- attr(U, "Size")
  grad.phi.mat <- matrix(NA, nrow = n, ncol = n)
  ind <- lower.tri(grad.phi.mat)
  grad.phi <- der.phi(as.numeric(U), phi, kappa)
  grad.phi.mat[ind] <-  grad.phi
  grad.phi.mat <- t(grad.phi.mat)
  grad.phi.mat[ind] <-  grad.phi
  diag(grad.phi.mat) <- rep(der.phi(0, phi, kappa), n)
  grad.phi.mat
}

##' @title Second Derivative with Respect to \eqn{\phi}
##' @description Computes the second derivative of the Matern correlation function with respect to \eqn{\phi}.
##' @param U A vector of distances between pairs of data locations.
##' @param phi The scale parameter \eqn{\phi}.
##' @param kappa The smoothness parameter \eqn{\kappa}.
##' @return A matrix with the values of the second derivative of the Matern function with respect to \eqn{\phi} for the given distances.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##' @export
matern.hessian.phi <- function(U, phi, kappa) {
  der2.phi <- function(u, phi, kappa) {
    u <- u + 10e-16
    if(kappa == 0.5) {
      out <- (u * (u - 2 * phi) * exp(-u / phi)) / phi^4
    } else {
      bk <- besselK(u / phi, kappa)
      bk.p1 <- besselK(u / phi, kappa + 1)
      bk.p2 <- besselK(u / phi, kappa + 2)
      bk.m1 <- besselK(u / phi, kappa - 1)
      bk.m2 <- besselK(u / phi, kappa - 2)
      out <- (2^(-kappa - 1) * phi^(-kappa - 4) * u^kappa * (bk.p2 * u^2 + 2 * bk * u^2 +
                                                               bk.m2 * u^2 - 4 * kappa * bk.p1 * phi * u - 4 *
                                                               bk.p1 * phi * u - 4 * kappa * bk.m1 * phi * u - 4 * bk.m1 * phi * u +
                                                               4 * kappa^2 * bk * phi^2 + 4 * kappa * bk * phi^2)) / (gamma(kappa))
    }
    out
  }
  n <- attr(U, "Size")
  hess.phi.mat <- matrix(NA, nrow = n, ncol = n)
  ind <- lower.tri(hess.phi.mat)
  hess.phi <- der2.phi(as.numeric(U), phi, kappa)
  hess.phi.mat[ind] <-  hess.phi
  hess.phi.mat <- t(hess.phi.mat)
  hess.phi.mat[ind] <-  hess.phi
  diag(hess.phi.mat) <- rep(der2.phi(0, phi, kappa), n)
  hess.phi.mat
}
##' @title Gaussian Process Model Specification
##' @description Specifies the terms, smoothness, and nugget effect for a Gaussian Process (GP) model.
##' @param ... Variables representing the spatial coordinates or covariates for the GP model.
##' @param kappa The smoothness parameter \eqn{\kappa}. Default is 0.5.
##' @param nugget The nugget effect, which represents the variance of the measurement error. Default is 0. A positive numeric value must be provided if not using the default.
##' @details The function constructs a list that includes the specified terms (spatial coordinates or covariates), the smoothness parameter \eqn{\kappa}, and the nugget effect. This list can be used as a specification for a Gaussian Process model.
##' @return A list of class \code{gp.spec} containing the following elements:
##' \item{term}{A character vector of the specified terms.}
##' \item{kappa}{The smoothness parameter \eqn{\kappa}.}
##' \item{nugget}{The nugget effect.}
##' \item{dim}{The number of specified terms.}
##' \item{label}{A character string representing the full call for the GP model.}
##' @note The nugget effect must be a positive real number if specified.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##' @export
gp <- function (..., kappa = 0.5, nugget = 0) {
  vars <- as.list(substitute(list(...)))[-1]
  d <- length(vars)
  term <- NULL

  if(length(nugget) > 0) {
    if(!is.numeric(nugget) |
       (is.numeric(nugget) & nugget <0)) stop("when 'nugget' is not NULL, this must be a positive
                                 real number")
  }

  if (d == 0) {
    term <- "sf"
  } else {
    if (d > 0) {
      for (i in 1:d) {
        term[i] <- deparse(vars[[i]], backtick = TRUE, width.cutoff = 500)
      }
    }

    for (i in 1:d) term[i] <- attr(terms(reformulate(term[i])),
                                   "term.labels")
  }
  full.call <- paste("gp(", term[1], sep = "")
  if (d > 1)
    for (i in 2:d) full.call <- paste(full.call, ",", term[i],
                                      sep = "")
  label <- gsub("sf", "", paste(full.call, ")", sep = ""))
  ret <- list(term = term, kappa = kappa, nugget = nugget, dim = d,
              label = label)
  class(ret) <- "gp.spec"
  ret
}
##' @title Random Effect Model Specification
##' @description Specifies the terms for a random effect model.
##' @param ... Variables representing the random effects in the model.
##' @details The function constructs a list that includes the specified terms for the random effects. This list can be used as a specification for a random effect model.
##' @return A list of class \code{re.spec} containing the following elements:
##' \item{term}{A character vector of the specified terms.}
##' \item{dim}{The number of specified terms.}
##' \item{label}{A character string representing the full call for the random effect model.}
##' @note At least one variable must be provided as input.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##' @export
re <- function (...) {
  vars <- as.list(substitute(list(...)))[-1]
  d <- length(vars)
  term <- NULL

  if (d == 0) {
    stop("You need to provide at least one variable.")
  } else {
    if (d > 0) {
      for (i in 1:d) {
        term[i] <- deparse(vars[[i]], backtick = TRUE, width.cutoff = 500)
      }
    }
    for (i in 1:d) term[i] <- attr(terms(reformulate(term[i])),
                                   "term.labels")
  }
  full.call <- paste("re(", term[1], sep = "")
  if (d > 1)
    for (i in 2:d) full.call <- paste(full.call, ",", term[i],
                                      sep = "")
  label <- gsub("sf", "", paste(full.call, ")", sep = ""))
  ret <- list(term = term, dim = d, label = label)
  class(ret) <- "re.spec"
  ret
}

interpret.formula <- function(formula) {
  p.env <- environment(formula)
  tf <- terms.formula(formula, specials = c("gp", "re"))
  terms <- attr(tf, "term.labels")
  nt <- length(terms)

  if (attr(tf, "response") > 0) {
    response <- as.character(attr(tf, "variables")[2])
  } else {
    response <- NULL
  }

  gp <- attr(tf, "specials")$gp
  re <- attr(tf, "specials")$re
  off <- attr(tf, "offset")
  vtab <- attr(tf, "factors")

  if (length(gp) > 0) {
    for (i in 1:length(gp)) {
      ind <- (1:nt)[as.logical(vtab[gp[i], ])]
      gp[i] <- ind
    }
  }

  if (length(re) > 0) {
    for (i in 1:length(re)) {
      ind <- (1:nt)[as.logical(vtab[re[i], ])]
      re[i] <- ind
    }
  }

  len.gp <- length(gp)
  len.re <- length(re)
  gp.spec <- eval(parse(text = terms[gp]), envir = p.env)
  re.spec <- eval(parse(text = terms[re]), envir = p.env)

  if (length(off) > 0) {
    offset <- as.character(attr(tf, "variables")[[off[i] + 1]])[2]
  } else {
    offset <- NULL
  }

  if (length(terms[-c(gp, re)]) > 0) {
    pf <- paste(response, "~", paste(terms[-c(gp, re)], collapse = " + "))
  } else if (length(terms[-c(gp, re)]) == 0) {
    pf <- paste(response, "~ 1")
  }

  if (attr(tf, "intercept") == 0) {
    pf <- paste(pf, "-1", sep = "")
  }

  ret <- list(
    pf = as.formula(pf, p.env),
    gp.spec = gp.spec,
    re.spec = re.spec,
    offset = offset,
    response = response
  )
  ret
}


##' @title Extract Parameter Estimates from a "RiskMap" Model Fit
##' @description This \code{coef} method for the "RiskMap" class extracts the
##' maximum likelihood estimates from model fits obtained from the \code{\link{glgpm}} function.
##' @param object An object of class "RiskMap" obtained as a result of a call to \code{\link{glgpm}}.
##' @param ... other parameters.
##' @return A list containing the maximum likelihood estimates:
##' \item{beta}{A vector of coefficient estimates.}
##' \item{sigma2}{The estimate for the variance parameter \eqn{\sigma^2}.}
##' \item{phi}{The estimate for the spatial range parameter \eqn{\phi}.}
##' \item{tau2}{The estimate for the nugget effect parameter \eqn{\tau^2}, if applicable.}
##' \item{sigma2_me}{The estimate for the measurement error variance \eqn{\sigma^2_{me}}, if applicable.}
##' \item{sigma2_re}{A vector of variance estimates for the random effects, if applicable.}
##' @details The function processes the \code{RiskMap} object to extract and name the estimated parameters appropriately, transforming them if necessary.
##' @note This function handles both Gaussian and non-Gaussian families, and accounts for fixed and random effects in the model.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##' @seealso \code{\link{glgpm}}
##' @method coef RiskMap
##' @export
##'
coef.RiskMap <- function(object,...) {
  n_re <- length(object$re)
  if(n_re > 0) {
    re_names <- names(object$re)
  }

  p <- ncol(as.matrix(object$D))
  ind_beta <- 1:p

  if(p==1) {
    object$D <- as.matrix(object$D)
    names(object$estimate)[ind_beta] <- "Intercept"
  } else {
    names(object$estimate)[ind_beta] <- colnames(object$D)
  }
  ind_sigma2 <- p+1
  names(object$estimate)[ind_sigma2] <- "sigma2"
  ind_phi <- p+2
  names(object$estimate)[ind_phi] <- "phi"

  if(is.null(object$fix_tau2)) {
    ind_tau2 <- p+3
    names(object$estimate)[ind_tau2] <- "tau2"
    object$estimate[ind_tau2] <- object$estimate[ind_tau2]+object$estimate[ind_sigma2]
    if(object$family=="gaussian") {
      if(is.null(object$fix_var_me)) {
        ind_sigma2_me <- p+4
        if(n_re>0) ind_sigma2_re <- (p+5):(p+4+n_re)
      } else {
        ind_sigma2_me <- NULL
        if(n_re>0) ind_sigma2_re <- (p+4):(p+3+n_re)
      }
    } else {
      ind_sigma2_me <- NULL
      if(n_re>0) ind_sigma2_re <- (p+4):(p+3+n_re)
    }
  } else {
    ind_tau2 <- NULL
    if(object$family=="gaussian") {
      if(is.null(object$fix_var_me)) {
        ind_sigma2_me <- p+3
        names(object$estimate)[ind_sigma2_me] <- "sigma2_me"
        if(n_re>0) {
          ind_sigma2_re <- (p+4):(p+3+n_re)
        }
      } else {
        ind_sigma2_me <- NULL
        if(n_re>0) {
          ind_sigma2_re <- (p+3):(p+2+n_re)
        }
      }
    } else {
      if(n_re>0) {
        ind_sigma2_re <- (p+3):(p+2+n_re)
      }
    }
  }
  ind_sp <- c(ind_sigma2, ind_phi, ind_tau2)
  object$estimate[ind_sp] <- exp(object$estimate[ind_sp])

  n_p <- length(object$estimate)

  if(n_re > 0) {
    for(i in 1:n_re) {
      names(object$estimate)[ind_sigma2_re[i]] <- paste(re_names[i],"_sigma2_re",sep="")
    }
  }


  res <- list()
  res$beta <- object$estimate[ind_beta]
  res$sigma2 <- as.numeric(object$estimate[ind_sigma2])
  res$phi <- as.numeric(object$estimate[ind_phi])
  if(object$family=="gaussian") {
    if(!is.null(ind_sigma2_me)) res$sigma2_me <- as.numeric(exp(object$estimate[ind_sigma2_me]))
  }
  if(!is.null(ind_tau2)) res$tau2 <- object$estimate[ind_tau2]
  if(n_re>0) res$sigma2_re <- as.numeric(object$estimate[ind_sigma2_re])
  dast_model <- !is.null(object$power_val)

  if(dast_model) {
    if(!is.null(ind_tau2)) {
      if(is.null(object$fix_alpha)) {
        ind_alpha <- p+n_re+4
        ind_gamma <- p+n_re+5
      } else {
        ind_gamma <- p+n_re+4
      }
    } else {
      if(is.null(object$fix_alpha)) {
        ind_alpha <- p+n_re+3
        ind_gamma <- p+n_re+4
      } else {
        ind_gamma <- p+n_re+3
      }
    }
    if(is.null(object$fix_alpha)) res$alpha <- as.numeric(1/(1+exp(-object$estimate[ind_alpha])))
    res$gamma <- as.numeric(exp(object$estimate[ind_gamma]))
  }
  if(object$sst) {
    ind_psi <- length(object$estimate)
    res$psi <- as.numeric(exp(object$estimate[ind_psi]))
  }
  return(res)
}

##' @title Summarize Model Fits
##' @description Provides a \code{summary} method for the "RiskMap" class that computes the standard errors and p-values for likelihood-based model fits.
##' @param object An object of class "RiskMap" obtained as a result of a call to \code{\link{glgpm}}.
##' @param ... other parameters.
##' @param conf_level The confidence level for the intervals (default is 0.95).
##' @return A list containing:
##' \item{reg_coef}{A matrix with the estimates, standard errors, z-values, p-values, and confidence intervals for the regression coefficients.}
##' \item{me}{A matrix with the estimates and confidence intervals for the measurement error variance, if applicable.}
##' \item{sp}{A matrix with the estimates and confidence intervals for the spatial process parameters.}
##' \item{tau2}{The fixed nugget variance, if applicable.}
##' \item{ranef}{A matrix with the estimates and confidence intervals for the random effects variances, if applicable.}
##' \item{conf_level}{The confidence level used for the intervals.}
##' \item{family}{The family of the model (e.g., "gaussian").}
##' \item{kappa}{The kappa parameter of the model.}
##' \item{log.lik}{The log-likelihood of the model fit.}
##' \item{cov_offset_used}{A logical indicating if a covariance offset was used.}
##' \item{aic}{The Akaike Information Criterion (AIC) for the model, if applicable.}
##' @details This function computes the standard errors and p-values for the parameters of a "RiskMap" model, adjusting for the covariance structure if needed.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##' @note Handles both Gaussian and non-Gaussian families, and accounts for fixed and random effects in the model.
##' @seealso \code{\link{glgpm}}, \code{\link{coef.RiskMap}}
##' @method summary RiskMap
##' @export
summary.RiskMap <- function(object, ..., conf_level = 0.95) {

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # --- Link reporting ----------------------------------------------------------
  link_name <- NULL
  inv_expr  <- NULL
  if (!is.null(object$linkf) && is.list(object$linkf) && is.function(object$linkf$inv)) {
    link_name <- object$linkf$name %||% "custom"
    inv_expr  <- tryCatch(
      paste0("Inverse link function = ", paste(deparse(body(object$linkf$inv), width.cutoff = 500L), collapse = " ")),
      error = function(e) "Inverse link function = <user-supplied function>"
    )
  } else {
    # Fallback: canonical links by family
    if (identical(object$family, "poisson")) {
      link_name <- "canonical (log)"
      inv_expr  <- "Inverse link function = exp(x)"
    } else if (identical(object$family, "binomial")) {
      link_name <- "canonical (logit)"
      inv_expr  <- "Inverse link function = 1 / (1 + exp(-x))"
    } else if (identical(object$family, "gaussian")) {
      link_name <- "identity"
      inv_expr  <- "Inverse link function = x"
    }
  }

  # --- existing code -----------------------------------------------------------
  n_re <- length(object$re)
  if (n_re > 0) {
    re_names <- names(object$re)
  }
  alpha <- 1 - conf_level

  p <- ncol(object$D)
  ind_beta <- 1:p

  names(object$estimate)[ind_beta] <- colnames(object$D)
  ind_sigma2 <- p + 1
  names(object$estimate)[ind_sigma2] <- "Spatial process var."
  ind_phi <- p + 2
  names(object$estimate)[ind_phi] <- "Spatial corr. scale"
  dast_model <- !is.null(object$power_val)
  sst <- object$sst

  if (sst) {
    ind_psi <- length(object$estimate)
  }
  if (is.null(object$fix_tau2)) {
    ind_tau2 <- p + 3
    names(object$estimate)[ind_tau2] <- "Variance of the nugget"
    object$estimate[ind_tau2] <- object$estimate[ind_tau2] + object$estimate[ind_sigma2]
    if (object$family == "gaussian") {
      if (is.null(object$fix_var_me)) {
        ind_sigma2_me <- p + 4
      } else {
        ind_sigma2_me <- NULL
      }
      if (n_re > 0) {
        ind_sigma2_re <- (p + 5):(p + 4 + n_re)
      }
    } else {
      ind_sigma2_re <- (p + 4):(p + 3 + n_re)
    }
    if (dast_model) {
      if (is.null(object$fix_alpha)) {
        ind_alpha <- p + n_re + 4
        ind_gamma <- p + n_re + 5
      } else {
        ind_gamma <- p + n_re + 4
      }

    } else {
      ind_alpha <- NULL
      ind_gamma <- NULL
    }
  } else {
    ind_tau2 <- NULL
    if (object$family == "gaussian") {
      if (is.null(object$fix_var_me)) {
        ind_sigma2_me <- p + 3
        names(object$estimate)[ind_sigma2_me] <- "Measuremment error var."
        object$estimate[ind_sigma2_me] <- exp(object$estimate[ind_sigma2_me])
      } else {
        ind_sigma2_me <- NULL
      }
      if (n_re > 0) {
        ind_sigma2_re <- (p + 4):(p + 3 + n_re)
      }
    } else {
      ind_sigma2_re <- (p + 3):(p + 2 + n_re)
    }
    if (is.null(object$fix_alpha)) {
      ind_alpha <- p + n_re + 3
      ind_gamma <- p + n_re + 4
    } else {
      ind_alpha <- NULL
      ind_gamma <- p + n_re + 3
    }
  }
  ind_sp <- c(ind_sigma2, ind_phi, ind_tau2)

  if (dast_model) {
    if (!is.null(object$fix_alpha)) {
      names(object$estimate)[ind_gamma] <- "Scale of the decay (gamma)"
    } else {
      names(object$estimate)[ind_alpha] <- "Drop (alpha)"
      names(object$estimate)[ind_gamma] <- "Scale of the decay (gamma)"
    }
  }
  n_p <- length(object$estimate)
  object$estimate[-c(ind_beta, ind_alpha, ind_gamma)] <-
    exp(object$estimate[-c(ind_beta, ind_alpha, ind_gamma)])

  if (n_re > 0) {
    for (i in 1:n_re) {
      names(object$estimate)[ind_sigma2_re[i]] <- paste(re_names[i], " (random eff. var.)", sep = "")
    }
  }

  J <- diag(1:n_p)
  if (length(ind_tau2) > 0) J[ind_tau2, ind_sigma2] <- 1

  H <- solve(-object$covariance)
  H_new <- t(J) %*% H %*% J

  covariance_new <- solve(-H_new)

  se_par <- sqrt(diag(covariance_new))
  res <- list()

  # Regression coefficients
  zval <- object$estimate[ind_beta] / se_par[ind_beta]
  res$reg_coef <- cbind(
    Estimate     = object$estimate[ind_beta],
    "Lower limit" = object$estimate[ind_beta] - se_par[ind_beta] * qnorm(1 - alpha / 2),
    "Upper limit" = object$estimate[ind_beta] + se_par[ind_beta] * qnorm(1 - alpha / 2),
    StdErr       = se_par[ind_beta],
    z.value      = zval,
    p.value      = 2 * pnorm(-abs(zval))
  )

  # Measurement error variance (linear model)
  if (object$family == "gaussian") {
    if (is.null(object$fix_var_me)) {
      res$me <- cbind(
        Estimate     = object$estimate[ind_sigma2_me],
        "Lower limit" = exp(log(object$estimate[ind_sigma2_me]) - qnorm(1 - alpha / 2) * se_par[ind_sigma2_me]),
        "Upper limit" = exp(log(object$estimate[ind_sigma2_me]) + qnorm(1 - alpha / 2) * se_par[ind_sigma2_me])
      )
    } else {
      res$me <- object$fix_var_me
    }
  }

  # Spatial process
  res$sp <- cbind(
    Estimate     = object$estimate[ind_sp],
    "Lower limit" = exp(log(object$estimate[ind_sp]) - qnorm(1 - alpha / 2) * se_par[ind_sp]),
    "Upper limit" = exp(log(object$estimate[ind_sp]) + qnorm(1 - alpha / 2) * se_par[ind_sp])
  )

  if (!is.null(object$fix_tau2)) res$tau2 <- object$fix_tau2

  # Random effects
  if (n_re > 0) {
    res$ranef <- cbind(
      Estimate     = object$estimate[ind_sigma2_re],
      "Lower limit" = exp(log(object$estimate[ind_sigma2_re]) - qnorm(1 - alpha / 2) * se_par[ind_sigma2_re]),
      "Upper limit" = exp(log(object$estimate[ind_sigma2_re]) + qnorm(1 - alpha / 2) * se_par[ind_sigma2_re])
    )
  }

  if (dast_model) {
    anti_logit <- function(x) 1 / (1 + exp(-x))

    if (is.null(object$fix_alpha)) {
      est_alpha   <- anti_logit(object$estimate[ind_alpha])
      lower_alpha <- anti_logit(object$estimate[ind_alpha] - qnorm(1 - alpha / 2) * se_par[ind_alpha])
      upper_alpha <- anti_logit(object$estimate[ind_alpha] + qnorm(1 - alpha / 2) * se_par[ind_alpha])
    } else {
      est_alpha <- lower_alpha <- upper_alpha <- NULL
      res$alpha <- object$fix_alpha
    }
    est_gamma   <- exp(object$estimate[ind_gamma])
    lower_gamma <- exp(object$estimate[ind_gamma] - qnorm(1 - alpha / 2) * se_par[ind_gamma])
    upper_gamma <- exp(object$estimate[ind_gamma] + qnorm(1 - alpha / 2) * se_par[ind_gamma])

    res$dast_par <- cbind(
      Estimate     = c(est_alpha, est_gamma),
      "Lower limit" = c(lower_alpha, lower_gamma),
      "Upper limit" = c(upper_alpha, upper_gamma)
    )
    res$power_val <- object$power_val
  }

  if (sst) {
    est_psi   <- exp(object$estimate[ind_psi])
    lower_psi <- exp(object$estimate[ind_psi] - qnorm(1 - alpha / 2) * se_par[ind_psi])
    upper_psi <- exp(object$estimate[ind_psi] + qnorm(1 - alpha / 2) * se_par[ind_psi])
    res$sp <- rbind(res$sp, c(est_psi, lower_psi, upper_psi))
    rownames(res$sp)[3] <- "Temporal corr. scale"
  }

  res$conf_level      <- conf_level
  res$sst             <- sst
  res$family          <- object$family
  res$dast            <- dast_model
  res$kappa           <- object$kappa
  res$log.lik         <- object$log.lik
  res$cov_offset_used <- !(is.null(object$cov_offset) || all(object$cov_offset==0))
  if (object$family == "gaussian") {
    # fix bug: use length(object$estimate), not res$estimate
    res$aic <- 2 * length(object$estimate) - 2 * res$log.lik
  }

  # New: include call and link info for printing
  res$call      <- object$call %||% NULL
  res$link_name <- link_name
  res$invlink_expression <- inv_expr

  class(res) <- "summary.RiskMap"
  return(res)
}


##' @title Print Summary of RiskMap Model
##' @description Provides a \code{print} method for the summary of "RiskMap" objects, detailing the model type, parameter estimates, and other relevant statistics.
##' @param x An object of class "summary.RiskMap".
##' @param ... other parameters.
##' @details This function prints a detailed summary of a fitted "RiskMap" model, including:
##' \itemize{
##'   \item The type of geostatistical model (e.g., Gaussian, Binomial, Poisson).
##'   \item Confidence intervals for parameter estimates.
##'   \item Regression coefficients with their standard errors and p-values.
##'   \item Measurement error variance, if applicable.
##'   \item Spatial process parameters, including the Matern covariance parameters.
##'   \item Variance of the nugget effect, if applicable.
##'   \item Unstructured random effects variances, if applicable.
##'   \item Log-likelihood of the model.
##'   \item Akaike Information Criterion (AIC) for Gaussian models.
##' }
##' @return This function is used for its side effect of printing to the console. It does not return a value.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##' @method print summary.RiskMap
##' @export
print.summary.RiskMap <- function(x, ...) {

  # First line: the call
  if (!is.null(x$call)) {
    cat("Call:\n")
    cat(paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  }

  # Model type
  if (x$family == "gaussian") {
    cat("Linear geostatistical model\n")
  } else if (x$family == "binomial") {
    if (x$dast) {
      cat("Decay adjusted spatio temporal model\n")
    } else {
      cat("Binomial geostatistical linear model\n")
    }
  } else if (x$family == "poisson") {
    cat("Poisson geostatistical linear model\n")
  }

  # Link report
  if (!is.null(x$link_name) || !is.null(x$invlink_expression)) {
    if (!is.null(x$link_name)) cat("Link:", x$link_name, "\n")
    if (!is.null(x$invlink_expression)) cat(x$invlink_expression, "\n \n")
  }

  cat("'Lower limit' and 'Upper limit' are the limits of the ",
      x$conf_level * 100, "% confidence level intervals\n", sep = "")

  cat("\nRegression coefficients\n")
  printCoefmat(x$reg_coef, P.values = TRUE, has.Pvalue = TRUE)
  if (x$cov_offset_used) cat("Offset included into the linear predictor\n")

  if (x$family == "gaussian") {
    if (length(x$me) > 1) {
      cat("\n")
      printCoefmat(x$me, Pvalues = FALSE)
    } else {
      cat("\nMeasurement error var. fixed at ", x$me, "\n", sep = "")
    }
  }

  if (!x$sst) {
    cat("\nSpatial Gaussian process\n")
    cat("Matern covariance parameters (kappa=", x$kappa, ")\n", sep = "")
  } else {
    cat("\nSpatio temporal Gaussian process\n")
    cat("Separable spatio temporal correlation function:\n")
    cat("Spatial Matern covariance parameters (kappa=", x$kappa, ")  x  Exponential time correlation function\n", sep = "")
  }
  printCoefmat(x$sp, Pvalues = FALSE)
  if (!is.null(x$tau2)) cat("Variance of the nugget effect fixed at ", x$tau2, "\n", sep = "")

  if (!is.null(x$dast) && isTRUE(x$dast)) {
    cat("\nMDA impact function\n")
    cat("f(v) = alpha * exp(-(v/gamma)^delta)\n")
    cat("The parameter delta is fixed to ", x$power_val, "\n", sep = "")
    if (!is.null(x$alpha)) {
      cat("The drop (alpha) parameter of the MDA impact function is fixed at ", x$alpha, "\n", sep = "")
    }
    printCoefmat(x$dast_par, Pvalues = FALSE)
  }

  if (!is.null(x$ranef)) {
    cat("\nUnstructured random effects\n")
    printCoefmat(x$ranef, Pvalues = FALSE)
  }
  cat("\nLog-likelihood: ", x$log.lik, "\n", sep = "")
  if (x$family == "gaussian" && !is.null(x$aic)) cat("AIC: ", x$aic, "\n", sep = "")
}



##' @title Generate LaTeX Tables from RiskMap Model Fits and Validation
##' @description Converts a fitted "RiskMap" model or cross-validation results into an \code{xtable} object, formatted for easy export to LaTeX or HTML.
##' @param object An object of class "RiskMap" resulting from a call to \code{\link{glgpm}}, or a summary object of class "summary.RiskMap.spatial.cv" containing cross-validation results.
##' @param ... Additional arguments to be passed to \code{\link[xtable]{xtable}} for customization.
##' @details This function creates a summary table from a fitted "RiskMap" model or cross-validation results for multiple models, returning it as an \code{xtable} object.
##'
##' When the input is a "RiskMap" model object, the table includes:
##' \itemize{
##'   \item Regression coefficients with their estimates, confidence intervals, and p-values.
##'   \item Parameters for the spatial process.
##'   \item Random effect variances.
##'   \item Measurement error variance, if applicable.
##' }
##'
##' When the input is a cross-validation summary object ("summary.RiskMap.spatial.cv"), the table includes:
##' \itemize{
##'   \item A row for each model being compared.
##'   \item Performance metrics such as CRPS and SCRPS for each model.
##' }
##'
##' The resulting \code{xtable} object can be further customized with additional formatting options and printed as a LaTeX or HTML table for reports or publications.
##' @return An object of class "xtable", which contains the formatted table as a \code{data.frame} and several attributes specifying table formatting options.
##' @importFrom xtable xtable
##' @export
##' @seealso \code{\link{glgpm}}, \code{\link[xtable]{xtable}}, \code{\link{summary.RiskMap.spatial.cv}}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
to_table <- function(object, ...) {
  summary_out <- summary(object)
  if(inherits(summary_out,
               what = "summary.RiskMap", which = FALSE)) {
    tab <- rbind(summary_out$reg_coef[,1:3], summary_out$sp, summary_out$ranef,
                 summary_out$me)
    out <- xtable(x = tab,...)
  } else if (inherits(summary_out,
                      what = "summary.RiskMap.spatial.cv", which = FALSE)) {
    n_models <- nrow(summary_out)
    n_metrics <- ncol(summary_out)
    model_names <- rownames(summary_out)
    metric_names <- toupper(colnames(summary_out))
    tab <- data.frame(Model = model_names)
    for(i in 1:n_metrics) {
      tab[[paste(metric_names[i])]] <- summary_out[,i]
    }
    out <- xtable(x = tab,...)
  }
  return(out)
}

##' @title Compute Unique Coordinate Identifiers
##'
##' @description
##' This function identifies unique coordinates from a `sf` (simple feature) object
##' and assigns an identifier to each coordinate occurrence. It returns a list
##' containing the identifiers for each row and a vector of unique identifiers.
##'
##' @param data_sf An `sf` object containing geometrical data from which coordinates are extracted.
##'
##' @return A list with the following elements:
##' \describe{
##'   \item{ID_coords}{An integer vector where each element corresponds to a row in the input,
##'   indicating the index of the unique coordinate in the full set of unique coordinates.}
##'   \item{s_unique}{An integer vector containing the unique identifiers of all distinct coordinates.}
##' }
##'
##' @details
##' The function extracts the coordinate pairs from the `sf` object and determines the unique
##' coordinates. It then assigns each row in the input data an identifier corresponding
##' to the unique coordinate it matches.
##'
##' @importFrom sf st_coordinates
##' @export
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##'
##'
compute_ID_coords <- function(data_sf) {
  if(!inherits(data_sf,
               what = c("sfc","sf"), which = FALSE)) {
    stop("The object passed to 'grid_pred' must be an object
         of class 'sfc'")
  }
  coords_o <- st_coordinates(data_sf)
  coords <- unique(coords_o)

  m <- nrow(coords_o)
  ID_coords <- sapply(1:m, function(i)
    which(coords_o[i,1]==coords[,1] &
            coords_o[i,2]==coords[,2]))
  out <- list()
  out$ID_coords <- ID_coords
  out$s_unique <- unique(ID_coords)
  return(out)
}


##' @title Summarize Cross-Validation Scores for Spatial RiskMap Models
##'
##' @description This function summarizes cross-validation scores for different spatial models obtained
##' from \code{\link{assess_pp}}.
##'
##' @param object A `RiskMap.spatial.cv` object containing cross-validation scores for each
##'               model, as obtained from \code{\link{assess_pp}}.
##' @param view_all Logical. If `TRUE`, stores the average scores across test sets for each
##'                 model alongside the overall average across all models. Defaults to `TRUE`.
##' @param ... Additional arguments passed to or from other methods.
##'
##' @details
##' The function computes and returns a matrix where rows correspond to models and columns
##' correspond to performance metrics (e.g., CRPS, SCRPS). Scores are weighted by subset sizes
##' to compute averages. Attributes of the returned object include:
##' \itemize{
##'   \item `test_set_means`: A list of average scores for each test set and model.
##'   \item `overall_averages`: Overall averages for each metric across all models.
##'   \item `view_all`: Indicates whether averages across test sets are available for visualization.
##' }
##'
##' @return A matrix of summary scores with models as rows and metrics as columns, with class
##' `"summary.RiskMap.spatial.cv"`.
##'
##' @seealso \code{\link{assess_pp}}
##'
##' @export
##' @method summary RiskMap.spatial.cv
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
summary.RiskMap.spatial.cv <- function(object, view_all = TRUE, ...) {
  model_names <- names(object$model)
  n_models <- length(model_names)

  metric_names <- names(object$model[[1]]$score)
  if (is.null(metric_names)) stop("No metrics of predictive performance were computed when running 'assess_pp'")
  n_metrics <- length(metric_names)

  res <- matrix(NA, ncol = n_metrics, nrow = n_models)
  colnames(res) <- metric_names
  rownames(res) <- model_names

  test_set_means <- list()

  n_subs <- length(object$model[[1]]$score[[1]])
  w <- unlist(lapply(object$model[[1]]$score[[1]], length))

  for (i in 1:n_models) {
    model_scores <- list()
    for (j in 1:n_metrics) {
      score_j <- rep(NA, n_subs)
      for (h in 1:n_subs) {
        score_j[h] <- mean(object$model[[i]]$score[[j]][[h]])
      }
      model_scores[[j]] <- score_j
      res[i, j] <- sum(w * score_j) / sum(w)
    }
    test_set_means[[model_names[i]]] <- model_scores
  }

  overall_averages <- colMeans(res, na.rm = TRUE)

  # Attach additional attributes for printing
  attr(res, "test_set_means") <- test_set_means
  attr(res, "overall_averages") <- overall_averages
  attr(res, "view_all") <- view_all

  class(res) <- "summary.RiskMap.spatial.cv"
  return(res)
}

##' @title Print Summary of RiskMap Spatial Cross-Validation Scores
##'
##' @description This function prints the matrix of cross-validation scores produced by
##' `summary.RiskMap.spatial.cv` in a readable format.
##'
##' @param x An object of class `"summary.RiskMap.spatial.cv"`, typically the output of
##'          `summary.RiskMap.spatial.cv`.
##' @param ... Additional arguments passed to or from other methods.
##'
##' @details
##' This method is primarily used to format and display the summary score matrix,
##' printing it to the console. It provides a clear view of the cross-validation performance
##' metrics across different spatial models.
##'
##' @return This function is used for its side effect of printing to the console. It does not
##'         return a value.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @export
##' @method print summary.RiskMap.spatial.cv
print.summary.RiskMap.spatial.cv <- function(x, ...) {
  # Extract attributes
  test_set_means <- attr(x, "test_set_means")
  overall_averages <- attr(x, "overall_averages")
  view_all <- attr(x, "view_all")

  cat("Summary of Cross-Validation Scores\n")
  cat("----------------------------------\n")

  for (model_name in names(test_set_means)) {
    cat(sprintf("Model: %s\n", model_name))

    if (view_all) {
      # Print scores for each test set
      model_test_set_means <- test_set_means[[model_name]]
      n_test_sets <- length(model_test_set_means[[1]])  # Number of test sets (assumes all metrics have same length)

      for (test_set_idx in seq_len(n_test_sets)) {
        cat(sprintf("  Test Set %d:\n", test_set_idx))
        for (metric_idx in seq_along(model_test_set_means)) {
          metric_name <- colnames(x)[metric_idx]
          test_set_value <- model_test_set_means[[metric_idx]][test_set_idx]
          cat(sprintf("    %s: %.4f\n", metric_name, test_set_value))
        }
      }
    }

    # Print overall average across test sets for the model
    cat("  Overall average across test sets:\n")
    for (metric_idx in seq_along(overall_averages)) {
      metric_name <- colnames(x)[metric_idx]
      overall_avg <- x[model_name, metric_idx]
      cat(sprintf("    %s: %.4f\n", metric_name, overall_avg))
    }
    cat("\n")
  }
}


##' @title Plot Calibration Curves (AnPIT / PIT) from Spatial Cross-Validation
##'
##' @description
##' Produce calibration plots from a \code{RiskMap.spatial.cv} object returned by
##' \code{\link{assess_pp}}.
##' * For Binomial or Poisson models the function visualises the
##'   \emph{Aggregated normalised Probability Integral Transform} (AnPIT)
##'   curves stored in \code{$AnPIT}.
##' * For Gaussian models it detects the list \code{$PIT} and instead plots
##'   the empirical \emph{Probability Integral Transform} curve
##'   (ECDF of PIT values) on the same \eqn{u}-grid.
##'
##' A 45° dashed red line indicates perfect calibration.
##'
##' @param object       A \code{RiskMap.spatial.cv} object.
##' @param mode         One of \code{"average"} (average curve across test sets),
##'                     \code{"single"} (a specific test set),
##'                     or \code{"all"} (every test set separately).
##' @param test_set     Integer; required when \code{mode = "single"}.
##' @param model_name   Optional character string; if supplied,
##'                     only that model is plotted.
##' @param combine_panels Logical; when \code{mode = "average"}, draw
##'                       all models in a single panel (\code{TRUE})
##'                       or one panel per model (\code{FALSE}, default).
##'
##' @return A \pkg{ggplot2} object (single plot) or a \pkg{ggpubr} grid.
##'
##' @importFrom ggplot2 ggplot aes geom_line geom_abline labs theme_minimal guides guide_legend
##' @importFrom ggpubr  ggarrange
##' @importFrom dplyr   filter group_by summarize %>%
##' @importFrom stats    ecdf
##' @export
plot_AnPIT <- function(object,
                       mode = "average",
                       test_set = NULL,
                       model_name = NULL,
                       combine_panels = FALSE) {

  if (!inherits(object, "RiskMap.spatial.cv"))
    stop("`object` must be a 'RiskMap.spatial.cv' produced by assess_pp().")

  all_models <- names(object$model)

  if (!is.null(model_name)) {
    if (!model_name %in% all_models)
      stop("Model name '", model_name, "' not found in `object$model`.")
    all_models <- model_name
  }

  make_df <- function(mname) {
    m <- object$model[[mname]]
    if (!is.null(m$AnPIT)) {
      lapply(seq_along(m$AnPIT), function(j) {
        curve_vals <- m$AnPIT[[j]]
        if (length(curve_vals) == 0) return(NULL)
        data.frame(
          u_val   = seq(0, 1, length.out = length(curve_vals)),
          value   = curve_vals,
          test_set = j,
          model    = mname,
          type     = "AnPIT"
        )
      })
    } else if (!is.null(m$PIT)) {
      u_grid <- seq(0, 1, length.out = 1000)
      lapply(seq_along(m$PIT), function(j) {
        pit_vec <- m$PIT[[j]]
        if (length(pit_vec) == 0) return(NULL)
        data.frame(
          u_val   = u_grid,
          value   = ecdf(pit_vec)(u_grid),
          test_set = j,
          model    = mname,
          type     = "PIT"
        )
      })
    } else {
      NULL
    }
  }

  plot_data <- do.call(rbind, unlist(lapply(all_models, make_df), recursive = FALSE))

  if (is.null(plot_data) || nrow(plot_data) == 0)
    stop("No AnPIT or PIT data available for plotting.")

  y_label <- unique(plot_data$type)
  if (length(y_label) > 1) y_label <- "Calibration curve"

  id_line <- geom_abline(intercept = 0, slope = 1,
                         linetype = "dashed", colour = "red")

  if (mode == "average" && combine_panels) {
    avg <- plot_data %>%
      dplyr::group_by(model, u_val) %>%
      dplyr::summarize(value = mean(value), .groups = "drop")

    return(
      ggplot(avg, aes(u_val, value, colour = model)) +
        geom_line() + id_line +
        labs(title = "Average calibration curves",
             x = "", y = y_label) +
        theme_minimal() +
        guides(colour = guide_legend(title = "Model"))
    )
  }

  build_plot <- function(df, title_suffix = "") {
    ggplot(df, aes(u_val, value,
                   colour = if (mode == "all") as.factor(test_set) else NULL)) +
      geom_line() + id_line +
      labs(title = title_suffix, x = "", y = unique(df$type)) +
      theme_minimal() +
      guides(colour = guide_legend(title = "Test set"))
  }

  plots <- list()
  for (mname in all_models) {
    df_model <- dplyr::filter(plot_data, model == mname)

    p <- switch(mode,
                average = {
                  avg <- df_model %>%
                    dplyr::group_by(u_val) %>%
                    dplyr::summarize(value = mean(value), .groups = "drop")
                  avg$type <- unique(df_model$type)
                  build_plot(avg, paste("Model", mname, ": average"))
                },
                single  = {
                  if (is.null(test_set))
                    stop("Provide `test_set` when mode = 'single'.")
                  df_ts <- dplyr::filter(df_model, test_set == test_set)
                  if (nrow(df_ts) == 0)
                    stop("No data for test_set ", test_set, " in model ", mname)
                  build_plot(df_ts,
                             paste("Model", mname, "- test set", test_set))
                },
                all     = build_plot(df_model,
                                     paste("Model", mname, "- all test sets")),
                stop("Invalid `mode`. Use 'average', 'single' or 'all'.")
    )

    plots[[mname]] <- p
  }

  if (length(plots) == 1) {
    plots[[1]]
  } else {
    ncol <- ifelse(length(plots) == 2, 2, 2)
    nrow <- ceiling(length(plots) / ncol)
    ggpubr::ggarrange(plotlist = plots, ncol = ncol, nrow = nrow)
  }
}


##' @title Plot Spatial Scores for a Specific Model and Metric
##'
##' @description This function visualizes spatial scores for a specified model and metric.
##' It combines test set data, handles duplicate locations by averaging scores,
##' and creates a customizable map using ggplot2.
##'
##' @param object A list containing test sets and model scores. The structure should include
##'   `object$test_set` (list of sf objects) and `object$model[[which_model]]$score[[which_score]]`.
##' @param which_score A string specifying the score to visualize. Must match a score computed in the model.
##' @param which_model A string specifying the model whose scores to visualize.
##' @param ... Additional arguments to customize ggplot, such as `scale_color_gradient` or `scale_color_manual`.
##' @return A ggplot object visualizing the spatial distribution of the specified score.
##' @export
plot_score <- function(object, which_score, which_model, ...) {
  geometry <- NULL
  score <- NULL

  # Check if "which_score" exists
  if (!which_score %in% names(object$model[[which_model]]$score)) {
    stop(paste("Error: The score", shQuote(which_score), "was not computed for model", shQuote(which_model)))
  }

  # Extract the test sets and number of test sets
  test_sets <- object$test_set
  n_test <- length(test_sets)

  # Combine the data and add the score variable
  data_full <- st_as_sf(test_sets[[1]])
  data_full$score <- object$model[[which_model]]$score[[which_score]][[1]]

  if (n_test > 1) {
    for (i in 1:n_test) {
      test_sets[[i]]$score <- object$model[[which_model]]$score[[which_score]][[i]]
      data_full <- rbind(data_full, test_sets[[i]])
    }
  }

  # Check for duplicate locations and average the score
  data_full <- data_full %>%
    mutate(geom_id = st_as_text(geometry)) %>%
    group_by(geom_id) %>%
    summarize(score = mean(score, na.rm = TRUE),
              geometry = first(geometry), .groups = "drop") %>%
    st_as_sf()


  # Create the base plot
  out <- ggplot(data = data_full) +
    geom_sf(aes(color = score), size = 2) +
    ggtitle(paste("Visualizing", which_score, "for model", which_model)) +
    theme_minimal()


  return(out)
}

##' @title Plot the estimated MDA impact function
##'
##' @description
##' Generate a plot of the estimated impact of mass drug administration (MDA)
##' on infection prevalence, based on a fitted decay-adjusted spatio-temporal (DAST) model.
##' The function simulates draws from the posterior distribution of model parameters,
##' propagates them through the MDA effect function, and produces uncertainty bands
##' around the estimated impact curve.
##'
##' @param object A fitted DAST model object, returned by \code{\link{dast}}.
##' @param mda_history Specification of the MDA schedule. This can be either:
##'   \itemize{
##'     \item A numeric vector of event times (integers starting at 0, e.g. \code{c(0,1,2,6)}),
##'     \item OR a 0/1 indicator vector on the yearly grid (e.g. \code{c(1,1,1,0,0,0,1)}),
##'     where position \code{i} corresponds to year \code{i-1}.
##'   }
##'   If omitted, the default is a single MDA at time 0.
##' @param n_sim Number of posterior draws used for uncertainty quantification (default: 1000).
##' @param x_min Minimum value for the x-axis (default: \code{1e-6}).
##' @param x_max Maximum value for the x-axis (default: \code{10}).
##' @param conf_level Confidence level for the pointwise uncertainty interval (default: 0.95).
##' @param lower_f Optional lower bound for the y-axis. If not provided, computed from the data.
##' @param upper_f Optional upper bound for the y-axis. If not provided, computed from the data.
##' @param mc_cores Number of CPU cores to use for parallel simulation. Default is 1 (serial).
##' @param parallel_backend Parallelisation backend to use. Options are \code{"none"} (default),
##'   \code{"fork"} (Unix-like systems), or \code{"psock"} (cross-platform).
##' @param ... Additional arguments (currently unused).
##'
##' @details
##' The time axis is assumed to start at 0 and increase in integer steps of 1 year.
##' The argument \code{mda_history} allows the user to specify when MDAs occurred either
##' by listing the years directly or by giving a binary indicator on the yearly grid.
##' The function then evaluates the cumulative relative reduction
##' \eqn{1 - \mathrm{effect}(t)} at a dense grid of time points between \code{x_min}
##' and \code{x_max}, using the fitted parameters from the supplied DAST model.
##'
##' @return
##' A \code{ggplot2} object showing the median estimated MDA impact function
##' and the pointwise uncertainty band at the chosen confidence level.
##'
##' @importFrom ggplot2 coord_cartesian geom_ribbon geom_line
##' @export
plot_mda <- function(object,
                     mda_history  = NULL,   # numeric event times (integers, starting at 0) OR 0/1 vector on yearly grid
                     n_sim        = 1000,
                     x_min        = 1e-6,
                     x_max        = 10,
                     conf_level   = 0.95,
                     lower_f      = NULL,
                     upper_f      = NULL,
                     mc_cores     = 1,      # safe default: serial
                     parallel_backend = c("none","fork","psock"),
                     ...) {

  parallel_backend <- match.arg(parallel_backend)

  # --- Time axis for evaluation (smooth curve), but MDA grid is YEARLY from 0 ---
  stopifnot(is.numeric(x_min), is.numeric(x_max), x_max > x_min)
  survey_times <- seq(x_min, x_max, length.out = 200)  # evaluation grid for the curve
  n_t <- length(survey_times)

  # --- MDA schedule (assumed on yearly integer grid: 0,1,2,...) ---
  if (is.null(mda_history)) {
    # Default: single MDA at time 0
    mda_times <- 0
  } else if (is.numeric(mda_history) && all(mda_history %in% c(0,1))) {
    # Indicator vector on the yearly grid: position i -> time (i-1)
    if (length(mda_history) == 0L) {
      mda_times <- numeric(0)
    } else {
      mda_times <- which(mda_history == 1) - 1
    }
  } else if (is.numeric(mda_history)) {
    # Event times given directly (assumed integers on the yearly grid starting at 0)
    mda_times <- sort(unique(as.numeric(mda_history)))
  } else {
    stop("`mda_history` must be numeric: either integer event times (0,1,2,...) or a 0/1 vector on that yearly grid.")
  }

  # --- Extract params from fitted object ---
  par_hat   <- coef(object)
  n_par     <- length(object$estimate)
  power_val <- object$power_val

  if (is.null(par_hat$alpha)) {
    ind_dast    <- n_par
    par_dast    <- log(par_hat$gamma)
    alpha_fixed <- object$fix_alpha
    has_alpha   <- FALSE
  } else {
    ind_dast    <- (n_par - 1):n_par
    par_dast    <- c(log(par_hat$alpha / (1 - par_hat$alpha)),
                     log(par_hat$gamma))
    alpha_fixed <- NA_real_
    has_alpha   <- TRUE
  }

  Sigma_par       <- as.matrix(object$covariance[ind_dast, ind_dast])
  Sigma_par_sroot <- t(chol(Sigma_par))
  par_hat_sim <- t(vapply(
    X   = seq_len(n_sim),
    FUN = function(i) par_dast + Sigma_par_sroot %*% stats::rnorm(length(ind_dast)),
    FUN.VALUE = numeric(length(ind_dast))
  ))

  alphas <- if (has_alpha) plogis(par_hat_sim[, 1]) else rep(alpha_fixed, n_sim)
  gammas <- if (has_alpha) exp(par_hat_sim[, 2]) else exp(par_hat_sim[, 1])

  # --- Build intervention & simulate through compute_mda_effect ---
  if (length(mda_times) == 0L) {
    effects_mat <- matrix(0, nrow = n_t, ncol = n_sim)
  } else {
    # One column per MDA event; one row per evaluation time
    intervention_mat <- matrix(1, nrow = n_t, ncol = length(mda_times))

    one_sim <- function(j) {
      eff <- compute_mda_effect(
        survey_times_data = survey_times,
        mda_times         = mda_times,
        intervention      = intervention_mat,
        alpha             = alphas[j],
        gamma             = gammas[j],
        kappa             = power_val
      )
      1 - eff  # cumulative relative reduction
    }

    # ---- EXECUTION BACKENDS ----
    if (parallel_backend == "none" || mc_cores <= 1L) {
      eff_list <- lapply(seq_len(n_sim), one_sim)

    } else if (parallel_backend == "fork" && .Platform$OS.type == "unix") {
      mc_cores <- as.integer(max(1L, min(mc_cores, parallel::detectCores(logical = TRUE) - 1L, n_sim)))
      eff_list <- parallel::mclapply(seq_len(n_sim), one_sim,
                                     mc.cores = mc_cores, mc.preschedule = TRUE)

    } else if (parallel_backend == "psock") {
      mc_cores <- as.integer(max(1L, min(mc_cores, parallel::detectCores(logical = TRUE), n_sim)))
      cl <- parallel::makeCluster(mc_cores, type = "PSOCK")
      on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
      parallel::clusterExport(cl,
                              varlist = c("survey_times","mda_times","intervention_mat","alphas","gammas","power_val","one_sim"),
                              envir = environment())
      eff_list <- parallel::parLapply(cl, seq_len(n_sim), one_sim)

    } else {
      eff_list <- lapply(seq_len(n_sim), one_sim)
    }

    effects_mat <- do.call(cbind, eff_list)
  }

  # --- Summaries & limits that respect x_min/x_max ---
  alpha_q <- (1 - conf_level) / 2
  med   <- apply(effects_mat, 1, stats::median,   na.rm = TRUE)
  lower <- apply(effects_mat, 1, stats::quantile, probs = alpha_q, na.rm = TRUE)
  upper <- apply(effects_mat, 1, stats::quantile, probs = 1 - alpha_q, na.rm = TRUE)

  in_view <- survey_times >= x_min & survey_times <= x_max
  if (!any(in_view)) in_view <- rep(TRUE, length(survey_times))

  if (is.null(lower_f)) lower_f <- min(lower[in_view], na.rm = TRUE)
  if (is.null(upper_f)) upper_f <- max(upper[in_view], na.rm = TRUE)

  plot_data <- data.frame(
    time   = survey_times,
    median = med,
    lower  = lower,
    upper  = upper
  )

  ggplot2::ggplot(plot_data, ggplot2::aes(x = time)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(y = median), color = "black", linewidth = 1) +
    ggplot2::labs(
      x = "Years since baseline",
      y = "Relative reduction in prevalence",
      title = "MDA Impact Over Time"
    ) +
    ggplot2::coord_cartesian(xlim = c(x_min, x_max), ylim = c(lower_f, upper_f)) +
    ggplot2::theme_minimal()
}
