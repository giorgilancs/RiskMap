#' @importFrom stats as.formula binomial coef complete.cases
#' @importFrom stats glm median model.frame model.matrix
#' @importFrom stats model.response na.fail na.omit nlminb pnorm
#' @importFrom stats poisson printCoefmat qnorm reformulate rnorm
#' @importFrom stats runif sd step terms terms.formula update

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
propose_utm <- function(data) {
  if(class(data)[1] != "sf") stop("'data' must be an object of class sf")
  if(is.na(st_crs(data))) stop("the CRS of the data is missing and must be specified; see ?st_crs")

  data <- st_transform(data, crs = 4326)
  utm_z <- floor((st_coordinates(data)[,2] + 180) / 6) + 1
  utm_z_u <- unique(utm_z)
  if(length(utm_z_u) > 1) {
    tab_utm <- table(utm_z)
    if(all(diff(tab_utm) == 0)) warning("an equal amount of locations falls in different UTM zones")
    utm_z_u <- as.numeric(names(which.max(tab_utm)))
  }
  ns <- sign(st_coordinates(data)[,1])
  ns_u <- unique(ns)
  if(length(ns_u) > 1) {
    tab_ns <- table(ns_u)
    if(all(diff(tab_ns) == 0)) warning("an equal amount of locations falls north and south of the Equator")
    ns_u <- as.numeric(names(which.max(tab_ns)))
  }
  if (ns_u == 1) {
    out <- as.numeric(paste(326, utm_z_u, sep = ""))
  } else if(ns_u == -1) {
    out <- as.numeric(paste(327, utm_z_u, sep = ""))
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
  if (length(gp) > 0)
    for (i in 1:length(gp)) {
      ind <- (1:nt)[as.logical(vtab[gp[i],])]
      gp[i] <- ind
    }
  if (length(re) > 0)
    for (i in 1:length(re)) {
      ind <- (1:nt)[as.logical(vtab[re[i],])]
      re[i] <- ind
    }
  len.gp <- length(gp)
  len.re <- length(re)
  ns <- len.gp + len.re
  gp.spec <- eval(parse(text = terms[gp]), envir = p.env)
  re.spec <- eval(parse(text = terms[re]), envir = p.env)
  if(length(terms[-c(gp, re)])>0) {
    pf <- paste(response, "~", paste(terms[-c(gp, re)], collapse = " + "))
  } else if (length(terms[-c(gp, re)])==0) {
    pf <- paste(response, "~ 1")
  }
  if (attr(tf, "intercept") == 0) {
    pf <- paste(pf, "-1", sep = "")
  }

  fake.formula <- pf
  fake.formula <- as.formula(fake.formula, p.env)
  ret <- list(pf = as.formula(pf, p.env),
              gp.spec = gp.spec,
              re.spec = re.spec,
              response = response)
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

  p <- ncol(object$D)
  ind_beta <- 1:p

  names(object$estimate)[ind_beta] <- colnames(object$D)
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

  n_p <- length(object$estimate)
  object$estimate[-ind_beta] <- exp(object$estimate[-ind_beta])

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
    if(!is.null(ind_sigma2_me)) res$sigma2_me <- as.numeric(object$estimate[ind_sigma2_me])
  }
  if(!is.null(ind_tau2)) res$tau2 <- object$estimate[ind_tau2]
  if(n_re>0) res$sigma2_re <- as.numeric(object$estimate[ind_sigma2_re])
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

  n_re <- length(object$re)
  if(n_re > 0) {
    re_names <- names(object$re)
  }
  alpha <- 1-conf_level

  p <- ncol(object$D)
  ind_beta <- 1:p

  names(object$estimate)[ind_beta] <- colnames(object$D)
  ind_sigma2 <- p+1
  names(object$estimate)[ind_sigma2] <- "Spatial process var."
  ind_phi <- p+2
  names(object$estimate)[ind_phi] <- "Spatial corr. scale"

  if(is.null(object$fix_tau2)) {
    ind_tau2 <- p+3
    names(object$estimate)[ind_tau2] <- "Variance of the nugget"
    object$estimate[ind_tau2] <- object$estimate[ind_tau2]+object$estimate[ind_sigma2]
    if(object$family=="gaussian") {
      if(is.null(object$fix_var_me)) {
        ind_sigma2_me <- p+4
      } else {
        ind_sigma2_me <- NULL
      }
      if(n_re>0) {
        ind_sigma2_re <- (p+5):(p+4+n_re)
      }
    } else {
      ind_sigma2_re <- (p+4):(p+3+n_re)
    }
  } else {
    ind_tau2 <- NULL
    if(object$family=="gaussian") {
      if(is.null(object$fix_var_me)) {
        ind_sigma2_me <- p+3
        names(object$estimate)[ind_sigma2_me] <- "Measuremment error var."
      } else {
        ind_sigma2_me <- NULL
      }
      if(n_re>0) {
        ind_sigma2_re <- (p+4):(p+3+n_re)
      }
    } else {
      ind_sigma2_re <- (p+3):(p+2+n_re)
    }
  }
  ind_sp <- c(ind_sigma2, ind_phi, ind_tau2)

  n_p <- length(object$estimate)
  object$estimate[-ind_beta] <- exp(object$estimate[-ind_beta])

  if(n_re > 0) {
    for(i in 1:n_re) {
      names(object$estimate)[ind_sigma2_re[i]] <- paste(re_names[i]," (random eff. var.)",sep="")
    }
  }

  J <- diag(1:n_p)
  if(length(ind_tau2)>0) J[ind_tau2,ind_sigma2] <- 1

  H <- solve(-object$covariance)
  H_new <- t(J)%*%H%*%J

  covariance_new <- solve(-H_new)

  se_par <- sqrt(diag(covariance_new))
  res <- list()
  # Reg. coefficeints
  zval <- object$estimate[ind_beta]/se_par[ind_beta]
  res$reg_coef <- cbind(Estimate = object$estimate[ind_beta],
                        'Lower limit' = object$estimate[ind_beta]-se_par[ind_beta]*qnorm(1-alpha/2),
                        'Upper limit' = object$estimate[ind_beta]+se_par[ind_beta]*qnorm(1-alpha/2),
                        StdErr = se_par[ind_beta],
                        z.value = zval, p.value = 2 * pnorm(-abs(zval)))

  # Measurement error variance (linear model)
  if(object$call$family=="gaussian") {
    if(is.null(object$fix_var_me)) {
      res$me <- cbind(Estimate = object$estimate[ind_sigma2_me],
                      'Lower limit' = exp(log(object$estimate[ind_sigma2_me])-qnorm(1-alpha/2)*se_par[ind_sigma2_me]),
                      'Upper limit' = exp(log(object$estimate[ind_sigma2_me])+qnorm(1-alpha/2)*se_par[ind_sigma2_me]))
    } else {
      res$me <- object$fix_var_me
    }
  }

  # Spatial process
  res$sp <- cbind(Estimate = object$estimate[ind_sp],
                  'Lower limit' = exp(log(object$estimate[ind_sp])-qnorm(1-alpha/2)*se_par[ind_sp]),
                  'Upper limit' = exp(log(object$estimate[ind_sp])+qnorm(1-alpha/2)*se_par[ind_sp]))

  if(!is.null(object$fix_tau2)) res$tau2 <- object$fix_tau2

  # Random effects
  if(n_re>0) {
    res$ranef <- cbind(Estimate = object$estimate[ind_sigma2_re],
                       'Lower limit' = exp(log(object$estimate[ind_sigma2_re])-qnorm(1-alpha/2)*se_par[ind_sigma2_re]),
                       'Upper limit' = exp(log(object$estimate[ind_sigma2_re])+qnorm(1-alpha/2)*se_par[ind_sigma2_re]))
  }

  res$conf_level <- conf_level
  res$family <- object$family
  res$kappa <- object$kappa
  res$log.lik <- object$log.lik
  res$cov_offset_used <- !is.null(object$cov_offset)
  if(object$family=="gaussian") res$aic <- 2*length(res$estimate)-2*res$log.lik

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
  if(x$family=="gaussian") {
    cat("Linear geostatsitical model \n")
  } else if(x$family=="binomial") {
    cat("Binomial geostatistical linear model \n")
  } else if(x$family=="poisson") {
    cat("Poisson geostatistical linear model \n")
  }

  cat("'Lower limit' and 'Upper limit' are the limits of the ",
      x$conf_level*100,
      "% confidence level intervals \n", sep="")


  cat("\n Regression coefficients \n")
  printCoefmat(x$reg_coef,P.values=TRUE,has.Pvalue=TRUE)
  if(x$cov_offset_used) cat("Offset included into the linear predictor \n")

  if(x$family=="gaussian") {
    if(length(x$me)>1) {
      cat("\n ")
      printCoefmat(x$me, Pvalues = FALSE)
    } else {
      cat("\n Measurement error var. fixed at", x$me,"\n")
    }
  }

  cat("\n Spatial Guassian process \n")
  cat("Matern covariance parameters (kappa=",x$kappa,") \n",sep="")
  printCoefmat(x$sp, Pvalues = FALSE)
  if(!is.null(x$tau2)) cat("Variance of the nugget effect fixed at",x$tau2,"\n")


  if(!is.null(x$ranef)) {
    cat("\n Unstructured random effects \n")
    printCoefmat(x$ranef, Pvalues = FALSE)
  }
  cat("\n Log-likelihood: ",x$log.lik,"\n",sep="")
  if(x$family=="gaussian") cat("\n AIC: ",x$aic,"\n \n",sep="")

}


##' @title Create LaTeX Table from Model Fit
##' @description Converts a "RiskMap" model fit into an \code{xtable} object, which can then be printed as a LaTeX or HTML table.
##' @param object An object of class "RiskMap" obtained as a result of a call to \code{\link{glgpm}}.
##' @param ... Additional arguments to be passed to \code{\link[xtable]{xtable}}.
##' @details This function takes a fitted "RiskMap" model and converts it into an \code{xtable} object. The resulting table includes:
##' \itemize{
##'   \item Regression coefficients with their estimates, confidence intervals, and p-values.
##'   \item Spatial process parameters.
##'   \item Random effects variances.
##'   \item Measurement error variance, if applicable.
##' }
##' The \code{xtable} object can be customized further using additional arguments and then printed as a LaTeX or HTML table.
##' @return An object of class "xtable" which inherits the \code{data.frame} class and contains several additional attributes specifying the table formatting options.
##' @importFrom xtable xtable
##' @export
##' @seealso \code{\link{glgpm}}, \code{\link[xtable]{xtable}}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##'
to_table <- function(object, ...) {
  summary_out <- summary(object)
  tab <- rbind(summary_out$reg_coef[,1:3], summary_out$sp, summary_out$ranef,
               summary_out$me)
  out <- xtable(x = tab,...)
  return(out)
}
