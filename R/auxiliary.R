
##' @title EPSG of the UTM zone
##' @description Suggests EPSG of the UTM in which the majority of data falls in.
##' @param data an object of class \code{sf} containing the variable for which the variogram
##' is to be computed and the coordinates
##' @details The function checks in which UTM zone and empishere the majority of the
##' data fall in and proposes an EPSG.
##'
##' @return an integer indicating the EPSG of the UTM zone.
##'
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##' @importFrom sf st_transform st_coordinates
##' @export
##'
propose_utm <- function(data) {
  if(class(data)[1] != "sf") stop("'data' must be an object of class sf")
  if(is.na(st_crs(data))) stop("the CRS of the data is missing
                               and must be specified; see ?st_crs")

  data <- st_transform(data,crs=4326)
  utm_z <- floor((st_coordinates(data)[,2] + 180) / 6) + 1
  utm_z_u <- unique(utm_z)
  if(length(utm_z_u) > 1) {
    tab_utm <- table(utm_z)
    if(all(diff(tab_utm)==0)) warning("an equal amount of locations falls
                                      in different UTM zones")
    utm_z_u <- as.numeric(names(which.max(tab_utm)))
  }
  ns <- sign(st_coordinates(data)[,1] )
  ns_u <- unique(ns)
  if(length(ns_u) > 1) {
    tab_ns <- table(ns_u)
    if(all(diff(tab_ns)==0)) warning("an equal amount of locations falls
                                      north and south of the Equator")
    ns_u <- as.numeric(names(which.max(tab_ns)))
  }
  if (ns_u==1) {
    out <- as.numeric(paste(326,utm_z_u,sep=""))
  } else if(ns_u==-1) {
    out <- as.numeric(paste(327,utm_z_u,sep=""))
  }
  return(out)
}

##' @title Matern correlation function
##' @description Computes the matern correlation function.
##' @param u 	a vector with values of the distances between pairs of data locations
##' @param phi 	value of the scale parameter \eqn{\phi}.
##' @param kappa 	value of the smoothness parameter \eqn{\kappa}.
##' @param return_sym_matrix a logical value which indicates whether a symmetric
##' correlation matrix should be return. By default \code{return_sym_matrix=FALSE}.
##' @details The Matern correlation function is defined as
##' \deqn{\rho(u; \phi; \kappa) = (2^{\kappa-1})^{-1}(u/\phi)^\kappa K_{\kappa}(u/\phi)}
##' where \eqn{\phi} and \eqn{\kappa} are the scale and smoothness parameters, and
##'  and \eqn{K_{\kappa}(\cdot)} denotes the the modified Bessel function of the third
##'  kind of order \eqn{\kappa}. The parameters \eqn{\phi} and \eqn{\kappa} must be
##'  positive.
##'
##' @return A vector of the same length of \code{u} with the value of the Matern
##' correlation function for the given distances, if \code{return_sym_matrix=FALSE}.
##' If \code{return_sym_matrix=TRUE}, a symmetric correlation matrix is returned.
##'
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##' @importFrom sf st_transform st_coordinates
##' @export
##'
matern_cor <- function(u, phi, kappa, return_sym_matrix = FALSE) {
  if (is.vector(u))
    names(u) <- NULL
  if (is.matrix(u))
    dimnames(u) <- list(NULL, NULL)
  uphi <- u/phi
  uphi <- ifelse(u > 0, (((2^(-(kappa - 1)))/ifelse(0, Inf,
                                                    gamma(kappa))) * (uphi^kappa) * besselK(x = uphi, nu = kappa)),
                 1)
  uphi[u > 600 * phi] <- 0

  if(return_sym_matrix) {
    n <- (1+sqrt(1+8*length(uphi)))/2
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

##' @title First derivative wtr phi
##' @description Computes the first derivative of the Matern correlation
##' function with respect to \eqn{\phi}.
##' @param U 	a vector with values of the distances between pairs of data locations.
##' @param phi 	value of the scale parameter \eqn{\phi}.
##' @param kappa 	value of the smoothness parameter \eqn{\kappa}.
##' @return A matrix with the values of the first derivative of the Matern function
##' with respect to \eqn{\phi} for the given distances
matern.grad.phi <- function(U,phi,kappa) {
  der.phi <- function(u,phi,kappa) {
    u <- u+10e-16
    if(kappa==0.5) {
      out <- (u*exp(-u/phi))/phi^2
    } else {
      out <- ((besselK(u/phi,kappa+1)+besselK(u/phi,kappa-1))*
                phi^(-kappa-2)*u^(kappa+1))/(2^kappa*gamma(kappa))-
        (kappa*2^(1-kappa)*besselK(u/phi,kappa)*phi^(-kappa-1)*
           u^kappa)/gamma(kappa)
    }
    out
  }

  n <- attr(U,"Size")
  grad.phi.mat <- matrix(NA,nrow=n,ncol=n)
  ind <- lower.tri(grad.phi.mat)
  grad.phi <- der.phi(as.numeric(U),phi,kappa)
  grad.phi.mat[ind] <-  grad.phi
  grad.phi.mat <- t(grad.phi.mat)
  grad.phi.mat[ind] <-  grad.phi
  diag(grad.phi.mat) <- rep(der.phi(0,phi,kappa),n)
  grad.phi.mat
}

##' @title Second derivative wtr phi
##' @description Computes the second derivative of the Matern correlation
##' function with respect to \eqn{\phi}.
##' @param U 	a vector with values of the distances between pairs of data locations.
##' @param phi 	value of the scale parameter \eqn{\phi}.
##' @param kappa 	value of the smoothness parameter \eqn{\kappa}.
##' @return A matrix with the values of the second derivative of the Matern function
##' with respect to \eqn{\phi} for the given distances
matern.hessian.phi <- function(U,phi,kappa) {
  der2.phi <- function(u,phi,kappa) {
    u <- u+10e-16
    if(kappa==0.5) {
      out <- (u*(u-2*phi)*exp(-u/phi))/phi^4
    } else {
      bk <- besselK(u/phi,kappa)
      bk.p1 <- besselK(u/phi,kappa+1)
      bk.p2 <- besselK(u/phi,kappa+2)
      bk.m1 <- besselK(u/phi,kappa-1)
      bk.m2 <- besselK(u/phi,kappa-2)
      out <- (2^(-kappa-1)*phi^(-kappa-4)*u^kappa*(bk.p2*u^2+2*bk*u^2+
                                                     bk.m2*u^2-4*kappa*bk.p1*phi*u-4*
                                                     bk.p1*phi*u-4*kappa*bk.m1*phi*u-4*bk.m1*phi*u+
                                                     4*kappa^2*bk*phi^2+4*kappa*bk*phi^2))/(gamma(kappa))
    }
    out
  }
  n <- attr(U,"Size")
  hess.phi.mat <- matrix(NA,nrow=n,ncol=n)
  ind <- lower.tri(hess.phi.mat)
  hess.phi <- der2.phi(as.numeric(U),phi,kappa)
  hess.phi.mat[ind] <-  hess.phi
  hess.phi.mat <- t(hess.phi.mat)
  hess.phi.mat[ind] <-  hess.phi
  diag(hess.phi.mat) <- rep(der2.phi(0,phi,kappa),n)
  hess.phi.mat
}

##' @export
##'
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

##' @export
##'
re <- function (...)
{
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


##' @title Extracting the parameter estimates from a model fit
##' @description \code{coef} method for the class "RiskMap" that extracts the
##' maximum likelihood estimates from the model fits obtained from the function
##' \code{\link{glgpm}}
##' @param object an object of class "RiskMap" obatained as result of a call to \code{\link{glgpm}}
##' @return A vector of the maximum likelihood estimates
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @method coef RiskMap
##' @export
##'
coef.RiskMap <- function(object) {
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

##' @title Summarizing model fits
##' @description \code{summary} method for the class "RiskMap" that computes the standard errors and p-values of likelihood-based model fits.
##' @param object an object of class "RiskMap" obatained as result of a call to \code{\link{glgpm}}
##' @return A list with the following components
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @method summary RiskMap
##' @export
##'

summary.RiskMap <- function(object, conf_level = 0.95) {

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

##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @method print summary.RiskMap
##' @export
print.summary.RiskMap <- function(x) {
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


##' @title Create table in Latex
##' @description Convert a model fit into an \code{xtable} object,
##' which can then be printed as a LaTeX or HTML table.
##' @param object an object of class "RiskMap" obatained as result of a call to \code{\link{glgpm}}
##' @param ... 	additional arguments to be passed to \code{\link{xtable}}.
##'
##' @return For most \code{xtable} methods, an object of class "xtable" which inherits the \code{data.frame}
##' class and contains several additional attributes specifying the table formatting options
##'
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##' @importFrom xtable xtable
##' @export
##'
to_table <- function(object, ...) {
  summary_out <- summary(object)
  tab <- rbind(summary_out$reg_coef[,1:3], summary_out$sp, summary_out$ranef,
               summary_out$me)
  out <- xtable(x = tab,...)
  return(out)
}
