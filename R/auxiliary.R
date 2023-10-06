
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


gp <- function (..., kappa = 0.5, nugget = NULL, crs = NULL) {
  vars <- as.list(substitute(list(...)))[-1]
  d <- length(vars)
  term <- NULL

  if(length(nugget) > 0) {
    if(!is.numeric(nugget) |
       (is.numeric(nugget) & nugget <0)) stop("when 'nugget' is not NULL, this must be a positive
                                 real number")
  }
  if(length(crs)>0) {
    if(!is.numeric(crs) |
       (is.numeric(crs) &
        (crs%%1!=0 | crs <0))) stop("'crs' must be a positive integer number")
  }

  if (d == 0) {
    term <- "sf"
  } else {
    if (d > 0) {
      for (i in 1:d) {
        term[i] <- deparse(vars[[i]], backtick = TRUE, width.cutoff = 500)
      }
      if(is.null(crs)) {
        warning("'crs' is set to 4326 (long/lat)")
        crs <- 4326
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
  ret <- list(term = term, kappa = kappa, nugget = nugget, dim = d, crs=crs,
              label = label)
  class(ret) <- "gp.spec"
  ret
}

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
