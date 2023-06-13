
##' @title Monte Carlo Maximum Likelihood estimation for the binomial logistic model
##' @description This function performs Monte Carlo maximum likelihood (MCML) estimation for the geostatistical binomial logistic model.
##' @param formula an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
##' @param units.m an object of class \code{\link{formula}} indicating the binomial denominators in the data.
##' @param coords an object of class \code{\link{formula}} indicating the spatial coordinates in the data.
##' @param times an object of class \code{\link{formula}} indicating the times in the data, used in the spatio-temporal model.
##' @param data a data frame containing the variables in the model.
##' @param ID_coords vector of ID values for the unique set of spatial coordinates obtained from \code{\link{create.ID_coords}}. These must be provided if, for example, spatial random effects are defined at household level but some of the covariates are at individual level. \bold{Warning}: the household coordinates must all be distinct otherwise see \code{jitterDupCoords}. Default is \code{NULL}.
##' @param par0 parameters of the importance sampling distribution: these should be given in the following order \code{c(beta,sigma2,phi,tau2)}, where \code{beta} are the regression coefficients, \code{sigma2} is the variance of the Gaussian process, \code{phi} is the scale parameter of the spatial correlation and \code{tau2} is the variance of the nugget effect (if included in the model).
##' @param control.mcmc output from \code{\link{control.mcmc.MCML}}.
##' @param kappa fixed value for the shape parameter of the Matern covariance function.
##' @param kappa.t fixed value for the shape parameter of the Matern covariance function in the separable double-Matern spatio-temporal model.
##' @param sst.model a character value that specifies the spatio-temporal correlation function.
##' \itemize{
##' \item \code{sst.model="DM"} separable double-Matern.
##' \item \code{sst.model="GN1"} separable correlation functions. Temporal correlation: \eqn{f(x) = 1/(1+x/\psi)}; Spatial correaltion: Matern function.
##' }
##' Deafault is \code{sst.model=NULL}, which is used when a purely spatial model is fitted.
##' @param fixed.rel.nugget fixed value for the relative variance of the nugget effect; \code{fixed.rel.nugget=NULL} if this should be included in the estimation. Default is \code{fixed.rel.nugget=NULL}.
##' @param start_cov_pars a vector of length two with elements corresponding to the starting values of \code{phi} and the relative variance of the nugget effect \code{nu2}, respectively, that are used in the optimization algorithm. If \code{nu2} is fixed through \code{fixed.rel.nugget}, then \code{start_cov_pars} represents the starting value for \code{phi} only.
##' @param method method of optimization. If \code{method="BFGS"} then the \code{\link{maxBFGS}} function is used; otherwise \code{method="nlminb"} to use the \code{\link{nlminb}} function. Default is \code{method="BFGS"}.
##' @param low.rank logical; if \code{low.rank=TRUE} a low-rank approximation of the Gaussian spatial process is used when fitting the model. Default is \code{low.rank=FALSE}.
##' @param SPDE logical; if \code{SPDE=TRUE} the SPDE approximation for the Gaussian spatial model is used. Default is \code{SPDE=FALSE}.
##' @param knots if \code{low.rank=TRUE}, \code{knots} is a matrix of spatial knots that are used in the low-rank approximation. Default is \code{knots=NULL}.
##' @param mesh an object obtained as result of a call to the function \code{inla.mesh.2d}.
##' @param messages logical; if \code{messages=TRUE} then status messages are printed on the screen (or output device) while the function is running. Default is \code{messages=TRUE}.
##' @param plot.correlogram logical; if \code{plot.correlogram=TRUE} the autocorrelation plot of the samples of the random effect is displayed after completion of conditional simulation. Default is \code{plot.correlogram=TRUE}.
##' @details
##' This function performs parameter estimation for a geostatistical binomial logistic model. Conditionally on a zero-mean stationary Gaussian process \eqn{S(x)} and mutually independent zero-mean Gaussian variables \eqn{Z} with variance \code{tau2}, the observations \code{y} are generated from a binomial distribution with probability \eqn{p} and binomial denominators \code{units.m}. A canonical logistic link is used, thus the linear predictor assumes the form
##' \deqn{\log(p/(1-p)) = d'\beta + S(x) + Z,}
##' where \eqn{d} is a vector of covariates with associated regression coefficients \eqn{\beta}. The Gaussian process \eqn{S(x)} has isotropic Matern covariance function (see \code{matern}) with variance \code{sigma2}, scale parameter \code{phi} and shape parameter \code{kappa}.
##' In the \code{binomial.logistic.MCML} function, the shape parameter is treated as fixed. The relative variance of the nugget effect, \code{nu2=tau2/sigma2}, can also be fixed through the argument \code{fixed.rel.nugget}; if \code{fixed.rel.nugget=NULL}, then the relative variance of the nugget effect is also included in the estimation.
##'
##' \bold{Monte Carlo Maximum likelihood.}
##' The Monte Carlo maximum likelihood method uses conditional simulation from the distribution of the random effect \eqn{T(x) = d(x)'\beta+S(x)+Z} given the data \code{y}, in order to approximate the high-dimensiional intractable integral given by the likelihood function. The resulting approximation of the likelihood is then maximized by a numerical optimization algorithm which uses analytic epression for computation of the gradient vector and Hessian matrix. The functions used for numerical optimization are \code{\link{maxBFGS}} (\code{method="BFGS"}), from the \pkg{maxLik} package, and \code{\link{nlminb}} (\code{method="nlminb"}).
##'
##' \bold{Using a two-level model to include household-level and individual-level information.}
##' When analysing data from household sruveys, some of the avilable information information might be at household-level (e.g. material of house, temperature) and some at individual-level (e.g. age, gender). In this case, the Gaussian spatial process \eqn{S(x)} and the nugget effect \eqn{Z} are defined at hosuehold-level in order to account for extra-binomial variation between and within households, respectively.
##'
##' \bold{Low-rank approximation.}
##' In the case of very large spatial data-sets, a low-rank approximation of the Gaussian spatial process \eqn{S(x)} might be computationally beneficial. Let \eqn{(x_{1},\dots,x_{m})} and \eqn{(t_{1},\dots,t_{m})} denote the set of sampling locations and a grid of spatial knots covering the area of interest, respectively. Then \eqn{S(x)} is approximated as \eqn{\sum_{i=1}^m K(\|x-t_{i}\|; \phi, \kappa)U_{i}}, where \eqn{U_{i}} are zero-mean mutually independent Gaussian variables with variance \code{sigma2} and \eqn{K(.;\phi, \kappa)} is the isotropic Matern kernel (see \code{\link{matern.kernel}}). Since the resulting approximation is no longer a stationary process (but only approximately), the parameter \code{sigma2} is then multiplied by a factor \code{constant.sigma2} so as to obtain a value that is closer to the actual variance of \eqn{S(x)}.
##' @return An object of class "PrevMap".
##' The function \code{\link{summary.PrevMap}} is used to print a summary of the fitted model.
##' The object is a list with the following components:
##' @return \code{estimate}: estimates of the model parameters; use the function \code{\link{coef.PrevMap}} to obtain estimates of covariance parameters on the original scale.
##' @return \code{covariance}: covariance matrix of the MCML estimates.
##' @return \code{log.lik}: maximum value of the log-likelihood.
##' @return \code{y}: binomial observations.
##' @return \code{units.m}: binomial denominators.
##' @return \code{D}: matrix of covariates.
##' @return \code{coords}: matrix of the observed sampling locations.
##' @return \code{method}: method of optimization used.
##' @return \code{ID_coords}: set of ID values defined through the argument \code{ID_coords}.
##' @return \code{kappa}: fixed value of the shape parameter of the Matern function.
##' @return \code{kappa.t}: fixed value for the shape parameter of the Matern covariance function in the separable double-Matern spatio-temporal model.
##' @return \code{knots}: matrix of the spatial knots used in the low-rank approximation.
##' @return \code{mesh}: the mesh used in the SPDE approximation.
##' @return \code{const.sigma2}: adjustment factor for \code{sigma2} in the low-rank approximation.
##' @return \code{h}: vector of the values of the tuning parameter at each iteration of the Langevin-Hastings MCMC algorithm; see \code{\link{Laplace.sampling}}, or \code{\link{Laplace.sampling.lr}} if a low-rank approximation is used.
##' @return \code{samples}: matrix of the random effects samples from the importance sampling distribution used to approximate the likelihood function.
##' @return \code{fixed.rel.nugget}: fixed value for the relative variance of the nugget effect.
##' @return \code{call}: the matched call.
##' @seealso \code{\link{Laplace.sampling}}, \code{\link{Laplace.sampling.lr}}, \code{\link{summary.PrevMap}}, \code{\link{coef.PrevMap}}, \code{matern}, \code{\link{matern.kernel}},  \code{\link{control.mcmc.MCML}}, \code{\link{create.ID_coords}}.
##' @references Diggle, P.J., Giorgi, E. (2019). \emph{Model-based Geostatistics for Global Public Health.} CRC/Chapman & Hall.
##' @references Giorgi, E., Diggle, P.J. (2017). \emph{PrevMap: an R package for prevalence mapping.} Journal of Statistical Software. 78(8), 1-29. doi: 10.18637/jss.v078.i08
##' @references Christensen, O. F. (2004). \emph{Monte carlo maximum likelihood in model-based geostatistics.} Journal of Computational and Graphical Statistics 13, 702-718.
##' @references Higdon, D. (1998). \emph{A process-convolution approach to modeling temperatures in the North Atlantic Ocean.} Environmental and Ecological Statistics 5, 173-190.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Peter J. Diggle \email{p.diggle@@lancaster.ac.uk}
##' @export
glgm <- function(formula,
                 m_offset = NULL,
                 geo_re = "GP",
                 hr_re = NULL,
                 data,
                 family,
                 convert_to_crs = NULL,
                 scale_to_km = TRUE,
                 control_MCMC = NULL,
                 kappa = 0.5,
                 start_beta = NULL,
                 start_cov_pars = NULL,
                 start_vars.re = NULL,
                 method="BFGS",
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

  if(method != "BFGS" & method != "nlminb") stop("'method' must be either 'BFGS' or 'nlminb'.")

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
      if(length(start_cov_pars)!=3) stop("When passing a value 'fix_tau2',
                                        'start_cov_pars' should contain three values
                                         corresponding to the starting values for
                                         sigma2, phi and the variance of the measurement error")
    } else {
      start_cov_pars <- c(1, quantile(dist(coords),0.25), 1)
    }

  } else {
    if(!is.null(start_cov_pars)) {
      if(length(start_cov_pars)!=4) stop("The argument 'start_cov_pars' should contain four values
                                         corresponding to the starting values for
                                         sigma2, phi, tau2 and the variance of the measurement error")
    } else {
      start_cov_pars <- c(1, quantile(dist(coords),0.25), 1, 1)
    }
  }

  if(!is.null(start_beta)) {
    if(length(start_beta)!=ncol(D)) stop("The values passed to 'start_beta' do not match
                                  the covariates passed to the 'formula'.")
  } else {
    start_beta <- as.numeric(solve(t(D)%*%D)%*%t(D)%*%y)
  }


  if(family=="gaussian") {
    glgm_lm(y, D, coords, ID_coords, ID_re, re_unique,
            fix_tau2, method, start_beta, start_cov_pars)
  } else if(family=="binomial" | family=="poisson") {

  }
  res$call <- match.call()
  return(res)
}


##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @importFrom maxLik maxBFGS
glgm_lm <- function(y, D, coords, ID_coords, ID_re, re_unique,
                    fix_tau2, method, start_beta, start_cov_pars) {

  m <- length(y)
  p <- ncol(D)
  U <- dist(coords)
  if(is.null(ID_re)) {
    n_re <- 0
  } else {
    n_re <- ncol(ID_re)
  }


  ID_g <- as.matrix(cbind(ID_coords, ID_re))

  n.coords <- list()
  for(i in 1:(n_re+1)) {
    n.coords[[i]] <- as.numeric(tapply(ID_g[,i],ID_g[,i], length))
  }


  log.lik <- function(par) {
    beta <- par[1:p]
    sigma2 <- exp(par[p+1])
    phi <- exp(par[p+2])
    if(!is.null(fix_tau2)) {
      nu2.1 <- fix_tau2/sigma2
      nu2.2 <- exp(par[p+3])
    } else {
      nu2.1 <- exp(par[p+3])
      nu2.2 <- exp(par[p+4])
    }

    if(n_re>0) {
      nu2_re <- exp(par[(p+4+1):(p+4+n_re)])
    }

    mu <- as.numeric(D%*%beta)
    diff.y <- y-mu
    diff.y.tilde <- list()

    for(i in 1:(n_re+1)) {
      diff.y.tilde[[i]] <- as.numeric(tapply(diff.y, ID_g[,i], sum))
    }

    R <- matern_cor(U,phi = phi, kappa=kappa,
                    return_sym_matrix = TRUE)
    diag(R) <- diag(R)+nu2.1
    R.inv <- solve(R)
    Omega <- R.inv
    diag(Omega) <- (diag(Omega)+n.coords[[1]]/nu2.2)
    Omega <- Omega*(nu2.2^2)
    Omega.inv <- solve(Omega)

    q.f <- as.numeric(sum(diff.y^2)/nu2.2-
                        t(diff.y.tilde[[1]])%*%Omega.inv%*%diff.y.tilde[[1]])

    q.f_re <- 0
    log.det_re <- 0
    for(i in 1:n_re) {
      o.inv_re <- (1/nu2_re[i]+n.coords[[i+1]]/nu2.2)*(nu2.2^2)
      q.f_re <- q.f_re +  sum(diff.y.tilde[[i+1]]^2/o.inv_re)
      log.det_re <- log.det_re
    }

    M.det <- t(R*n.coords[[1]]/nu2.2)
    for(i in 1:n_re) {
      diag(M.det) <- diag(M.det)+ nu2_re[i]*n.coords[[i+1]]/nu2.2
    }
    diag(M.det) <- diag(M.det)+1
    log.det.tot <- m*log(sigma2)+m*log(nu2.2)+
      as.numeric(determinant(M.det)$modulus)
    out <- -0.5*(log.det.tot+(q.f+q.f_re)/sigma2)
    return(out)
  }

  grad.log.lik <- function(par) {
    beta <- par[1:p]
    sigma2 <- exp(par[p+1])
    phi <- exp(par[p+2])
    if(!is.null(fix_tau2)) {
      nu2.1 <-  fix_tau2/sigma2
      nu2.2 <- exp(par[p+3])
    } else {
      nu2.1 <- exp(par[p+3])
      nu2.2 <- exp(par[p+4])
    }
    R <- matern_cor(U,phi = phi, kappa=kappa,
                    return_sym_matrix = TRUE)
    diag(R) <- diag(R)+nu2.1
    R.inv <- solve(R)
    Omega.star <- R.inv
    diag(Omega.star) <- (diag(Omega.star)+n.coords/nu2.2)
    Omega.star.inv <- solve(Omega.star)


    mu <- as.numeric(D%*%beta)
    diff.y <- y-mu
    diff.y.tilde <- as.numeric(tapply(diff.y,ID_coords,sum))

    v <- as.numeric(Omega.star.inv%*%diff.y.tilde)
    q.f1 <- sum(diff.y^2)
    q.f2 <- as.numeric(t(diff.y.tilde)%*%v)
    q.f <- q.f1/nu2.2-q.f2/(nu2.2^2)

    grad.beta <- as.numeric(t(D)%*%(diff.y/nu2.2-
                                      v[ID_coords]/(nu2.2^2)))/sigma2
    grad.sigma2 <- -0.5*(m-q.f/sigma2)

    M.det <- t(R*n.coords/nu2.2)
    diag(M.det) <- diag(M.det)+1
    M.det.inv <- solve(M.det)

    R1.phi <- matern.grad.phi(U,phi,kappa)*phi
    Omega.star.R.inv <- Omega.star.inv%*%R.inv
    der.Omega.star.inv.phi <- Omega.star.R.inv%*%R1.phi%*%t(Omega.star.R.inv)
    der.M.phi <- t(R1.phi*n.coords/nu2.2)
    M1.trace.phi <- M.det.inv*t(der.M.phi)
    trace1.phi <- sum(M1.trace.phi)
    v.beta.phi <- as.numeric(der.Omega.star.inv.phi%*%diff.y.tilde)
    q.f.phi <-  as.numeric(t(diff.y.tilde)%*%v.beta.phi)
    grad.phi <- -0.5*(trace1.phi-
                        q.f.phi/((nu2.2^2)*sigma2))

    if(is.null(fix_tau2)) {
      der.Omega.star.inv.nu2.1 <- -Omega.star.R.inv%*%t(Omega.star.R.inv)*nu2.1
      trace1.nu2.1 <- sum(diag(M.det.inv)*n.coords*nu2.1/nu2.2)
      v.beta.nu2.1 <- as.numeric(der.Omega.star.inv.nu2.1%*%
                                   diff.y.tilde)
      q.f.nu2.1 <- as.numeric(t(diff.y.tilde)%*%v.beta.nu2.1)
      grad.nu2.1 <- -0.5*(trace1.nu2.1+
                            q.f.nu2.1/((nu2.2^2)*sigma2))
    }

    der.Omega.star.inv.nu2.2 <- -t(Omega.star.inv*n.coords/nu2.2)%*%
      Omega.star.inv
    der.M.nu2.2 <- -t(R*n.coords/nu2.2)
    M1.trace.nu2.2 <- M.det.inv*t(der.M.nu2.2)
    trace1.nu2.2 <- sum(M1.trace.nu2.2)
    v.beta.nu2.2 <- as.numeric(der.Omega.star.inv.nu2.2%*%
                                 diff.y.tilde)
    q.f.nu2.2 <- as.numeric(t(diff.y.tilde)%*%v.beta.nu2.2)
    grad.nu2.2 <- -0.5*(-q.f1/nu2.2+2*q.f2/(nu2.2^2))/sigma2+
      -0.5*(m+trace1.nu2.2+
              q.f.nu2.2/((nu2.2^2)*sigma2))

    ind.beta <- 1:p
    ind.sigma2 <- p+1
    ind.phi <- p+2
    if(is.null(fix_tau2)) {
      g <- c(grad.beta,grad.sigma2,grad.phi,grad.nu2.1,grad.nu2.2)
    } else {
      g <- c(grad.beta,grad.sigma2,grad.phi,grad.nu2.2)
    }

    return(g)
  }

  DtD <- t(D)%*%D
  Dty <- as.numeric(t(D)%*%y)
  D.tilde <- sapply(1:p,function(i) as.numeric(tapply(D[,i],ID_coords,sum)))
  y.tilde <- tapply(y,ID_coords,sum)

  hessian.log.lik <- function(par) {
    beta <- par[1:p]
    sigma2 <- exp(par[p+1])
    phi <- exp(par[p+2])
    if(!is.null(fix_tau2)) {
      nu2.1 <- fix_tau2/sigma2
      nu2.2 <- exp(par[p+3])
    } else {
      nu2.1 <- exp(par[p+3])
      nu2.2 <- exp(par[p+4])
    }
    R <- matern_cor(U,phi = phi, kappa=kappa,
                    return_sym_matrix = TRUE)
    diag(R) <- diag(R)+nu2.1
    R.inv <- solve(R)
    Omega.star <- R.inv
    diag(Omega.star) <- (diag(Omega.star)+n.coords/nu2.2)
    Omega.star.inv <- solve(Omega.star)

    mu <- as.numeric(D%*%beta)
    diff.y <- y-mu
    diff.y.tilde <- as.numeric(tapply(diff.y,ID_coords,sum))

    v <- as.numeric(Omega.star.inv%*%diff.y.tilde)
    q.f1 <- sum(diff.y^2)
    q.f2 <- as.numeric(t(diff.y.tilde)%*%v)
    q.f <- q.f1/nu2.2-q.f2/(nu2.2^2)

    grad.beta <- as.numeric(t(D)%*%(diff.y/nu2.2-
                                      v[ID_coords]/(nu2.2^2)))/sigma2
    grad.sigma2 <- -0.5*(m-q.f/sigma2)

    M.det <- t(R*n.coords/nu2.2)
    diag(M.det) <- diag(M.det)+1
    M.det.inv <- solve(M.det)

    R1.phi <- matern.grad.phi(U,phi,kappa)*phi
    Omega.star.R.inv <- Omega.star.inv%*%R.inv
    der.Omega.star.inv.phi <- Omega.star.R.inv%*%R1.phi%*%t(Omega.star.R.inv)
    der.M.phi <- t(R1.phi*n.coords/nu2.2)
    M1.trace.phi <- M.det.inv*t(der.M.phi)
    trace1.phi <- sum(M1.trace.phi)
    v.beta.phi <- as.numeric(der.Omega.star.inv.phi%*%diff.y.tilde)
    q.f.phi <-  as.numeric(t(diff.y.tilde)%*%v.beta.phi)
    grad.phi <- -0.5*(trace1.phi-
                        q.f.phi/((nu2.2^2)*sigma2))

    if(is.null(fix_tau2)) {
      der.Omega.star.inv.nu2.1 <- -Omega.star.R.inv%*%t(Omega.star.R.inv)*nu2.1
      trace1.nu2.1 <- sum(diag(M.det.inv)*n.coords*nu2.1/nu2.2)
      v.beta.nu2.1 <- as.numeric(der.Omega.star.inv.nu2.1%*%
                                   diff.y.tilde)
      q.f.nu2.1 <- as.numeric(t(diff.y.tilde)%*%v.beta.nu2.1)
      grad.nu2.1 <- -0.5*(trace1.nu2.1+
                            q.f.nu2.1/((nu2.2^2)*sigma2))
    }

    der.Omega.star.inv.nu2.2 <- -t(Omega.star.inv*n.coords/nu2.2)%*%
      Omega.star.inv
    der.M.nu2.2 <- -t(R*n.coords/nu2.2)
    M1.trace.nu2.2 <- M.det.inv*t(der.M.nu2.2)
    trace1.nu2.2 <- sum(M1.trace.nu2.2)
    v.beta.nu2.2 <- as.numeric(der.Omega.star.inv.nu2.2%*%
                                 diff.y.tilde)
    q.f.nu2.2 <- as.numeric(t(diff.y.tilde)%*%v.beta.nu2.2)
    grad.nu2.2 <- -0.5*(-q.f1/nu2.2+2*q.f2/(nu2.2^2))/sigma2+
      -0.5*(m+trace1.nu2.2+
              q.f.nu2.2/((nu2.2^2)*sigma2))

    ind.beta <- 1:p
    ind.sigma2 <- p+1
    ind.phi <- p+2
    if(is.null(fix_tau2)) {
      g <- c(grad.beta,grad.sigma2,grad.phi,grad.nu2.1,grad.nu2.2)
      H <- matrix(NA,p+4,p+4)
      ind.nu2.1 <- p+3
      ind.nu2.2 <- p+4
    } else {
      g <- c(grad.beta,grad.sigma2,grad.phi,grad.nu2.2)
      H <- matrix(NA,p+3,p+3)
      ind.nu2.2 <- p+3
    }



    H[ind.beta,ind.beta] <- (-DtD/nu2.2+
                               t(D.tilde)%*%Omega.star.inv%*%D.tilde/(nu2.2^2))/sigma2

    H[ind.beta,ind.sigma2] <- H[ind.sigma2,ind.beta] <- -grad.beta

    H[ind.beta,ind.phi] <- H[ind.phi,ind.beta] <- -t(D)%*%v.beta.phi[ID_coords]/((nu2.2^2)*sigma2)

    if(is.null(fix_tau2)) {

      H[ind.beta,ind.nu2.1] <- H[ind.nu2.1,ind.beta] <- t(D)%*%v.beta.nu2.1[ID_coords]/((nu2.2^2)*sigma2)

    }

    H[ind.beta,ind.nu2.2] <-
      H[ind.nu2.2,ind.beta] <- as.numeric(t(D)%*%(-diff.y/nu2.2+2*v[ID_coords]/(nu2.2^2))/sigma2)+
      t(D)%*%v.beta.nu2.2[ID_coords]/((nu2.2^2)*sigma2)

    H[ind.sigma2,ind.sigma2] <- -0.5*q.f/sigma2
    H[ind.sigma2,ind.phi] <- H[ind.phi,ind.sigma2] <- -(grad.phi+0.5*(trace1.phi))
    if(is.null(fix_tau2)) {
      H[ind.sigma2,ind.nu2.1] <- H[ind.nu2.1,ind.sigma2] <- -(grad.nu2.1+0.5*(trace1.nu2.1))
    }
    H[ind.sigma2,ind.nu2.2] <- H[ind.nu2.2,ind.sigma2] <- -(grad.nu2.2+0.5*(trace1.nu2.2+m))

    R2.phi <- matern.hessian.phi(U,phi,kappa)*(phi^2)+
      R1.phi
    V1.phi <- -R.inv%*%R1.phi%*%R.inv
    V2.phi <- R.inv%*%(2*R1.phi%*%R.inv%*%R1.phi-R2.phi)%*%R.inv
    der2.Omega.star.inv.phi <- Omega.star.inv%*%(2*V1.phi%*%
                                                   Omega.star.inv%*%V1.phi-V2.phi)%*%
      Omega.star.inv
    M2.trace.phi <- M.det.inv*((R2.phi*n.coords/nu2.2))
    B2.phi <- M.det.inv%*%der.M.phi

    trace2.phi <- sum(M2.trace.phi)-
      sum(B2.phi*t(B2.phi))
    H[ind.phi,ind.phi] <- -0.5*(trace2.phi-
                                  as.numeric(t(diff.y.tilde)%*%
                                               der2.Omega.star.inv.phi%*%diff.y.tilde)/((nu2.2^2)*sigma2))

    if(is.null(fix_tau2)) {

      V1.nu2.1 <- -R.inv%*%R.inv*nu2.1
      V2.phi.nu2.1 <- nu2.1*R.inv%*%(R1.phi%*%R.inv+R.inv%*%R1.phi)%*%R.inv
      der2.Omega.star.inv.phi.nu2.1 <- Omega.star.inv%*%(
        V1.phi%*%Omega.star.inv%*%V1.nu2.1+
          V1.nu2.1%*%Omega.star.inv%*%V1.phi-
          V2.phi.nu2.1)%*%
        Omega.star.inv
      B2.nu2.1 <- t(t(M.det.inv)*n.coords*nu2.1/nu2.2)
      trace2.phi.nu2.1 <- -sum(B2.phi*t(B2.nu2.1))

      H[ind.phi,ind.nu2.1] <- H[ind.nu2.1,ind.phi] <- -0.5*(trace2.phi.nu2.1-
                                                              as.numeric(t(diff.y.tilde)%*%
                                                                           der2.Omega.star.inv.phi.nu2.1%*%diff.y.tilde)/((nu2.2^2)*sigma2))

    }

    der2.Omega.star.inv.phi.nu2.2 <- Omega.star.inv%*%(
      V1.phi%*%t(-Omega.star.inv*
                   n.coords/nu2.2)+
        (-Omega.star.inv*
           n.coords/nu2.2)%*%V1.phi)%*%
      Omega.star.inv
    B2.nu2.2 <- M.det.inv%*%der.M.nu2.2
    trace2.phi.nu2.2 <- -trace1.phi-sum(B2.phi*t(B2.nu2.2))
    H[ind.phi,ind.nu2.2] <- H[ind.nu2.2,ind.phi] <- -0.5*(trace2.phi.nu2.2+
                                                            -as.numeric(t(diff.y.tilde)%*%
                                                                          der2.Omega.star.inv.phi.nu2.2%*%diff.y.tilde)/((nu2.2^2)*sigma2)+
                                                            +2*q.f.phi/((nu2.2^2)*sigma2))

    if(is.null(fix_tau2)) {
      aux.nu2.1 <- 2*R.inv*(nu2.1^2)
      diag(aux.nu2.1) <- diag(aux.nu2.1)-nu2.1
      V2.nu2.1 <- R.inv%*%aux.nu2.1%*%R.inv
      der2.Omega.star.inv.nu2.1 <- Omega.star.inv%*%(2*V1.nu2.1%*%
                                                       Omega.star.inv%*%V1.nu2.1-V2.nu2.1)%*%
        Omega.star.inv

      trace2.nu2.1 <- trace1.nu2.1-
        sum(B2.nu2.1*t(B2.nu2.1))

      H[ind.nu2.1,ind.nu2.1] <- -0.5*(trace2.nu2.1-
                                        as.numeric(t(diff.y.tilde)%*%
                                                     der2.Omega.star.inv.nu2.1%*%diff.y.tilde)/((nu2.2^2)*sigma2))

      der2.Omega.star.inv.nu2.1.nu2.2 <- Omega.star.inv%*%(
        V1.nu2.1%*%t(-Omega.star.inv*
                       n.coords/nu2.2)+
          (-Omega.star.inv*
             n.coords/nu2.2)%*%V1.nu2.1)%*%
        Omega.star.inv

      trace2.nu2.1.nu2.2 <- -trace1.nu2.1-sum(B2.nu2.1*t(B2.nu2.2))
      H[ind.nu2.1,ind.nu2.2] <- H[ind.nu2.2,ind.nu2.1] <- -0.5*(trace2.nu2.1.nu2.2+
                                                                  -as.numeric(t(diff.y.tilde)%*%
                                                                                der2.Omega.star.inv.nu2.1.nu2.2%*%diff.y.tilde)/((nu2.2^2)*sigma2)+
                                                                  -2*q.f.nu2.1/((nu2.2^2)*sigma2))

    }

    aux.nu2.2 <- 2*t(Omega.star.inv*n.coords)*n.coords/
      (nu2.2^2)
    diag(aux.nu2.2) <- diag(aux.nu2.2)-n.coords/nu2.2
    der2.Omega.star.inv.nu2.2 <- -Omega.star.inv%*%(
      aux.nu2.2)%*%
      Omega.star.inv

    trace2.nu2.2 <- -trace1.nu2.2-
      sum(B2.nu2.2*t(B2.nu2.2))
    H[ind.nu2.2,ind.nu2.2] <- -0.5*(q.f1/nu2.2-4*q.f2/(nu2.2^2)-4*q.f.nu2.2/(nu2.2^2))/sigma2+
      -0.5*(trace2.nu2.2+as.numeric(t(diff.y.tilde)%*%
                                      der2.Omega.star.inv.nu2.2%*%diff.y.tilde)/((nu2.2^2)*sigma2))


    return(H)
  }
  estim <- list()

  start_cov_pars[-(1:2)] <- start_cov_pars[-(1:2)]/start_cov_pars[1]
  start_par <- c(start_beta, log(start_cov_pars))
  if(method=="nlminb") {
    estNLMINB <- nlminb(start_par,
                        function(x) -log.lik(x),
                        function(x) -grad.log.lik(x),
                        function(x) -hessian.log.lik(x),
                        control=list(trace=1*messages))

    estim$estimate <- estNLMINB$par
    hess.MLE <- hessian.log.lik(estim$estimate)
    if(messages) {
      grad.MLE <- grad.log.lik(estim$estimate)
      cat("\n","Gradient at MLE:",
          paste(round(grad.MLE,10)),"\n")
    }
    estim$covariance <- solve(-hess.MLE)
    estim$log.lik <- -estNLMINB$objective

  } else if(method=="BFGS") {
    estimBFGS <- maxBFGS(log.lik,grad.log.lik,hessian.log.lik,
                         start_par,print.level=1*messages)
    estim$estimate <- estimBFGS$estimate
    if(messages) {
      cat("\n","Gradient at MLE:",
          paste(round(estimBFGS$gradient,10)),"\n")
    }
    estim$covariance <- solve(-estimBFGS$hessian)
    estim$log.lik <- estimBFGS$maximum
  }
  class(estim) <- "PrevMap"
  return(estim)
}
