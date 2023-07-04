
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


##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @importFrom maxLik maxBFGS
##' @importFrom Matrix Matrix forceSymmetric
glgm_lm <- function(y, D, coords, ID_coords, ID_re, s_unique, re_unique,
                    fix_tau2, method, start_beta, start_cov_pars) {

  m <- length(y)
  n <- nrow(coords)
  p <- ncol(D)
  U <- dist(coords)
  if(is.null(ID_re)) {
    n_re <- 0
  } else {
    n_re <- ncol(ID_re)
  }


  ID_g <- as.matrix(cbind(ID_coords, ID_re))

  n_dim_re <- sapply(1:(n_re+1), function(i) length(unique(ID_g[,i])))
  C_g <- matrix(0, nrow = m, ncol = sum(n_dim_re))

  for(i in 1:m) {
    ind_s_i <- which(s_unique==ID_g[i,1])
    C_g[i,1:n_dim_re[1]][ind_s_i] <- 1
  }

  if(n_re>0) {
    for(j in 1:n_re) {
      select_col <- sum(n_dim_re[1:j])

      for(i in 1:m) {
        ind_re_j_i <- which(re_unique[[j]]==ID_g[i,j+1])
        C_g[i,select_col+1:n_dim_re[j+1]][ind_re_j_i] <- 1
      }
    }
  }
  C_g <- Matrix(C_g, sparse = TRUE)

  C_g_m <- t(C_g)%*%C_g
  C_g_m <- forceSymmetric(C_g_m)

  log.lik <- function(par) {
    beta <- par[1:p]
    sigma2 <- exp(par[p+1])
    phi <- exp(par[p+2])
    if(!is.null(fix_tau2)) {
      nu2 <- fix_tau2/sigma2
      omega2 <- exp(par[p+3])
      if(n_re>0) {
        ind_sigma2_re <- (p+3+1):(p+3+n_re)
        sigma2_re <- exp(par[ind_sigma2_re])
      }
    } else {
      nu2 <- exp(par[p+3])
      omega2 <- exp(par[p+4])
      if(n_re>0) {
        ind_sigma2_re <- (p+4+1):(p+4+n_re)
        sigma2_re <- exp(par[ind_sigma2_re])
      }
    }

    R <- matern_cor(U,phi = phi, kappa=kappa,return_sym_matrix = TRUE)
    diag(R) <- diag(R)+nu2

    Sigma_g <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
    Sigma_g_inv <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
    Sigma_g[1:n_dim_re[1], 1:n_dim_re[1]] <- sigma2*R
    Sigma_g_inv[1:n_dim_re[1], 1:n_dim_re[1]] <-
    solve(R)/sigma2
    if(n_re > 0) {
      for(j in 1:n_re) {
        select_col <- sum(n_dim_re[1:j])

        diag(Sigma_g[select_col+1:n_dim_re[j+1], select_col+1:n_dim_re[j+1]]) <-
          sigma2_re[j]

        diag(Sigma_g_inv[select_col+1:n_dim_re[j+1], select_col+1:n_dim_re[j+1]]) <-
        1/sigma2_re[j]

      }
    }

    Sigma_g_inv <- Matrix(Sigma_g_inv, sparse = TRUE)
    Sigma_g_inv <- forceSymmetric(Sigma_g_inv)

    mu <- as.numeric(D%*%beta)
    diff.y <- y-mu
    diff.y.tilde <- as.numeric(t(C_g)%*%diff.y)
    Sigma_star <- Sigma_g_inv+C_g_m/omega2
    Sigma_star_inv <- forceSymmetric(solve(Sigma_star))

    q.f.y <- as.numeric(sum(diff.y^2)/omega2)
    q.f.y_tilde <- as.numeric(t(diff.y.tilde)%*%Sigma_star_inv%*%diff.y.tilde/
                                (omega2^2))
    Sigma_g_C_g_m <- Sigma_g%*%C_g_m
    Sigma_tilde <- Sigma_g_C_g_m/omega2
    diag(Sigma_tilde) <- diag(Sigma_tilde) + 1
    log_det <- as.numeric(m*log(omega2)+determinant(Sigma_tilde)$modulus)

    out <- -0.5*(log_det+q.f.y-q.f.y_tilde)
    return(out)
  }

  D.tilde <- t(D)%*%C_g
  U <- dist(coords)
  grad.log.lik <- function(par) {
    beta <- par[1:p]
    sigma2 <- exp(par[p+1])
    phi <- exp(par[p+2])
    if(!is.null(fix_tau2)) {
      nu2 <- fix_tau2/sigma2
      omega2 <- exp(par[p+3])
      if(n_re>0) {
        ind_sigma2_re <- (p+3+1):(p+3+n_re)
        ind_omega2 <- p+3
        sigma2_re <- exp(par[ind_sigma2_re])
      }
      g <- rep(0, p+3+n_re)
    } else {
      nu2 <- exp(par[p+3])
      omega2 <- exp(par[p+4])
      if(n_re>0) {
        ind_omega2 <- p+4
        ind_sigma2_re <- (p+4+1):(p+4+n_re)
        sigma2_re <- exp(par[ind_sigma2_re])
      }
      g <- rep(0, p+4+n_re)
    }

    R <- matern_cor(U,phi = phi, kappa=kappa,return_sym_matrix = TRUE)
    diag(R) <- diag(R)+nu2

    Sigma_g <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
    Sigma_g_inv <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
    Sigma_g[1:n_dim_re[1], 1:n_dim_re[1]] <- sigma2*R
    R.inv <- solve(R)
    Sigma_g_inv[1:n_dim_re[1], 1:n_dim_re[1]] <-
      R.inv/sigma2
    if(n_re > 0) {
      for(j in 1:n_re) {
        select_col <- sum(n_dim_re[1:j])

        diag(Sigma_g[select_col+1:n_dim_re[j+1], select_col+1:n_dim_re[j+1]]) <-
          sigma2_re[j]

        diag(Sigma_g_inv[select_col+1:n_dim_re[j+1], select_col+1:n_dim_re[j+1]]) <-
          1/sigma2_re[j]

      }
    }

    Sigma_g_inv <- Matrix(Sigma_g_inv, sparse = TRUE)
    Sigma_g_inv <- forceSymmetric(Sigma_g_inv)

    mu <- as.numeric(D%*%beta)
    diff.y <- y-mu
    diff.y.tilde <- as.numeric(t(C_g)%*%diff.y)
    Sigma_star <- Sigma_g_inv+C_g_m/omega2
    Sigma_star_inv <- forceSymmetric(solve(Sigma_star))
    M_aux <- D.tilde%*%Sigma_star_inv


    g[1:p] <- t(D)%*%diff.y/omega2-M_aux%*%diff.y.tilde/(omega2^2)

    der_Sigma_g_inv_sigma2 <- matrix(0, nrow = sum(n_dim_re),
                                     ncol = sum(n_dim_re))

    der_Sigma_g_inv_sigma2[1:n_dim_re[1], 1:n_dim_re[1]] <-
      -R.inv/sigma2^2
    der_sigma2_aux <- Sigma_star_inv%*%der_Sigma_g_inv_sigma2
    Sigma_g_C_g_m <- Sigma_g%*%C_g_m
    Sigma_tilde <- Sigma_g_C_g_m/omega2
    diag(Sigma_tilde) <- diag(Sigma_tilde) + 1
    Sigma_tilde_inv <- solve(Sigma_tilde)
    der_sigma2_Sigma_g <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
    der_sigma2_Sigma_g[1:n_dim_re[1], 1:n_dim_re[1]] <- R
    der_sigma2_Sigma_g <- der_sigma2_Sigma_g%*%C_g_m/omega2
    der_sigma2_trace <- sum(diag(Sigma_tilde_inv%*%der_sigma2_Sigma_g))
    g[p+1] <- (-0.5*der_sigma2_trace-0.5*t(diff.y.tilde)%*%
                      der_sigma2_aux%*%Sigma_star_inv%*%
              diff.y.tilde/(omega2^2))*sigma2

    der_Sigma_g_inv_phi <- matrix(0, nrow = sum(n_dim_re),
                                     ncol = sum(n_dim_re))
    M.der.phi <- matern.grad.phi(U, phi, kappa)
    der_Sigma_g_inv_phi[1:n_dim_re[1], 1:n_dim_re[1]] <-
      -M.der.phi*sigma2
    der_Sigma_g_inv_phi <- Sigma_g_inv%*%der_Sigma_g_inv_phi%*%Sigma_g_inv
    der_phi_aux <- Sigma_star_inv%*%der_Sigma_g_inv_phi
    der_phi_Sigma_g <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
    der_phi_Sigma_g[1:n_dim_re[1], 1:n_dim_re[1]] <- sigma2*M.der.phi
    der_phi_Sigma_g <- der_phi_Sigma_g%*%C_g_m/omega2
    der_phi_trace <- sum(diag(Sigma_tilde_inv%*%der_phi_Sigma_g))
    g[p+2] <- (-0.5*der_phi_trace-0.5*t(diff.y.tilde)%*%
                 der_phi_aux%*%Sigma_star_inv%*%
                 diff.y.tilde/(omega2^2))*phi
    if(is.null(fix_tau2)) {
      der_Sigma_g_inv_nu2 <- matrix(0, nrow = sum(n_dim_re),
                                    ncol = sum(n_dim_re))
      diag(der_Sigma_g_inv_nu2[1:n_dim_re[1], 1:n_dim_re[1]]) <-
        -sigma2
      der_Sigma_g_inv_nu2 <- Sigma_g_inv%*%der_Sigma_g_inv_nu2%*%Sigma_g_inv
      der_nu2_aux <- Sigma_star_inv%*%der_Sigma_g_inv_nu2
      der_nu2_Sigma_g <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
      diag(der_nu2_Sigma_g[1:n_dim_re[1], 1:n_dim_re[1]]) <- sigma2
      der_nu2_Sigma_g <- der_nu2_Sigma_g%*%C_g_m/omega2
      der_nu2_trace <- sum(diag(Sigma_tilde_inv%*%der_nu2_Sigma_g))
      g[p+3] <- (-0.5*der_nu2_trace-0.5*t(diff.y.tilde)%*%
                   der_nu2_aux%*%Sigma_star_inv%*%
                   diff.y.tilde/(omega2^2))*nu2
    }

    #q.f.y <- as.numeric(sum(diff.y^2)/omega2)
    #q.f.y_tilde <- as.numeric(t(diff.y.tilde)%*%Sigma_star_inv%*%diff.y.tilde/
    #                            (omega2^2))
    #Sigma_tilde <- Sigma_g_C_g_m/omega2
    #diag(Sigma_tilde) <- diag(Sigma_tilde) + 1
    #log_det <- as.numeric(m*log(omega2)+determinant(Sigma_tilde)$modulus)

    #out <- -0.5*(log_det+q.f.y-q.f.y_tilde)

    der_omega2_q.f.y <- -as.numeric(sum(diff.y^2)/omega2^2)
    der_omega2_Sigma_star <- -C_g_m/omega2^2
    der_omega2_Sigma_tilde <- -Sigma_g_C_g_m/omega2^2
    der_omega2_trace <- sum(diag(Sigma_tilde_inv%*%der_omega2_Sigma_tilde))
    der_omega2_q.f.y_tilde <- -2*as.numeric(t(diff.y.tilde)%*%Sigma_star_inv%*%diff.y.tilde/
                                (omega2^3))-
                              as.numeric(t(diff.y.tilde)%*%Sigma_star_inv%*%
                                           der_omega2_Sigma_star%*%
                                           Sigma_star_inv%*%diff.y.tilde/
                              (omega2^2))
    g[ind_omega2] <- (-0.5*(m/omega2+der_omega2_trace+
                           der_omega2_q.f.y-der_omega2_q.f.y_tilde))*omega2

    if(n_re>0) {
      for(i in 1:n_re) {
        select_col <- sum(n_dim_re[1:i])

        der_Sigma_g_inv_sigma2_re <- matrix(0, nrow = sum(n_dim_re),
                                            ncol = sum(n_dim_re))
        diag(der_Sigma_g_inv_sigma2_re[select_col+1:n_dim_re[i+1],
                                       select_col+1:n_dim_re[i+1]])  <-
          -1/sigma2_re[i]^2
        der_sigma2_re_aux <- Sigma_star_inv%*%der_Sigma_g_inv_sigma2_re
        der_sigma2_re_Sigma_g <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
        diag(der_sigma2_re_Sigma_g[select_col+1:n_dim_re[i+1],
                                   select_col+1:n_dim_re[i+1]])  <- 1
        der_sigma2_re_Sigma_g <- der_sigma2_re_Sigma_g%*%C_g_m/omega2
        der_sigma2_re_trace <- sum(diag(Sigma_tilde_inv%*%der_sigma2_re_Sigma_g))
        g[ind_sigma2_re[i]] <- (-0.5*der_sigma2_re_trace-0.5*t(diff.y.tilde)%*%
                                  der_sigma2_re_aux%*%Sigma_star_inv%*%
                                  diff.y.tilde/(omega2^2))*sigma2_re[i]
      }
    }
    return(g)
  }

  DtD <- t(D)%*%D
  Dty <- as.numeric(t(D)%*%y)
  D.tilde <- list()
  y.tilde <- list()
  for(i in 1:(1+n_re)) {
    D.tilde[[i]] <- sapply(1:p,function(j) as.numeric(tapply(D[,j],ID_g[,i],sum)))
    y.tilde[[i]] <- tapply(y,ID_coords,sum)
  }


  hessian.log.lik <- function(par) {
    beta <- par[1:p]
    sigma2 <- exp(par[p+1])
    phi <- exp(par[p+2])
    if(!is.null(fix_tau2)) {
      nu2.1 <- fix_tau2/sigma2
      omega2 <- exp(par[p+3])
      if(n_re>0) {
        nu2_re <- exp(par[(p+3+1):(p+3+n_re)])
      }
    } else {
      nu2.1 <- exp(par[p+3])
      omega2 <- exp(par[p+4])
      if(n_re>0) {
        nu2_re <- exp(par[(p+4+1):(p+4+n_re)])
      }
    }

    R <- matern_cor(U,phi = phi, kappa=kappa,
                    return_sym_matrix = TRUE)
    diag(R) <- diag(R)+nu2.1
    R.inv <- solve(R)
    Omega.star <- R.inv
    diag(Omega.star) <- (diag(Omega.star)+n.coords[[1]]/omega2)
    Omega.star.inv <- solve(Omega.star)

    mu <- as.numeric(D%*%beta)
    diff.y <- y-mu
    q.f1 <- sum(diff.y^2)

    diff.y.tilde <- list()
    q.f2 <- list()
    v <- list()
    v.beta.omega2 <- list()
    q.f <- q.f1/omega2
    grad.beta <- t(D)%*%(diff.y/omega2)/sigma2


    der.Omega.star.inv.omega2 <- -t(Omega.star.inv*n.coords[[1]]/omega2)%*%
      Omega.star.inv
    for(i in 1:(1+n_re)) {
      diff.y.tilde[[i]] <- as.numeric(tapply(diff.y, ID_g[,i], sum))
      if(i==1) {
        v[[i]] <- as.numeric(Omega.star.inv%*%diff.y.tilde[[i]])
        v.beta.omega2[[i]] <- as.numeric(der.Omega.star.inv.omega2%*%
                                          diff.y.tilde[[i]])
      } else {
        v[[i]] <- diff.y.tilde[[i]]/(1/nu2_re[i-1]+n.coords[[i]]/omega2)

        v.beta.omega2[[i]] <- (-n.coords[[i]]/omega2)*diff.y.tilde[[i]]

      }
      q.f2[[i]] <- as.numeric(t(diff.y.tilde[[i]])%*%v[[i]])
      q.f <- q.f-q.f2[[i]]/(omega2^2)

      grad.beta <- grad.beta - as.numeric(t(D)%*%v[[i]][ID_g[,i]]/
                                            (omega2^2))/sigma2
    }

    grad.sigma2 <- -0.5*(m-q.f/sigma2)

    M.det <- t(R*n.coords[[1]]/omega2)
    diag(M.det) <- diag(M.det)+1
    M.det.inv <- solve(M.det)

    R1.phi <- matern.grad.phi(U,phi,kappa)*phi
    Omega.star.R.inv <- Omega.star.inv%*%R.inv
    der.Omega.star.inv.phi <- Omega.star.R.inv%*%R1.phi%*%t(Omega.star.R.inv)
    der.M.phi <- t(R1.phi*n.coords[[1]]/omega2)
    M1.trace.phi <- M.det.inv*t(der.M.phi)
    trace1.phi <- sum(M1.trace.phi)
    v.beta.phi <- as.numeric(der.Omega.star.inv.phi%*%diff.y.tilde[[1]])
    q.f.phi <-  as.numeric(t(diff.y.tilde[[1]])%*%v.beta.phi)
    grad.phi <- -0.5*(trace1.phi-
                        q.f.phi/((omega2^2)*sigma2))

    der.M.omega2 <- -t(R*n.coords[[1]]/omega2)
    M1.trace.omega2 <- M.det.inv*t(der.M.omega2)
    trace1.omega2 <- sum(M1.trace.omega2)
    v.beta.omega2 <- as.numeric(der.Omega.star.inv.omega2%*%
                                 diff.y.tilde[[1]])
    q.f.omega2 <- as.numeric(t(diff.y.tilde[[1]])%*%v.beta.omega2)
    grad.omega2 <- -0.5*(-q.f1/omega2+2*q.f2[[1]]/(omega2^2))/sigma2+
      -0.5*(m+trace1.omega2+
              q.f.omega2/((omega2^2)*sigma2))

    if(n_re>0) {
      grad.nu2_re <- rep(NA, n_re)
      for(i in 1:n_re) {

        o.inv_re <- ((1/nu2_re[i]+n.coords[[i+1]]/omega2)*(omega2^2))
        der_nu2_re_o.inv_re <- -(omega2^2)/(nu2_re[i]^2)

        der_nu2_re_log.det_re_i <- sum((n.coords[[i+1]]/omega2)/
                                         (1+n.coords[[i+1]]*nu2_re[i]/omega2))

        grad.nu2_re[i] <-  (-0.5*(der_nu2_re_log.det_re_i+
                                    sum(der_nu2_re_o.inv_re*
                                          diff.y.tilde[[i+1]]^2/(sigma2*o.inv_re^2))))*nu2_re[i]

        der_omega2_re_o.inv_re <- (2*omega2/nu2_re[i]+n.coords[[i+1]])

        der_omega2_log.det_re_i <- -sum((n.coords[[i+1]]*nu2_re[i]/omega2^2)/
                                         (1+n.coords[[i+1]]*nu2_re[i]/omega2))

        grad.omega2 <- grad.omega2 + -0.5*(sum(der_omega2_re_o.inv_re*
                                               diff.y.tilde[[i+1]]^2/
                                               (sigma2*o.inv_re^2))*omega2+
                                           der_omega2_log.det_re_i*omega2)
      }
    }

    if(is.null(fix_tau2)) {
      der.Omega.star.inv.nu2.1 <- -Omega.star.R.inv%*%t(Omega.star.R.inv)*nu2.1
      trace1.nu2.1 <- sum(diag(M.det.inv)*n.coords[[1]]*nu2.1/omega2)
      v.beta.nu2.1 <- as.numeric(der.Omega.star.inv.nu2.1%*%
                                   diff.y.tilde[[1]])
      q.f.nu2.1 <- as.numeric(t(diff.y.tilde[[1]])%*%v.beta.nu2.1)
      grad.nu2.1 <- -0.5*(trace1.nu2.1+
                            q.f.nu2.1/((omega2^2)*sigma2))
    }


    ind.beta <- 1:p
    ind.sigma2 <- p+1
    ind.phi <- p+2
    if(is.null(fix_tau2)) {
      if(n_re>0) {
        g <- c(grad.beta,grad.sigma2,grad.phi,grad.nu2.1,grad.omega2,
               grad.nu2_re)
        H <- matrix(0,p+4+n_re,p+4+n_re)
        ind.nu2.1 <- p+3
        ind.omega2 <- p+4
        ind.nu2_re <- (p+4+1):(p+4+n_re)
      } else {
        g <- c(grad.beta,grad.sigma2,grad.phi,grad.nu2.1,grad.omega2)
        H <- matrix(0,p+4,p+4)
        ind.nu2.1 <- p+3
        ind.omega2 <- p+4
      }
    } else {
      if(n_re>0) {
        g <- c(grad.beta,grad.sigma2,grad.phi,grad.omega2,grad.nu2_re)
        H <- matrix(0,p+3+n_re,p+3+n_re)
        ind.omega2 <- p+3
        ind.nu2_re <- (p+3+1):(p+3+n_re)
      } else {
        g <- c(grad.beta,grad.sigma2,grad.phi,grad.omega2)
        H <- matrix(0,p+3,p+3)
        ind.omega2 <- p+3
      }
    }

    H[ind.beta,ind.beta] <- -DtD/(sigma2*omega2)
    for(i in 1:(n_re+1)) {
      if(i==1) {
        H[ind.beta,ind.beta] <- H[ind.beta,ind.beta] +
          (t(D.tilde[[i]])%*%Omega.star.inv%*%D.tilde[[i]]/(omega2^2))/sigma2
      } else {
        o.inv_re_aux <- ((1/nu2_re[i-1]+n.coords[[i]]/omega2))
        H[ind.beta,ind.beta] <- H[ind.beta,ind.beta] +
          (t(D.tilde[[i]])%*%((1/o.inv_re_aux)*D.tilde[[i]])/(omega2^2))/sigma2
      }
    }


    H[ind.beta,ind.sigma2] <- H[ind.sigma2,ind.beta] <- -grad.beta

    H[ind.beta,ind.phi] <- H[ind.phi,ind.beta] <- -t(D)%*%v.beta.phi[ID_g[,1]]/((omega2^2)*sigma2)

    if(is.null(fix_tau2)) {

      H[ind.beta,ind.nu2.1] <- H[ind.nu2.1,ind.beta] <- t(D)%*%v.beta.nu2.1[ID_g[,1]]/((omega2^2)*sigma2)

    }

    H[ind.beta,ind.omega2] <-
    H[ind.omega2,ind.beta] <- as.numeric(t(D)%*%(-diff.y/omega2+
                            2*v[[1]][ID_g[,1]]/(omega2^2))/sigma2)+
    t(D)%*%v.beta.omega2[ID_g[,1]]/((omega2^2)*sigma2)

    if(n_re>0) {
      for(i in 1:n_re) {
        o.inv_re <- ((1/nu2_re[i]+n.coords[[i+1]]/omega2)*(omega2^2))
        der_nu2_re_o.inv_re <- -(omega2^2)/(nu2_re[i]^2)


        der_omega2_re_o.inv_re <- (2*omega2/nu2_re[i]+n.coords[[i+1]])

        H[ind.beta,ind.omega2] <-
          H[ind.omega2,ind.beta] <- H[ind.omega2,ind.beta] +
          t(D.tilde[[i+1]]*(((der_omega2_re_o.inv_re/
                                (sigma2*o.inv_re^2))*omega2)))%*%
          (diff.y.tilde[[i+1]])

        H[ind.beta, ind.nu2_re[i]] <-
          H[ind.nu2_re[i], ind.beta] <-  t(D.tilde[[i+1]]*((der_nu2_re_o.inv_re/
                                                              (sigma2*o.inv_re^2))*nu2_re[i]))%*%
          (diff.y.tilde[[i+1]])

      }
    }

    H[ind.sigma2,ind.sigma2] <- -0.5*q.f/sigma2

    H[ind.sigma2,ind.phi] <- H[ind.phi,ind.sigma2] <- -(grad.phi+0.5*(trace1.phi))
    if(is.null(fix_tau2)) {
      H[ind.sigma2,ind.nu2.1] <- H[ind.nu2.1,ind.sigma2] <- -(grad.nu2.1+0.5*(trace1.nu2.1))
    }
    H[ind.sigma2,ind.omega2] <-
      0.5*(-q.f1/omega2+2*q.f2[[1]]/(omega2^2))/sigma2+
      0.5*(q.f.omega2/((omega2^2)*sigma2))

    if(n_re>0) {
      for(i in 1:n_re) {


        o.inv_re <- ((1/nu2_re[i]+n.coords[[i+1]]/omega2)*(omega2^2))
        der_nu2_re_o.inv_re <- -(omega2^2)/(nu2_re[i]^2)

        der_omega2_re_o.inv_re <- (2*omega2/nu2_re[i]+n.coords[[i+1]])


        H[ind.sigma2,ind.omega2] <- H[ind.sigma2,ind.omega2] +
          0.5*sum(der_omega2_re_o.inv_re*diff.y.tilde[[i+1]]^2/
         (sigma2*o.inv_re^2))*omega2

        H[ind.sigma2,ind.nu2_re[i]] <-
          H[ind.nu2_re[i],ind.sigma2] <- (0.5*(sum(der_nu2_re_o.inv_re*
                                                     diff.y.tilde[[i+1]]^2/
                                                     (sigma2*o.inv_re^2))))*
                                         nu2_re[i]

      }
    }
    H[ind.omega2,ind.sigma2] <- H[ind.sigma2,ind.omega2]

    R2.phi <- matern.hessian.phi(U,phi,kappa)*(phi^2)+
              R1.phi
    V1.phi <- -R.inv%*%R1.phi%*%R.inv
    V2.phi <- R.inv%*%(2*R1.phi%*%R.inv%*%R1.phi-R2.phi)%*%R.inv
    der2.Omega.star.inv.phi <- Omega.star.inv%*%(2*V1.phi%*%
                                                   Omega.star.inv%*%V1.phi-V2.phi)%*%
      Omega.star.inv
    M2.trace.phi <- M.det.inv*((R2.phi*n.coords[[1]]/omega2))
    B2.phi <- M.det.inv%*%der.M.phi

    trace2.phi <- sum(M2.trace.phi)-
      sum(B2.phi*t(B2.phi))
    H[ind.phi,ind.phi] <- -0.5*(trace2.phi-
                                  as.numeric(t(diff.y.tilde[[1]])%*%
                                               der2.Omega.star.inv.phi%*%
                                               diff.y.tilde[[1]])/((omega2^2)*sigma2))
    if(is.null(fix_tau2)) {

      V1.nu2.1 <- -R.inv%*%R.inv*nu2.1
      V2.phi.nu2.1 <- nu2.1*R.inv%*%(R1.phi%*%R.inv+R.inv%*%R1.phi)%*%R.inv
      der2.Omega.star.inv.phi.nu2.1 <- Omega.star.inv%*%(
        V1.phi%*%Omega.star.inv%*%V1.nu2.1+
          V1.nu2.1%*%Omega.star.inv%*%V1.phi-
          V2.phi.nu2.1)%*%
        Omega.star.inv
      B2.nu2.1 <- t(t(M.det.inv)*n.coords[[1]]*nu2.1/omega2)
      trace2.phi.nu2.1 <- -sum(B2.phi*t(B2.nu2.1))

      H[ind.phi,ind.nu2.1] <- H[ind.nu2.1,ind.phi] <- -0.5*(trace2.phi.nu2.1-
                                                              as.numeric(t(diff.y.tilde[[1]])%*%
                                                                           der2.Omega.star.inv.phi.nu2.1%*%
                                                                           diff.y.tilde[[1]])/((omega2^2)*sigma2))

    }

    der2.Omega.star.inv.phi.omega2 <- Omega.star.inv%*%(
      V1.phi%*%t(-Omega.star.inv*
                   n.coords[[1]]/omega2)+
        (-Omega.star.inv*
           n.coords[[1]]/omega2)%*%V1.phi)%*%
      Omega.star.inv
    B2.omega2 <- M.det.inv%*%der.M.omega2
    trace2.phi.omega2 <- -trace1.phi-sum(B2.phi*t(B2.omega2))
    H[ind.phi,ind.omega2] <- H[ind.omega2,ind.phi] <- -0.5*(trace2.phi.omega2+
                                                            -as.numeric(t(diff.y.tilde[[1]])%*%
                                                             der2.Omega.star.inv.phi.omega2%*%
                                                               diff.y.tilde[[1]])/((omega2^2)*sigma2)+
                                                            +2*q.f.phi/((omega2^2)*sigma2))


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
                                        as.numeric(t(diff.y.tilde[[1]])%*%
                                                     der2.Omega.star.inv.nu2.1%*%
                                                     diff.y.tilde[[1]])/((omega2^2)*sigma2))

      der2.Omega.star.inv.nu2.1.omega2 <- Omega.star.inv%*%(
        V1.nu2.1%*%t(-Omega.star.inv*
                       n.coords[[1]]/omega2)+
          (-Omega.star.inv*
             n.coords[[1]]/omega2)%*%V1.nu2.1)%*%
        Omega.star.inv

      trace2.nu2.1.omega2 <- -trace1.nu2.1-sum(B2.nu2.1*t(B2.omega2))
      H[ind.nu2.1,ind.omega2] <- H[ind.omega2,ind.nu2.1] <- -0.5*(trace2.nu2.1.omega2+
                                                                  -as.numeric(t(diff.y.tilde[[1]])%*%
                                                          der2.Omega.star.inv.nu2.1.omega2%*%
                                                            diff.y.tilde[[1]])/((omega2^2)*sigma2)+
                                                                  -2*q.f.nu2.1/((omega2^2)*sigma2))

    }

    aux.omega2 <- 2*t(Omega.star.inv*n.coords[[1]])*n.coords[[1]]/
      (omega2^2)
    diag(aux.omega2) <- diag(aux.omega2)-n.coords[[1]]/omega2
    der2.Omega.star.inv.omega2 <- -Omega.star.inv%*%(
      aux.omega2)%*%
      Omega.star.inv

    trace2.omega2 <- -trace1.omega2-
      sum(B2.omega2*t(B2.omega2))
    H[ind.omega2,ind.omega2] <- -0.5*(q.f1/omega2-4*q.f2[[1]]/(omega2^2)-4*q.f.omega2/(omega2^2))/sigma2+
      -0.5*(trace2.omega2+as.numeric(t(diff.y.tilde[[1]])%*%
                                      der2.Omega.star.inv.omega2%*%
                                      diff.y.tilde[[1]])/
              ((omega2^2)*sigma2))
    if(n_re>0) {
      for(i in 1:n_re) {

        o.inv_re <- ((1/nu2_re[i]+n.coords[[i+1]]/omega2)*(omega2^2))

        der_omega2_re_o.inv_re <- (2*omega2/nu2_re[i]+n.coords[[i+1]])
        der2_omega2_re_o.inv_re <- 2/nu2_re[i]
        der_nu2_re_o.inv_re <- -(omega2^2)/(nu2_re[i]^2)
        der2_nu2_re_o.inv_re <- 2*(omega2^2)/(nu2_re[i]^3)
        der_nu2_re_omega2_o.inv_re <- -2*(omega2)/(nu2_re[i]^2)

        den_aux <- 1+n.coords[[i+1]]*nu2_re[i]/omega2
        a_aux <- -n.coords[[i+1]]*nu2_re[i]/omega2^2
        b_aux <- 2*n.coords[[i+1]]*nu2_re[i]/omega2^3
        der_omega2_log.det_re_i <- sum(a_aux/den_aux)
        der2_omega2_log.det_re_i <- sum((b_aux*den_aux-a_aux^2)/
                                         (den_aux^2))
        der_nu2_re_omega2_log.det_re_i <- sum(((-n.coords[[i+1]]/(omega2^2))*den_aux+
                                                -(n.coords[[i+1]]/omega2)*a_aux)/
                                               (den_aux^2))
        der_nu2_re_log.det_re_i <- sum((n.coords[[i+1]]/omega2)/
                                         (1+n.coords[[i+1]]*nu2_re[i]/omega2))

        der2_nu2_re_log.det_re_i <- -sum((n.coords[[i+1]]/omega2)^2/
                                         (1+n.coords[[i+1]]*nu2_re[i]/omega2)^2)

        H[ind.omega2,ind.omega2] <-
          H[ind.omega2,ind.omega2] + -0.5*(sum(der_omega2_re_o.inv_re*
                                               diff.y.tilde[[i+1]]^2/
                                               (sigma2*o.inv_re^2))*omega2+
                                         sum((der2_omega2_re_o.inv_re*o.inv_re+
                                            -2*der_omega2_re_o.inv_re^2)*
                                               diff.y.tilde[[i+1]]^2/
                                                 (sigma2*o.inv_re^3))*(omega2^2)+
                                           der_omega2_log.det_re_i*omega2+
                                           der2_omega2_log.det_re_i*omega2^2)

        #der_omega2_log.det_re_i <- -sum((n.coords[[i+1]]*nu2_re[i]/omega2^2)/
        #                                 (1+n.coords[[i+1]]*nu2_re[i]/omega2))


        #der_nu2_re_o.inv_re <- -(omega2^2)/(nu2_re[i]^2)

        #grad.nu2_re[i] <-  (-0.5*(der_nu2_re_log.det_re_i+
        #                            sum(der_nu2_re_o.inv_re*
        #                                  diff.y.tilde[[i+1]]^2/(sigma2*o.inv_re^2))))*nu2_re[i]

        #der_nu2_re_log.det_re_i <- sum((n.coords[[i+1]]/omega2)/
        #                                 (1+n.coords[[i+1]]*nu2_re[i]/omega2))


        H[ind.omega2,ind.nu2_re[i]] <-
        H[ind.nu2_re[i],ind.omega2] <-
          (-0.5*(der_nu2_re_omega2_log.det_re_i+
          sum((der_nu2_re_omega2_o.inv_re*o.inv_re+
          -2*der_nu2_re_o.inv_re*der_omega2_re_o.inv_re)*
          diff.y.tilde[[i+1]]^2/(sigma2*o.inv_re^3))))*
          omega2*nu2_re[i]

        #der_nu2_re_o.inv_re <- -(omega2^2)/(nu2_re[i]^2)
        #der2_nu2_re_o.inv_re <- 2*(omega2^2)/(nu2_re[i]^3)

        #der_nu2_re_log.det_re_i <- sum((n.coords[[i+1]]/omega2)/
        #                                 (1+n.coords[[i+1]]*nu2_re[i]/omega2))
        #der2_nu2_re_log.det_re_i <- -sum((n.coords[[i+1]]/omega2)^2/
        #                                 (1+n.coords[[i+1]]*nu2_re[i]/omega2)^2)

        H[ind.nu2_re[i], ind.nu2_re[i]] <-
          (-0.5*(der_nu2_re_log.det_re_i+
          sum(der_nu2_re_o.inv_re*
          diff.y.tilde[[i+1]]^2/(sigma2*o.inv_re^2))))*nu2_re[i]+
          (-0.5*(der2_nu2_re_log.det_re_i+
                sum((der2_nu2_re_o.inv_re*o.inv_re-
                     2*der_nu2_re_o.inv_re^2)*
                diff.y.tilde[[i+1]]^2/(sigma2*o.inv_re^3))))*(nu2_re[i]^2)

      }
    }
    return(H)
  }

  start_cov_pars[-(1:2)] <- start_cov_pars[-(1:2)]/start_cov_pars[1]
  start_par <- c(start_beta, log(start_cov_pars))

  estNLMINB <- nlminb(start_par,
                        function(x) -log.lik(x),
                        control=list(trace=1*messages))

  estim$estimate <- estNLMINB$par
  hess.MLE <- hessian.log.lik(estim$estimate)

  grad.MLE <- grad.log.lik(estim$estimate)
      cat("\n","Gradient at MLE:",
          paste(round(grad.MLE,10)),"\n")

  estim$covariance <- solve(-hess.MLE)
  estim$log.lik <- -estNLMINB$objective
  estim$covariance <- solve(-estimBFGS$hessian)
  estim$log.lik <- estimBFGS$maximum

  class(estim) <- "RiskMap"
  return(estim)
}
