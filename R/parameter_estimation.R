##' @title Estimation of Generalized Linear Gaussian Process Models
##' @description Fits generalized linear Gaussian process models to spatial data, incorporating spatial Gaussian processes with a Matern correlation function. Supports Gaussian, binomial, and Poisson response families.
##' @param formula A formula object specifying the model to be fitted. The formula should include fixed effects, random effects (specified using \code{re()}), and spatial effects (specified using \code{gp()}).
##' @param data A data frame or sf object containing the variables in the model.
##' @param family A character string specifying the distribution of the response variable. Must be one of "gaussian", "binomial", or "poisson".
##' @param den Optional offset for binomial or Poisson distributions. If not provided, defaults to 1 for binomial.
##' @param cov_offset Optional numeric vector for covariate offset.
##' @param crs Optional integer specifying the Coordinate Reference System (CRS) if data is not an sf object. Defaults to 4326 (long/lat).
##' @param convert_to_crs Optional integer specifying a CRS to convert the spatial coordinates.
##' @param scale_to_km Logical indicating whether to scale coordinates to kilometers. Defaults to TRUE.
##' @param control_mcmc Control parameters for MCMC sampling. Must be an object of class "mcmc.RiskMap" as returned by \code{\link{set_control_sim}}.
##' @param par0 Optional list of initial parameter values for the MCMC algorithm.
##' @param S_samples Optional matrix of pre-specified sample paths for the spatial random effect.
##' @param return_samples Logical indicating whether to return MCMC samples when fitting a Binomial or Poisson model. Defaults to FALSE.
##' @param messages Logical indicating whether to print progress messages. Defaults to TRUE.
##' @param fix_var_me Optional fixed value for the measurement error variance.
##' @param start_pars Optional list of starting values for model parameters: beta (regression coefficients), sigma2 (spatial process variance), tau2 (nugget effect variance), phi (spatial correlation scale), sigma2_me (measurement error variance), and sigma2_re (random effects variances).
##' @details
##' Generalized linear Gaussian process models extend generalized linear models (GLMs) by incorporating spatial Gaussian processes to account for spatial correlation in the data. This function fits GLGPMs using maximum likelihood methods, allowing for Gaussian, binomial, and Poisson response families.
##' In the case of the Binomial and Poisson families, a Monte Carlo maximum likelihood algorithm is used.
##'
##' The spatial Gaussian process is modeled with a Matern correlation function, which is flexible and commonly used in geostatistical modeling. The function supports both spatial covariates and unstructured random effects, providing a comprehensive framework to analyze spatially correlated data across different response distributions.
##'
##' Additionally, the function allows for the inclusion of unstructured random effects, specified through the \code{re()} term in the model formula. These random effects can capture unexplained variability at specific locations beyond the fixed and spatial covariate effects, enhancing the model's flexibility in capturing complex spatial patterns.
##'
##' The \code{convert_to_crs} argument can be used to reproject the spatial coordinates to a different CRS. The \code{scale_to_km} argument scales the coordinates to kilometers if set to TRUE.
##'
##' The \code{control_mcmc} argument specifies the control parameters for MCMC sampling. This argument must be an object returned by \code{\link{set_control_sim}}.
##'
##' The \code{start_pars} argument allows for specifying starting values for the model parameters. If not provided, default starting values are used.
##'
##' @return An object of class "RiskMap" containing the fitted model and relevant information:
##' \item{y}{Response variable.}
##' \item{D}{Covariate matrix.}
##' \item{coords}{Unique spatial coordinates.}
##' \item{ID_coords}{Index of coordinates.}
##' \item{re}{Random effects.}
##' \item{ID_re}{Index of random effects.}
##' \item{fix_tau2}{Fixed nugget effect variance.}
##' \item{fix_var_me}{Fixed measurement error variance.}
##' \item{formula}{Model formula.}
##' \item{family}{Response family.}
##' \item{crs}{Coordinate Reference System.}
##' \item{scale_to_km}{Indicator if coordinates are scaled to kilometers.}
##' \item{data_sf}{Original data as an sf object.}
##' \item{kappa}{Spatial correlation parameter.}
##' \item{units_m}{Distribution offset for binomial/Poisson.}
##' \item{cov_offset}{Covariate offset.}
##' \item{call}{Matched call.}
##' @seealso \code{\link{set_control_sim}}, \code{\link{summary.RiskMap}}, \code{\link{to_table}}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##' @importFrom sf st_crs st_as_sf st_drop_geometry
##' @export
glgpm <- function(formula,
                 data,
                 family,
                 den = NULL,
                 cov_offset = NULL,
                 crs = NULL, convert_to_crs = NULL,
                 scale_to_km = TRUE,
                 control_mcmc = set_control_sim(),
                 par0=NULL,
                 S_samples = NULL,
                 return_samples = TRUE,
                 messages = TRUE,
                 fix_var_me = NULL,
                 start_pars = list(beta = NULL,
                                   sigma2 = NULL,
                                   tau2 = NULL,
                                   phi = NULL,
                                   sigma2_me = NULL,
                                   sigma2_re = NULL)) {

  nong <- family=="binomial" | family=="poisson"


  if(!inherits(formula,
               what = "formula", which = FALSE)) {
    stop("'formula' must be a 'formula'
         object indicating the variables of the
         model to be fitted")
  }

  inter_f <- interpret.formula(formula)

  if(length(crs)>0) {
    if(!is.numeric(crs) |
       (is.numeric(crs) &
        (crs%%1!=0 | crs <0))) stop("'crs' must be a positive integer number")
  }
  if(class(data)[1]=="data.frame") {
    if(is.null(crs)) {
      warning("'crs' is set to 4326 (long/lat)")
      crs <- 4326
    }
    if(length(inter_f$gp.spec$term)==2) {
      new_x <- paste(inter_f$gp.spec$term[1],"_sf",sep="")
      new_y <- paste(inter_f$gp.spec$term[2],"_sf",sep="")
      data[[new_x]] <-  data[[inter_f$gp.spec$term[1]]]
      data[[new_y]] <-  data[[inter_f$gp.spec$term[2]]]
      data <- st_as_sf(data,
                       coords = c(new_x, new_y),
                       crs = crs)
    }
  }

  if(length(inter_f$gp.spec$term) == 1 & inter_f$gp.spec$term[1]=="sf" &
     class(data)[1]!="sf") stop("'data' must be an object of class 'sf'")


  if(class(data)[1]=="sf") {
    if(is.na(st_crs(data)) & is.null(crs)) {
      stop("the CRS of the sf object passed to 'data' is missing and and is not specified through 'crs'")
    } else if(is.na(st_crs(data))) {
      data <- st_as_sf(data, crs = crs)
    }
  }


  kappa <- inter_f$gp.spec$kappa
  if(kappa < 0) stop("kappa must be positive.")

  if(family != "gaussian" & family != "binomial" &
     family != "poisson") stop("'family' must be either 'gaussian', 'binomial'
                               or 'poisson'")


  mf <- model.frame(inter_f$pf,data=data, na.action = na.fail)

  # Extract outcome data
  y <- as.numeric(model.response(mf))
  n <- length(y)

  # Extract covariates matrix
  D <- as.matrix(model.matrix(attr(mf,"terms"),data=data))
  if(is.null(cov_offset)) {
    cov_offset <- 0
  } else {
    if(!is.numeric(cov_offset)) stop("the variable passed to 'cov_offset'
                                     must be numeric vector")
    if(any(is.na(cov_offset))) stop("missing values not accepted in the offset")
    if(length(cov_offset)!=n) stop("the offset values do not match the number of observations in the data")
  }

  # Define distributional offset for Binomial and Poisson distributions
  if(nong) {
    do_name <- deparse(substitute(den))
    if(do_name=="NULL") {
      units_m <- 1
      if(family=="binomial") warning("'den' is assumed to be 1 for all observations \n")
    } else {
      units_m <- data[[do_name]]
    }
    if(is.integer(units_m)) units_m <- as.numeric(units_m)
    if(!is.numeric(units_m)) stop("the variable passed to `den` must be numeric")
    if(family=="binomial" & any(y > units_m)) stop("The counts identified by the outcome variable cannot be larger
                              than `den` in the case of a Binomial distribution")
    if(!inherits(control_mcmc,
                 what = "mcmc.RiskMap", which = FALSE)) {
      stop ("the argument passed to 'control_mcmc' must be an output
                                                  from the function set_control_sim; see ?set_control_sim
                                                  for more details")

    }

  }

  if(length(inter_f$re.spec) > 0) {
    hr_re <- inter_f$re.spec$term
    re_names <- inter_f$re.spec$term
  } else {
    hr_re <- NULL
  }


  if(!is.null(hr_re)) {
    # Define indices of random effects
    re_mf <- st_drop_geometry(data[hr_re])
    re_mf_n <- re_mf

    if(any(is.na(re_mf))) stop("Missing values in the variable(s) of the random effects specified through re() ")
    names_re <- colnames(re_mf)
    n_re <- ncol(re_mf)

    ID_re <- matrix(NA, nrow = n, ncol = n_re)
    re_unique <- list()
    re_unique_f <- list()
    for(i in 1:n_re) {
      if(is.factor(re_mf[,i])) {
        re_mf_n[,i] <- as.numeric(re_mf[,i])
        re_unique[[names_re[i]]] <- 1:length(levels(re_mf[,i]))
        ID_re[, i] <- sapply(1:n,
                             function(j) which(re_mf_n[j,i]==re_unique[[names_re[i]]]))
        re_unique_f[[names_re[i]]] <-levels(re_mf[,i])
      } else if(is.numeric(re_mf[,i])) {
        re_unique[[names_re[i]]] <- unique(re_mf[,i])
        ID_re[, i] <- sapply(1:n,
                             function(j) which(re_mf_n[j,i]==re_unique[[names_re[i]]]))
        re_unique_f[[names_re[i]]] <- re_unique[[names_re[i]]]
      }
    }
    ID_re <- data.frame(ID_re)
    colnames(ID_re) <- re_names
  } else {
    n_re <- 0
    re_unique <- NULL
    ID_re <- NULL
  }


  # Extract coordinates
  if(!is.null(convert_to_crs)) {
    if(!is.numeric(convert_to_crs)) stop("'convert_to_utm' must be a numeric object")
    data <- st_transform(data, crs = convert_to_crs)
    crs <- convert_to_crs
  }
  if(messages) message("The CRS used is", as.list(st_crs(data))$input, "\n")

  coords_o <- st_coordinates(data)
  coords <- unique(coords_o)

  m <- nrow(coords_o)
  ID_coords <- sapply(1:m, function(i)
               which(coords_o[i,1]==coords[,1] &
                     coords_o[i,2]==coords[,2]))
  s_unique <- unique(ID_coords)

  fix_tau2 <- inter_f$gp.spec$nugget

  if(all(table(ID_coords)==1) &
    is.null(family=="gaussian" && is.null(fix_tau2)) & is.null(fix_var_me)) {
    stop("When there is only one observation per location, both the nugget and measurement error cannot
         be estimate. Consider removing either one of them. ")
  }

  if(scale_to_km) {
    coords_o <- coords_o/1000
    coords <- coords/1000
    if(messages) message("Distances between locations are computed in kilometers \n")
  } else {
    if(messages) message("Distances between locations are computed in meters \n")
  }


  if(is.null(start_pars$beta)) {
    if(family=="gaussian") {
      start_pars$beta <- as.numeric(solve(t(D)%*%D)%*%t(D)%*%y)
    } else if(family=="binomial") {
      aux_data <- data.frame(y=y, units_m = units_m, D[,-1])
      if(length(cov_offset)==1) cov_offset_aux <- rep(cov_offset, n)
      glm_fitted <- glm(cbind(y, units_m - y) ~ ., offset = cov_offset_aux,
                        data = aux_data, family = binomial)
      start_pars$beta <- stats::coef(glm_fitted)
    } else if(family=="poisson") {
      pf_aux <- stats::update(inter_f$pf, . ~ . + offset(log(units_m)) + offset(cov_offset))
      data_aux <- data
      data_aux$units_m <- units_m; data_aux$cov_offset <- cov_offset
      glm_fitted <- glm(pf_aux, data = data_aux, family = poisson)
      start_pars$beta <- stats::coef(glm_fitted)
    }
  } else {
    if(length(start_pars$beta)!=ncol(D)) stop("number of starting values provided
                                              for 'beta' do not match the number of
                                              covariates specified in the model,
                                              including the intercept")
  }

  if(is.null(start_pars$sigma2)) {
    start_pars$sigma2 <- 1
  } else {
    if(start_pars$sigma2<0) stop("the starting value for sigma2 must be positive")
  }

  if(is.null(start_pars$phi)) {
    start_pars$phi <- quantile(dist(coords),0.1)
  } else {
    if(start_pars$phi<0) stop("the starting value for phi must be positive")
  }

  if(is.null(fix_tau2)) {
    if(is.null(start_pars$tau2)) {
      start_pars$tau2 <- 1
    } else {
      if(start_pars$tau2<0) stop("the starting value for tau2 must be positive")
    }
  }

  if(n_re > 0) {
    if(is.null(start_pars$sigma2_re)) {
      start_pars$sigma2_re <- rep(1,n_re)
    } else {
      if(length(start_pars$sigma2_re)!=n_re) stop("starting values for 'sigma2_re' do not
                                       match the number of specified unstructured
                                       random effects")
      if(any(start_pars$sigma2_re<0)) stop("all the starting values for sigma2_re must be positive")
    }
  }


  if(!is.null(start_pars$beta)) {
    if(length(start_pars$beta)!=ncol(D)) stop("The values passed to 'start_beta' do not match
                                  the covariates passed to the 'formula'.")
  } else {
    start_pars$beta <- as.numeric(solve(t(D)%*%D)%*%t(D)%*%y)
  }



  if(!nong) {
    if(is.null(fix_var_me)) {
      if(is.null(start_pars$sigma2_me)) {
        start_pars$sigma2_me <- 1
      } else {
        if(start_pars$sigma2_me<0) stop("the starting value for sigma2_me must be positive")
      }
    }
    res <- glgpm_lm(y = y-cov_offset, D, coords, kappa = inter_f$gp.spec$kappa,
            ID_coords, ID_re, s_unique, re_unique,
            fix_var_me, fix_tau2,
            start_beta = start_pars$beta,
            start_cov_pars = c(start_pars$sigma2,
                               start_pars$phi,
                               start_pars$tau2,
                               start_pars$sigma2_re,
                               start_pars$sigma2_me),
            messages = messages)
  } else if(nong) {
    if(is.null(par0)) {
      par0 <- start_pars
    } else {
      if(length(par0$beta)!=ncol(D)) stop("the values passed to `beta` in par0 do not match the
                                          variables specified in the formula")
    }
    res <- glgpm_nong(y = y, D, coords, units_m, kappa = inter_f$gp.spec$kappa,
                        ID_coords, ID_re, s_unique, re_unique,
                        fix_tau2, family = family,
                        return_samples = return_samples,
                        par0 = par0, cov_offset = cov_offset,
                        start_beta = start_pars$beta,
                        start_cov_pars = c(start_pars$sigma2,
                                           start_pars$phi,
                                           start_pars$tau2,
                                           start_pars$sigma2_re),
                        control_mcmc = control_mcmc,
                        messages = messages)
  }

  res$y <- y
  res$D <- D
  res$coords <- coords
  res$ID_coords <- ID_coords
  if(n_re>0) {
    res$re <- re_unique_f
    res$ID_re <- as.data.frame(ID_re)
    colnames(res$ID_re) <- names_re
  }
  res$fix_tau2 <- fix_tau2
  res$fix_var_me <- fix_var_me
  res$formula <- formula
  res$family <- family
  if(!is.null(convert_to_crs)) {
    crs <- convert_to_crs
  } else {
    crs <- sf::st_crs(data)$input
  }
  res$crs <- crs
  res$scale_to_km <- scale_to_km
  res$data_sf <- data
  res$kappa <- kappa
  if(nong) res$units_m <- units_m
  res$cov_offset <- cov_offset
  res$call <- match.call()
  return(res)
}


##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @importFrom Matrix Matrix forceSymmetric
glgpm_lm <- function(y, D, coords, kappa, ID_coords, ID_re, s_unique, re_unique,
                     fix_var_me, fix_tau2, start_beta, start_cov_pars, messages) {

  m <- length(y)
  p <- ncol(D)
  U <- dist(coords)
  if(is.null(ID_re)) {
    n_re <- 0
  } else {
    n_re <- ncol(ID_re)
  }

  if(!is.null(fix_var_me)) {
    if(fix_var_me==0) {
      fix_var_me <- 10e-10
    }
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
  C_g <- Matrix(C_g, sparse = TRUE, doDiag = FALSE)


  C_g_m <- Matrix::t(C_g)%*%C_g
  C_g_m <- forceSymmetric(C_g_m)

  ind_beta <- 1:p

  ind_sigma2 <- p+1

  ind_phi <- p+2

  if(!is.null(fix_tau2)) {
    ind_omega2 <- p+3
    if(n_re>0) {
      ind_sigma2_re <- (p+3+1):(p+3+n_re)
    }
  } else {
    ind_nu2 <- p+3
    ind_omega2 <- p+4
    if(n_re>0) {
      ind_omega2 <- p+4
      ind_sigma2_re <- (p+4+1):(p+4+n_re)
    }
  }


  log.lik <- function(par) {
    beta <- par[ind_beta]
    sigma2 <- exp(par[p+1])
    phi <- exp(par[ind_phi])
    if(!is.null(fix_tau2)) {
      nu2 <- fix_tau2/sigma2
    } else {
      nu2 <- exp(par[ind_nu2])
    }
    if(n_re>0) {
      sigma2_re <- exp(par[ind_sigma2_re])
    }

    if(!is.null(fix_var_me)) {
      omega2 <- fix_var_me
    } else {
      omega2 <- exp(par[ind_omega2])
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
    diff.y.tilde <- as.numeric(Matrix::t(C_g)%*%diff.y)
    Sigma_star <- Sigma_g_inv+C_g_m/omega2
    Sigma_star_inv <- forceSymmetric(solve(Sigma_star))

    q.f.y <- as.numeric(sum(diff.y^2)/omega2)
    q.f.y_tilde <- as.numeric(t(diff.y.tilde)%*%Sigma_star_inv%*%diff.y.tilde/
                                (omega2^2))
    Sigma_g_C_g_m <- Sigma_g%*%C_g_m
    Sigma_tilde <- Sigma_g_C_g_m/omega2
    Matrix::diag(Sigma_tilde) <- Matrix::diag(Sigma_tilde) + 1
    log_det <- as.numeric(m*log(omega2)+Matrix::determinant(Sigma_tilde)$modulus)

    out <- -0.5*(log_det+q.f.y-q.f.y_tilde)
    return(out)
  }

  D.tilde <- t(D)%*%C_g
  U <- dist(coords)

  grad.log.lik <- function(par) {
    beta <- par[ind_beta]
    sigma2 <- exp(par[ind_sigma2])
    phi <- exp(par[ind_phi])
    if(!is.null(fix_tau2)) {
      nu2 <- fix_tau2/sigma2
    } else {
      nu2 <- exp(par[ind_nu2])
    }
    if(n_re>0) {
      sigma2_re <- exp(par[ind_sigma2_re])
    }
    if(!is.null(fix_var_me)) {
      omega2 <- fix_var_me
    } else {
      omega2 <- exp(par[ind_omega2])
    }

    n_p <- length(par)
    g <- rep(0, n_p)

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
    diff.y.tilde <- as.numeric(Matrix::t(C_g)%*%diff.y)
    Sigma_star <- Sigma_g_inv+C_g_m/omega2
    Sigma_star_inv <- forceSymmetric(solve(Sigma_star))
    M_aux <- D.tilde%*%Sigma_star_inv


    g[ind_beta] <- t(D)%*%diff.y/omega2-M_aux%*%diff.y.tilde/(omega2^2)

    der_Sigma_g_inv_sigma2 <- matrix(0, nrow = sum(n_dim_re),
                                     ncol = sum(n_dim_re))

    der_Sigma_g_inv_sigma2[1:n_dim_re[1], 1:n_dim_re[1]] <-
      -R.inv/sigma2^2
    der_sigma2_aux <- Sigma_star_inv%*%der_Sigma_g_inv_sigma2
    Sigma_g_C_g_m <- Sigma_g%*%C_g_m
    Sigma_tilde <- Sigma_g_C_g_m/omega2
    Matrix::diag(Sigma_tilde) <- Matrix::diag(Sigma_tilde) + 1
    Sigma_tilde_inv <- solve(Sigma_tilde)
    der_sigma2_Sigma_g <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
    der_sigma2_Sigma_g[1:n_dim_re[1], 1:n_dim_re[1]] <- R
    der_sigma2_Sigma_g <- der_sigma2_Sigma_g%*%C_g_m/omega2
    der_sigma2_trace <- sum(Matrix::diag(Sigma_tilde_inv%*%der_sigma2_Sigma_g))
    g[ind_sigma2] <- (-0.5*der_sigma2_trace-0.5*t(diff.y.tilde)%*%
                        der_sigma2_aux%*%Sigma_star_inv%*%
                        diff.y.tilde/(omega2^2))*sigma2

    der_R_phi <- matrix(0, nrow = sum(n_dim_re),
                        ncol = sum(n_dim_re))
    M.der.phi <- matern.grad.phi(U, phi, kappa)
    der_R_phi[1:n_dim_re[1], 1:n_dim_re[1]] <-
      M.der.phi*sigma2
    der_Sigma_g_inv_phi <- Sigma_g_inv%*%der_R_phi%*%Sigma_g_inv
    der_phi_aux <- -Sigma_star_inv%*%der_Sigma_g_inv_phi
    der_phi_Sigma_g <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
    der_phi_Sigma_g[1:n_dim_re[1], 1:n_dim_re[1]] <- sigma2*M.der.phi
    der_phi_Sigma_g <- der_phi_Sigma_g%*%C_g_m/omega2
    der_phi_trace <- sum(Matrix::diag(Sigma_tilde_inv%*%der_phi_Sigma_g))
    g[p+2] <- (-0.5*der_phi_trace-0.5*t(diff.y.tilde)%*%
                 der_phi_aux%*%Sigma_star_inv%*%
                 diff.y.tilde/(omega2^2))*phi
    if(is.null(fix_tau2)) {
      der_R_nu2 <- matrix(0, nrow = sum(n_dim_re),
                          ncol = sum(n_dim_re))
      diag(der_R_nu2[1:n_dim_re[1], 1:n_dim_re[1]]) <-
        sigma2
      der_Sigma_g_inv_nu2_aux <- Sigma_g_inv%*%der_R_nu2%*%Sigma_g_inv
      der_nu2_aux <- -Sigma_star_inv%*%der_Sigma_g_inv_nu2_aux
      der_nu2_Sigma_g <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
      diag(der_nu2_Sigma_g[1:n_dim_re[1], 1:n_dim_re[1]]) <- sigma2
      der_nu2_Sigma_g <- der_nu2_Sigma_g%*%C_g_m/omega2
      der_nu2_trace <- sum(Matrix::diag(Sigma_tilde_inv%*%der_nu2_Sigma_g))
      g[p+3] <- (-0.5*der_nu2_trace-0.5*t(diff.y.tilde)%*%
                   der_nu2_aux%*%Sigma_star_inv%*%
                   diff.y.tilde/(omega2^2))*nu2
    }


    if(is.null(fix_var_me)) {
      der_omega2_q.f.y <- -as.numeric(sum(diff.y^2)/omega2^2)
      der_omega2_Sigma_star <- -C_g_m/omega2^2
      der_omega2_Sigma_tilde <- -Sigma_g_C_g_m/omega2^2
      der_omega2_trace <- sum(Matrix::diag(Sigma_tilde_inv%*%der_omega2_Sigma_tilde))
      M_beta_omega2 <- -Sigma_star_inv%*%der_omega2_Sigma_star%*%Sigma_star_inv

      num1_omega2 <- as.numeric(t(diff.y.tilde)%*%Sigma_star_inv%*%diff.y.tilde)
      der_num1_omega2 <- as.numeric(t(diff.y.tilde)%*%
                                      M_beta_omega2%*%diff.y.tilde)
      der_omega2_q.f.y_tilde <- -2*num1_omega2/(omega2^3)+
        der_num1_omega2/(omega2^2)

      g[ind_omega2] <- (-0.5*(m/omega2+der_omega2_trace+
                                der_omega2_q.f.y-der_omega2_q.f.y_tilde))*omega2
    }

    if(n_re>0) {
      der_sigma2_re_trace <- list()
      sigma2_re_trace_aux <- list()
      der_sigma2_re_Sigma_tilde <- list()
      der_sigma2_re_Sigma_g <- list()
      der_Sigma_g_inv_sigma2_re_aux <- list()
      der_sigma2_re_aux <- list()
      M_beta_sigma2_re <- list()
      for(i in 1:n_re) {

        select_col <- sum(n_dim_re[1:i])
        der_Sigma_g_inv_sigma2_re_aux[[i]] <- matrix(0, nrow = sum(n_dim_re),
                                                     ncol = sum(n_dim_re))
        diag(der_Sigma_g_inv_sigma2_re_aux[[i]][select_col+1:n_dim_re[i+1],
                                                select_col+1:n_dim_re[i+1]])  <-
          -1/sigma2_re[i]^2
        der_sigma2_re_aux[[i]] <- Sigma_star_inv%*%der_Sigma_g_inv_sigma2_re_aux[[i]]

        M_beta_sigma2_re[[i]] <- der_sigma2_re_aux[[i]]%*%Sigma_star_inv

        der_sigma2_re_Sigma_g[[i]] <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
        diag(der_sigma2_re_Sigma_g[[i]][select_col+1:n_dim_re[i+1],
                                        select_col+1:n_dim_re[i+1]]) <- 1
        der_sigma2_re_Sigma_tilde[[i]] <- der_sigma2_re_Sigma_g[[i]]%*%C_g_m/omega2
        sigma2_re_trace_aux[[i]] <- Sigma_tilde_inv%*%der_sigma2_re_Sigma_tilde[[i]]
        der_sigma2_re_trace[[i]] <- sum(Matrix::diag(sigma2_re_trace_aux[[i]]))
        g[ind_sigma2_re[i]] <- (-0.5*der_sigma2_re_trace[[i]]-0.5*t(diff.y.tilde)%*%
                                  M_beta_sigma2_re[[i]]%*%
                                  diff.y.tilde/(omega2^2))*sigma2_re[i]

      }
    }
    return(g)
  }

  DtD <- t(D)%*%D



  hessian.log.lik <- function(par) {
    beta <- par[ind_beta]
    sigma2 <- exp(par[p+1])
    phi <- exp(par[ind_phi])
    if(!is.null(fix_tau2)) {
      nu2 <- fix_tau2/sigma2
    } else {
      nu2 <- exp(par[ind_nu2])
    }
    if(n_re>0) {
      sigma2_re <- exp(par[ind_sigma2_re])
    }
    if(!is.null(fix_var_me)) {
      omega2 <- fix_var_me
    } else {
      omega2 <- exp(par[ind_omega2])
    }
    n_p <- length(par)
    g <- rep(0, n_p)

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

    Sigma_g_C_g_m <- Sigma_g%*%C_g_m
    Sigma_tilde <- Sigma_g_C_g_m/omega2
    Matrix::diag(Sigma_tilde) <- Matrix::diag(Sigma_tilde) + 1
    Sigma_tilde_inv <- solve(Sigma_tilde)


    mu <- as.numeric(D%*%beta)
    diff.y <- y-mu
    diff.y.tilde <- as.numeric(Matrix::t(C_g)%*%diff.y)
    Sigma_star <- Sigma_g_inv+C_g_m/omega2
    Sigma_star_inv <- forceSymmetric(solve(Sigma_star))
    M_aux <- D.tilde%*%Sigma_star_inv

    H <- matrix(0, n_p, n_p)

    # beta - beta
    H[ind_beta, ind_beta] <- as.matrix(-DtD/omega2+
                                         +M_aux%*%Matrix::t(D.tilde)/(omega2^2))

    # beta - sigma2
    der_Sigma_g_inv_sigma2_aux <- matrix(0, nrow = sum(n_dim_re),
                                         ncol = sum(n_dim_re))
    der_Sigma_g_inv_sigma2_aux[1:n_dim_re[1], 1:n_dim_re[1]] <-
      -R.inv/sigma2^2
    der_sigma2_aux <- Sigma_star_inv%*%der_Sigma_g_inv_sigma2_aux
    M_beta_sigma2 <- der_sigma2_aux%*%Sigma_star_inv
    H[ind_beta, ind_sigma2] <-
      H[ind_sigma2, ind_beta] <- as.numeric(D.tilde%*%M_beta_sigma2%*%
                                              diff.y.tilde/(omega2^2))*sigma2

    # beta - phi

    # Derivatives for phi
    der_R_phi <- matrix(0, nrow = sum(n_dim_re),
                        ncol = sum(n_dim_re))
    M.der.phi <- matern.grad.phi(U, phi, kappa)
    der_R_phi[1:n_dim_re[1], 1:n_dim_re[1]] <-
      M.der.phi*sigma2
    der_Sigma_g_inv_phi_aux <- -Sigma_g_inv%*%der_R_phi%*%Sigma_g_inv
    der_phi_aux <- Sigma_star_inv%*%der_Sigma_g_inv_phi_aux
    M_beta_phi <- der_phi_aux%*%Sigma_star_inv

    H[ind_beta, ind_phi] <-
      H[ind_phi, ind_beta] <- as.numeric(D.tilde%*%M_beta_phi%*%
                                           diff.y.tilde/(omega2^2))*phi

    # beta - nu2
    if(is.null(fix_tau2)) {
      # Derivatives for nu2
      der_R_nu2 <- matrix(0, nrow = sum(n_dim_re),
                          ncol = sum(n_dim_re))
      diag(der_R_nu2[1:n_dim_re[1], 1:n_dim_re[1]]) <-
        sigma2
      der_Sigma_g_inv_nu2_aux <- -Sigma_g_inv%*%der_R_nu2%*%Sigma_g_inv
      der_nu2_aux <- Sigma_star_inv%*%der_Sigma_g_inv_nu2_aux
      M_beta_nu2 <- der_nu2_aux%*%Sigma_star_inv

      H[ind_beta, ind_nu2] <-
        H[ind_nu2, ind_beta] <- as.numeric(D.tilde%*%M_beta_nu2%*%
                                             diff.y.tilde/(omega2^2))*nu2
    }


    if(is.null(fix_var_me)) {
      # beta - omega2
      der_omega2_Sigma_star <- -C_g_m/omega2^2
      M_beta_omega2 <- -Sigma_star_inv%*%
        der_omega2_Sigma_star%*%
        Sigma_star_inv
      H[ind_beta, ind_omega2] <-
        H[ind_omega2, ind_beta] <-
        -(t(D)%*%diff.y/(omega2^2)+
            -2*as.numeric(D.tilde%*%Sigma_star_inv%*%diff.y.tilde/
                            (omega2^3))+
            as.numeric(D.tilde%*%M_beta_omega2%*%diff.y.tilde/
                         (omega2^2)))*omega2
    }

    # beta - sigma2_re
    if(n_re > 0) {
      M_beta_sigma2_re <- list()
      der_sigma2_re_aux <- list()
      der_Sigma_g_inv_sigma2_re_aux <- list()
      for(i in 1:n_re) {
        select_col <- sum(n_dim_re[1:i])

        der_Sigma_g_inv_sigma2_re_aux[[i]] <- matrix(0, nrow = sum(n_dim_re),
                                                     ncol = sum(n_dim_re))
        diag(der_Sigma_g_inv_sigma2_re_aux[[i]][select_col+1:n_dim_re[i+1],
                                                select_col+1:n_dim_re[i+1]])  <-
          -1/sigma2_re[i]^2
        der_sigma2_re_aux[[i]] <- Sigma_star_inv%*%der_Sigma_g_inv_sigma2_re_aux[[i]]

        M_beta_sigma2_re[[i]] <- der_sigma2_re_aux[[i]]%*%Sigma_star_inv
        H[ind_beta, ind_sigma2_re[i]] <-
          H[ind_sigma2_re[i], ind_beta] <-
          as.numeric(D.tilde%*%M_beta_sigma2_re[[i]]%*%
                       diff.y.tilde/(omega2^2))*sigma2_re[i]
      }
    }

    # sigma2 - sigma2
    # Derivatives for sigma2
    der2_Sigma_g_inv_sigma2_aux <- matrix(0, nrow = sum(n_dim_re),
                                          ncol = sum(n_dim_re))
    der2_Sigma_g_inv_sigma2_aux[1:n_dim_re[1], 1:n_dim_re[1]] <-
      2*R.inv/sigma2^3
    der_R_sigma2 <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
    der_R_sigma2[1:n_dim_re[1], 1:n_dim_re[1]] <- R
    der_sigma2_Sigma_g <- der_R_sigma2%*%C_g_m/omega2
    sigma2_trace_aux <- Sigma_tilde_inv%*%der_sigma2_Sigma_g
    der_sigma2_trace <- sum(Matrix::diag(sigma2_trace_aux))
    der2_sigma2_trace <- sum(Matrix::diag(-sigma2_trace_aux%*%sigma2_trace_aux))
    der.sigma2 <- (-0.5*der_sigma2_trace-0.5*t(diff.y.tilde)%*%
                     M_beta_sigma2%*%
                     diff.y.tilde/(omega2^2))*sigma2
    M2_sigma2 <- Sigma_star_inv%*%(2*der_Sigma_g_inv_sigma2_aux%*%
                                     der_sigma2_aux-der2_Sigma_g_inv_sigma2_aux)%*%
      Sigma_star_inv
    H[ind_sigma2, ind_sigma2] <-as.numeric(
      der.sigma2+
        (-0.5*der2_sigma2_trace+0.5*t(diff.y.tilde)%*%
           M2_sigma2%*%
           diff.y.tilde/(omega2^2))*sigma2^2)

    # sigma2 - phi
    der_R_sigma2_phi <- der_R_phi/sigma2
    der_phi_Sigma_g <- der_R_phi%*%C_g_m/omega2
    der_sigma2_phi_Sigma_g <- der_phi_Sigma_g/sigma2
    phi_trace_aux <- Sigma_tilde_inv%*%der_phi_Sigma_g
    der_phi_trace <- sum(Matrix::diag(phi_trace_aux))
    der_sigma2_phi_trace <- sum(Matrix::diag(Sigma_tilde_inv%*%der_sigma2_phi_Sigma_g-
                                               sigma2_trace_aux%*%phi_trace_aux))

    der_Sigma_g_inv_sigma2_phi_aux <-
      -Sigma_g_inv%*%(
        der_R_sigma2%*%Sigma_g_inv%*%der_R_phi+
          der_R_phi%*%Sigma_g_inv%*%der_R_sigma2-
          der_R_sigma2_phi
      )%*%
      Sigma_g_inv

    M2_sigma2_phi <- -Sigma_star_inv%*%(
      -der_Sigma_g_inv_sigma2_aux%*%der_phi_aux+
        -der_Sigma_g_inv_phi_aux%*%der_sigma2_aux-
        der_Sigma_g_inv_sigma2_phi_aux
    )%*%Sigma_star_inv

    H[ind_sigma2, ind_phi] <-
      H[ind_phi, ind_sigma2] <- as.numeric(
        (-0.5*der_sigma2_phi_trace+0.5*t(diff.y.tilde)%*%
           M2_sigma2_phi%*%
           diff.y.tilde/(omega2^2))*sigma2*phi)

    # sigma2 - nu2
    if(is.null(fix_tau2)) {
      der_R_sigma2_nu2 <- der_R_nu2/sigma2
      der_nu2_Sigma_g <- der_R_nu2%*%C_g_m/omega2
      der_sigma2_nu2_Sigma_g <- der_nu2_Sigma_g/sigma2
      nu2_trace_aux <- Sigma_tilde_inv%*%der_nu2_Sigma_g
      der_nu2_trace <- sum(Matrix::diag(nu2_trace_aux))
      der_sigma2_nu2_trace <- sum(Matrix::diag(Sigma_tilde_inv%*%der_sigma2_nu2_Sigma_g-
                                                 sigma2_trace_aux%*%nu2_trace_aux))

      der_Sigma_g_inv_sigma2_nu2_aux <-
        -Sigma_g_inv%*%(
          der_R_sigma2%*%Sigma_g_inv%*%der_R_nu2+
            der_R_nu2%*%Sigma_g_inv%*%der_R_sigma2-
            der_R_sigma2_nu2
        )%*%
        Sigma_g_inv

      M2_sigma2_nu2 <- -Sigma_star_inv%*%(
        -der_Sigma_g_inv_sigma2_aux%*%der_nu2_aux+
          -der_Sigma_g_inv_nu2_aux%*%der_sigma2_aux-
          der_Sigma_g_inv_sigma2_nu2_aux
      )%*%Sigma_star_inv

      H[ind_sigma2, ind_nu2] <-
        H[ind_nu2, ind_sigma2] <- as.numeric(
          (-0.5*der_sigma2_nu2_trace+0.5*t(diff.y.tilde)%*%
             M2_sigma2_nu2%*%
             diff.y.tilde/(omega2^2))*sigma2*nu2)
    }

    if(is.null(fix_var_me)) {
      # sigma2 - omega2
      M2_sigma2_omega2 <- -Sigma_star_inv%*%
        (der_omega2_Sigma_star%*%Sigma_star_inv%*%der_Sigma_g_inv_sigma2_aux+
           der_Sigma_g_inv_sigma2_aux%*%Sigma_star_inv%*%der_omega2_Sigma_star)%*%
        Sigma_star_inv
      der_omega2_Sigma_tilde <- -Sigma_g_C_g_m/omega2^2
      omega2_trace_aux <- Sigma_tilde_inv%*%der_omega2_Sigma_tilde
      der_sigma2_omega2_Sigma_g <- -der_sigma2_Sigma_g/omega2
      der_sigma2_omega2_trace <- sum(Matrix::diag(Sigma_tilde_inv%*%der_sigma2_omega2_Sigma_g-
                                                    sigma2_trace_aux%*%omega2_trace_aux))
      der_sigma2_omega2_q.f.y_tilde <- -2*as.numeric(t(diff.y.tilde)%*%
                                                       M_beta_sigma2%*%
                                                       diff.y.tilde/
                                                       (omega2^3))+
        as.numeric(t(diff.y.tilde)%*%
                     M2_sigma2_omega2%*%diff.y.tilde/
                     (omega2^2))

      H[ind_sigma2, ind_omega2] <-
        H[ind_omega2, ind_sigma2] <-
        (-0.5*(der_sigma2_omega2_trace+
                 der_sigma2_omega2_q.f.y_tilde))*omega2*sigma2
    }


    # sigma2 - sigma2_re
    if(n_re>0) {
      sigma2_re_trace_aux <- list()
      der_sigma2_re_Sigma_g <- list()
      der_sigma2_re_Sigma_tilde <- list()
      for(i in 1:n_re) {
        select_col <- sum(n_dim_re[1:i])


        der_sigma2_re_Sigma_g[[i]] <- matrix(0, nrow = sum(n_dim_re), ncol = sum(n_dim_re))
        diag(der_sigma2_re_Sigma_g[[i]][select_col+1:n_dim_re[i+1],
                                        select_col+1:n_dim_re[i+1]]) <- 1
        der_sigma2_re_Sigma_tilde[[i]] <- der_sigma2_re_Sigma_g[[i]]%*%C_g_m/omega2
        sigma2_re_trace_aux[[i]] <- Sigma_tilde_inv%*%der_sigma2_re_Sigma_tilde[[i]]
        der_sigma2_sigma2_re_trace <- sum(Matrix::diag(-sigma2_trace_aux%*%sigma2_re_trace_aux[[i]]))

        M2_sigma2_sigma2_re <- -Sigma_star_inv%*%(
          -der_Sigma_g_inv_sigma2_aux%*%der_sigma2_re_aux[[i]]+
            -der_Sigma_g_inv_sigma2_re_aux[[i]]%*%der_sigma2_aux
        )%*%Sigma_star_inv

        H[ind_sigma2, ind_sigma2_re[i]] <-
          H[ind_sigma2_re[i], ind_sigma2] <- as.numeric(
            (-0.5*der_sigma2_sigma2_re_trace+0.5*t(diff.y.tilde)%*%
               M2_sigma2_sigma2_re%*%
               diff.y.tilde/(omega2^2))*sigma2*sigma2_re[i])
      }
    }

    # phi - phi
    der2_R_phi <- matrix(0, nrow = sum(n_dim_re),
                         ncol = sum(n_dim_re))
    M.der2.phi <- matern.hessian.phi(U, phi, kappa)
    der2_R_phi[1:n_dim_re[1], 1:n_dim_re[1]] <-
      M.der2.phi*sigma2

    der2_Sigma_g_inv_phi_aux <- Sigma_g_inv%*%(
      2*der_R_phi%*%Sigma_g_inv%*%der_R_phi-
        der2_R_phi)%*%Sigma_g_inv
    der2_phi_Sigma_g <- der2_R_phi%*%C_g_m/omega2
    der2_phi_trace <- sum(Matrix::diag(Sigma_tilde_inv%*%der2_phi_Sigma_g
                                       -phi_trace_aux%*%phi_trace_aux))
    der.phi <- (-0.5*der_phi_trace-0.5*t(diff.y.tilde)%*%
                  M_beta_phi%*%
                  diff.y.tilde/(omega2^2))*phi
    M2_phi <- Sigma_star_inv%*%(2*der_Sigma_g_inv_phi_aux%*%
                                  der_phi_aux-
                                  der2_Sigma_g_inv_phi_aux)%*%
      Sigma_star_inv
    H[ind_phi, ind_phi] <-as.numeric(
      der.phi+
        (-0.5*der2_phi_trace+0.5*t(diff.y.tilde)%*%
           M2_phi%*%
           diff.y.tilde/(omega2^2))*phi^2)

    # phi - nu2
    if(is.null(fix_tau2)) {
      der_phi_nu2_trace <- sum(Matrix::diag(-phi_trace_aux%*%nu2_trace_aux))

      der_Sigma_g_inv_phi_nu2_aux <-
        -Sigma_g_inv%*%(
          der_R_phi%*%Sigma_g_inv%*%der_R_nu2+
            der_R_nu2%*%Sigma_g_inv%*%der_R_phi
        )%*%
        Sigma_g_inv

      M2_phi_nu2 <- -Sigma_star_inv%*%(
        -der_Sigma_g_inv_phi_aux%*%der_nu2_aux+
          -der_Sigma_g_inv_nu2_aux%*%der_phi_aux-
          der_Sigma_g_inv_phi_nu2_aux
      )%*%Sigma_star_inv

      H[ind_phi, ind_nu2] <-
        H[ind_nu2, ind_phi] <- as.numeric(
          (-0.5*der_phi_nu2_trace+0.5*t(diff.y.tilde)%*%
             M2_phi_nu2%*%
             diff.y.tilde/(omega2^2))*phi*nu2)
    }

    if(is.null(fix_var_me)) {
      # phi - omega2
      M2_phi_omega2 <- -Sigma_star_inv%*%
        (der_omega2_Sigma_star%*%Sigma_star_inv%*%der_Sigma_g_inv_phi_aux+
           der_Sigma_g_inv_phi_aux%*%Sigma_star_inv%*%der_omega2_Sigma_star)%*%
        Sigma_star_inv
      der_phi_omega2_Sigma_g <- -der_phi_Sigma_g/omega2
      der_phi_omega2_trace <- sum(Matrix::diag(Sigma_tilde_inv%*%der_phi_omega2_Sigma_g-
                                                 phi_trace_aux%*%omega2_trace_aux))
      der_phi_omega2_q.f.y_tilde <- -2*as.numeric(t(diff.y.tilde)%*%
                                                    M_beta_phi%*%
                                                    diff.y.tilde/
                                                    (omega2^3))+
        as.numeric(t(diff.y.tilde)%*%
                     M2_phi_omega2%*%diff.y.tilde/
                     (omega2^2))

      H[ind_phi, ind_omega2] <-
        H[ind_omega2, ind_phi] <-
        (-0.5*(der_phi_omega2_trace+
                 der_phi_omega2_q.f.y_tilde))*omega2*phi
    }


    #phi - sigma2_re
    if(n_re>0) {
      for(i in 1:n_re) {
        der_phi_sigma2_re_trace <- sum(Matrix::diag(-phi_trace_aux%*%sigma2_re_trace_aux[[i]]))

        M2_phi_sigma2_re <- -Sigma_star_inv%*%(
          -der_Sigma_g_inv_phi_aux%*%der_sigma2_re_aux[[i]]+
            -der_Sigma_g_inv_sigma2_re_aux[[i]]%*%der_phi_aux
        )%*%Sigma_star_inv

        H[ind_phi, ind_sigma2_re[i]] <-
          H[ind_sigma2_re[i], ind_phi] <- as.numeric(
            (-0.5*der_phi_sigma2_re_trace+0.5*t(diff.y.tilde)%*%
               M2_phi_sigma2_re%*%
               diff.y.tilde/(omega2^2))*phi*sigma2_re[i])
      }
    }

    if(is.null(fix_tau2)) {
      # nu2 - nu2
      der2_Sigma_g_inv_nu2_aux <- Sigma_g_inv%*%(
        2*der_R_nu2%*%Sigma_g_inv%*%der_R_nu2)%*%Sigma_g_inv

      der2_nu2_trace <- sum(Matrix::diag(-nu2_trace_aux%*%nu2_trace_aux))
      der.nu2 <- (-0.5*der_nu2_trace-0.5*t(diff.y.tilde)%*%
                    M_beta_nu2%*%
                    diff.y.tilde/(omega2^2))*nu2
      M2_nu2 <- Sigma_star_inv%*%(2*der_Sigma_g_inv_nu2_aux%*%
                                    der_nu2_aux-
                                    der2_Sigma_g_inv_nu2_aux)%*%
        Sigma_star_inv
      H[ind_nu2, ind_nu2] <-as.numeric(
        der.nu2+
          (-0.5*der2_nu2_trace+0.5*t(diff.y.tilde)%*%
             M2_nu2%*%
             diff.y.tilde/(omega2^2))*nu2^2)

      if(is.null(fix_var_me)) {
        # nu2 - omega2
        M2_nu2_omega2 <- -Sigma_star_inv%*%
          (der_omega2_Sigma_star%*%Sigma_star_inv%*%der_Sigma_g_inv_nu2_aux+
             der_Sigma_g_inv_nu2_aux%*%Sigma_star_inv%*%der_omega2_Sigma_star)%*%
          Sigma_star_inv
        der_nu2_omega2_Sigma_g <- -der_nu2_Sigma_g/omega2
        der_nu2_omega2_trace <- sum(Matrix::diag(Sigma_tilde_inv%*%der_nu2_omega2_Sigma_g-
                                                   nu2_trace_aux%*%omega2_trace_aux))
        der_nu2_omega2_q.f.y_tilde <- -2*as.numeric(t(diff.y.tilde)%*%
                                                      M_beta_nu2%*%
                                                      diff.y.tilde/
                                                      (omega2^3))+
          as.numeric(t(diff.y.tilde)%*%
                       M2_nu2_omega2%*%diff.y.tilde/
                       (omega2^2))

        H[ind_nu2, ind_omega2] <-
          H[ind_omega2, ind_nu2] <-
          (-0.5*(der_nu2_omega2_trace+
                   der_nu2_omega2_q.f.y_tilde))*omega2*nu2
      }


      #nu2 - sigma2_re
      if(n_re>0) {
        for(i in 1:n_re) {
          der_nu2_sigma2_re_trace <- sum(Matrix::diag(-nu2_trace_aux%*%sigma2_re_trace_aux[[i]]))

          M2_nu2_sigma2_re <- -Sigma_star_inv%*%(
            -der_Sigma_g_inv_nu2_aux%*%der_sigma2_re_aux[[i]]+
              -der_Sigma_g_inv_sigma2_re_aux[[i]]%*%der_nu2_aux
          )%*%Sigma_star_inv

          H[ind_nu2, ind_sigma2_re[i]] <-
            H[ind_sigma2_re[i], ind_nu2] <- as.numeric(
              (-0.5*der_nu2_sigma2_re_trace+0.5*t(diff.y.tilde)%*%
                 M2_nu2_sigma2_re%*%
                 diff.y.tilde/(omega2^2))*nu2*sigma2_re[i])
        }
      }
    }

    if(is.null(fix_var_me)) {
      #omega2 - omega2
      der_omega2_q.f.y <- -as.numeric(sum(diff.y^2)/omega2^2)
      der2_omega2_q.f.y <- 2*as.numeric(sum(diff.y^2)/omega2^3)

      omega2_trace_aux <- Sigma_tilde_inv%*%der_omega2_Sigma_tilde
      der_omega2_trace <- sum(Matrix::diag(omega2_trace_aux))
      der2_omega2_Sigma_tilde <- 2*Sigma_g_C_g_m/omega2^3
      der2_omega2_trace <- sum(Matrix::diag(Sigma_tilde_inv%*%der2_omega2_Sigma_tilde-
                                              omega2_trace_aux%*%omega2_trace_aux))

      num1_omega2 <- as.numeric(t(diff.y.tilde)%*%Sigma_star_inv%*%diff.y.tilde)
      der_num1_omega2 <- as.numeric(t(diff.y.tilde)%*%
                                      M_beta_omega2%*%diff.y.tilde)
      der2_omega2_Sigma_star <- 2*C_g_m/omega2^3
      M2_beta_omega2 <- Sigma_star_inv%*%
        (2*der_omega2_Sigma_star%*%Sigma_star_inv%*%der_omega2_Sigma_star-
           der2_omega2_Sigma_star)%*%
        Sigma_star_inv
      der2_num1_omega2 <- as.numeric(t(diff.y.tilde)%*%
                                       M2_beta_omega2%*%diff.y.tilde)

      der_omega2_q.f.y_tilde <- -2*num1_omega2/(omega2^3)+
        der_num1_omega2/(omega2^2)
      der2_omega2_q.f.y_tilde <- -2*(der_num1_omega2*omega2-3*num1_omega2)/
        (omega2^4)+
        (der2_num1_omega2*omega2-
           2*der_num1_omega2)/(omega2^3)

      #     g[ind_omega2] <- (-0.5*(m/omega2+der_omega2_trace+
      #     der_omega2_q.f.y-der_omega2_q.f.y_tilde))*omega2

      der.omega2 <- (-0.5*(m/omega2+der_omega2_trace+
                             der_omega2_q.f.y-der_omega2_q.f.y_tilde))*omega2

      H[ind_omega2, ind_omega2] <- der.omega2+
        -0.5*(-m/omega2^2+der2_omega2_trace+
                der2_omega2_q.f.y-der2_omega2_q.f.y_tilde)*(omega2^2)

      #omega2 - sigma2_re
      if(n_re>0) {
        for(i in 1:n_re) {
          der_omega2_sigma2_re_Sigma_g <- -der_sigma2_re_Sigma_tilde[[i]]/omega2
          der_omega2_sigma2_re_trace <- sum(Matrix::diag(Sigma_tilde_inv%*%der_omega2_sigma2_re_Sigma_g-
                                                           omega2_trace_aux%*%sigma2_re_trace_aux[[i]]))

          M2_omega2_sigma2_re <- -Sigma_star_inv%*%
            (der_omega2_Sigma_star%*%Sigma_star_inv%*% der_Sigma_g_inv_sigma2_re_aux[[i]]+
               der_Sigma_g_inv_sigma2_re_aux[[i]]%*%Sigma_star_inv%*%der_omega2_Sigma_star)%*%
            Sigma_star_inv

          der_sigma2_omega2_q.f.y_tilde <- -2*as.numeric(t(diff.y.tilde)%*%
                                                           M_beta_sigma2_re[[i]]%*%
                                                           diff.y.tilde/
                                                           (omega2^3))+
            as.numeric(t(diff.y.tilde)%*%
                         M2_omega2_sigma2_re%*%diff.y.tilde/
                         (omega2^2))

          H[ind_omega2, ind_sigma2_re[i]] <-
            H[ind_sigma2_re[i], ind_omega2] <- as.numeric(
              -0.5*(der_omega2_sigma2_re_trace+der_sigma2_omega2_q.f.y_tilde)*
                omega2*sigma2_re[i])

        }
      }
    }
    # sigma2_re - sigma2_re
    if(n_re > 0) {
      der_sigma2_re_trace <- list()
      der.sigma2_re <- list()
      M2_sigma2_re <- list()
      der2_Sigma_g_inv_sigma2_re_aux <- list()
      der2_sigma2_re_trace <- list()
      for(i in 1:n_re) {
        select_col <- sum(n_dim_re[1:i])
        der2_Sigma_g_inv_sigma2_re_aux[[i]] <- matrix(0, nrow = sum(n_dim_re),
                                                      ncol = sum(n_dim_re))
        diag(der2_Sigma_g_inv_sigma2_re_aux[[i]][select_col+1:n_dim_re[i+1],
                                                 select_col+1:n_dim_re[i+1]]) <-
          2/sigma2_re[i]^3
        der_sigma2_re_trace[[i]] <- sum(Matrix::diag(sigma2_re_trace_aux[[i]]))
        der.sigma2_re[[i]] <- (-0.5*der_sigma2_re_trace[[i]]-0.5*t(diff.y.tilde)%*%
                                 M_beta_sigma2_re[[i]]%*%
                                 diff.y.tilde/(omega2^2))*sigma2_re[i]
        M2_sigma2_re[[i]] <- Sigma_star_inv%*%(
          2*der_Sigma_g_inv_sigma2_re_aux[[i]]%*%
            der_sigma2_re_aux[[i]]-der2_Sigma_g_inv_sigma2_re_aux[[i]])%*%
          Sigma_star_inv
        der2_sigma2_re_trace[[i]] <- sum(Matrix::diag(-sigma2_re_trace_aux[[i]]%*%
                                                        sigma2_re_trace_aux[[i]]))
        H[ind_sigma2_re[i], ind_sigma2_re[i]] <-as.numeric(
          der.sigma2_re[[i]]+
            (-0.5*der2_sigma2_re_trace[[i]]+0.5*t(diff.y.tilde)%*%
               M2_sigma2_re[[i]]%*%
               diff.y.tilde/(omega2^2))*sigma2_re[i]^2)
        if(i < n_re) {
          for(j in (i+1):n_re) {

            der_sigma2_re_ij_trace <- sum(Matrix::diag(-sigma2_re_trace_aux[[i]]%*%
                                                         sigma2_re_trace_aux[[j]]))

            M2_sigma2_re_ij <- -Sigma_star_inv%*%(
              -der_Sigma_g_inv_sigma2_re_aux[[i]]%*%der_sigma2_re_aux[[j]]+
                -der_Sigma_g_inv_sigma2_re_aux[[j]]%*%der_sigma2_re_aux[[i]]
            )%*%Sigma_star_inv

            H[ind_sigma2_re[i], ind_sigma2_re[j]] <-
              H[ind_sigma2_re[j], ind_sigma2_re[i]] <- as.numeric(
                (-0.5*der_sigma2_re_ij_trace+0.5*t(diff.y.tilde)%*%
                   M2_sigma2_re_ij%*%
                   diff.y.tilde/(omega2^2))*sigma2_re[i]*sigma2_re[j])
          }
        }
      }
    }


    return(H)
  }

  start_cov_pars[-(1:2)] <- start_cov_pars[-(1:2)]/start_cov_pars[1]
  start_par <- c(start_beta, log(start_cov_pars))

  out <- list()
  estim <- nlminb(start_par,
                  function(x) -log.lik(x),
                  function(x) -grad.log.lik(x),
                  function(x) -hessian.log.lik(x),
                  control=list(trace=1*messages))

  out$estimate <- estim$par
  out$grad.MLE <- grad.log.lik(estim$par)
  hess.MLE <- hessian.log.lik(estim$par)
  out$covariance <- solve(-hess.MLE)
  out$log.lik <- -estim$objective


  class(out) <- "RiskMap"
  return(out)
}


##' Simulation from Generalized Linear Gaussian Process Models
##'
##' Simulates data from a fitted Generalized Linear Gaussian Process Model (GLGPM) or a specified model formula and data.
##'
##' @param n_sim Number of simulations to perform.
##' @param model_fit Fitted GLGPM model object of class 'RiskMap'. If provided, overrides 'formula', 'data', 'family', 'crs', 'convert_to_crs', 'scale_to_km', and 'control_mcmc' arguments.
##' @param formula Model formula indicating the variables of the model to be simulated.
##' @param data Data frame or 'sf' object containing the variables in the model formula.
##' @param family Distribution family for the response variable. Must be one of 'gaussian', 'binomial', or 'poisson'.
##' @param den Required for 'binomial' to denote the denominator (i.e. number of trials) of the Binomial distribution.
##' For the 'poisson' family, the argument is optional and is used a multiplicative term to express the mean counts.
##' @param cov_offset Offset for the covariate part of the GLGPM.
##' @param crs Coordinate reference system (CRS) code for spatial data.
##' @param convert_to_crs CRS code to convert spatial data if different from 'crs'.
##' @param scale_to_km Logical; if TRUE, distances between locations are computed in kilometers; if FALSE, in meters.
##' @param control_mcmc Control parameters for MCMC simulation if applicable.
##' @param sim_pars List of simulation parameters including 'beta', 'sigma2', 'tau2', 'phi', 'sigma2_me', and 'sigma2_re'.
##' @param messages Logical; if TRUE, display progress and informative messages.
##'
##' @details
##' Generalized Linear Gaussian Process Models (GLGPMs) extend generalized linear models (GLMs) by incorporating spatial Gaussian processes to model spatial correlation. This function simulates data from GLGPMs using Markov Chain Monte Carlo (MCMC) methods. It supports Gaussian, binomial, and Poisson response families, utilizing a Matern correlation function to model spatial dependence.
##'
##' The simulation process involves generating spatially correlated random effects and simulating responses based on the fitted or specified model parameters. For 'gaussian' family, the function simulates response values by adding measurement error.
##'
##' Additionally, GLGPMs can incorporate unstructured random effects specified through the \code{re()} term in the model formula, allowing for capturing additional variability beyond fixed and spatial covariate effects.
##'
##' @return A list containing simulated data, simulated spatial random effects (if applicable), and other simulation parameters.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##' @export
glgpm_sim <- function(n_sim,
                      model_fit = NULL,
                      formula = NULL,
                      data = NULL,
                      family = NULL,
                      den = NULL,
                      cov_offset = NULL,
                      crs = NULL, convert_to_crs = NULL,
                      scale_to_km = TRUE,
                      control_mcmc = NULL,
                      sim_pars = list(beta = NULL,
                                      sigma2 = NULL,
                                      tau2 = NULL,
                                      phi = NULL,
                                      sigma2_me = NULL,
                                      sigma2_re = NULL),
                      messages = TRUE) {

  if(!is.null(model_fit)) {
    if(!inherits(model_fit,
                    what = "RiskMap", which = FALSE)) stop("'model_fit' must be of class 'RiskMap'")
    formula <- as.formula(model_fit$formula)
    data <- model_fit$data_sf
    family = model_fit$family
    crs <- model_fit$crs
    convert_to_crs <- model_fit$convert_to_crs
    scale_to_km <- model_fit$scale_to_km
  }
  inter_f <- interpret.formula(formula)

  if(family=="binomial" | family=="poisson") {
    if(is.null(control_mcmc)) stop("if family='binomial' or family='poisson'
                                   'control_mcmc' must be provided")
  }

  if(!inherits(formula,
               what = "formula", which = FALSE)) {
    stop("'formula' must be a 'formula'
                                     object indicating the variables of the
                                     model to be fitted")
  }


  if(length(crs)>0) {
    if(!is.numeric(crs) |
       (is.numeric(crs) &
        (crs%%1!=0 | crs <0))) stop("'crs' must be a positive integer number")
  }
  if(class(data)[1]=="data.frame") {
    if(is.null(crs)) {
      warning("'crs' is set to 4326 (long/lat)")
      crs <- 4326
    }
    if(length(inter_f$gp.spec$term)==2) {
      new_x <- paste(inter_f$gp.spec$term[1],"_sf",sep="")
      new_y <- paste(inter_f$gp.spec$term[2],"_sf",sep="")
      data[[new_x]] <-  data[[inter_f$gp.spec$term[1]]]
      data[[new_y]] <-  data[[inter_f$gp.spec$term[2]]]
      data <- st_as_sf(data,
                       coords = c(new_x, new_y),
                       crs = crs)
    }
  }

  if(length(inter_f$gp.spec$term) == 1 & inter_f$gp.spec$term[1]=="sf" &
     class(data)[1]!="sf") stop("'data' must be an object of class 'sf'")


  if(class(data)[1]=="sf") {
    if(is.na(st_crs(data)) & is.null(crs)) {
      stop("the CRS of the sf object passed to 'data' is missing and and is not specified through 'crs'")
    } else if(is.na(st_crs(data))) {
      data <- st_as_sf(data, crs = crs)
    }
  }

  kappa <- inter_f$gp.spec$kappa
  if(kappa < 0) stop("kappa must be positive.")

  if(family != "gaussian" & family != "binomial" &
     family != "poisson") stop("'family' must be either 'gaussian', 'binomial'
                               or 'poisson'")


  mf <- model.frame(inter_f$pf,data=data, na.action = na.fail)


  # Extract covariates matrix
  D <- as.matrix(model.matrix(attr(mf,"terms"),data=data))
  n <- nrow(D)

  if(length(inter_f$re.spec) > 0) {
    hr_re <- inter_f$re.spec$term
  } else {
    hr_re <- NULL
  }


  if(!is.null(hr_re)) {
    # Define indices of random effects
    re_mf <- st_drop_geometry(data[hr_re])
    re_mf_n <- re_mf

    if(any(is.na(re_mf))) stop("Missing values in the variable(s) of the random effects specified through re() ")
    names_re <- colnames(re_mf)
    n_re <- ncol(re_mf)

    ID_re <- matrix(NA, nrow = n, ncol = n_re)
    re_unique <- list()
    re_unique_f <- list()
    for(i in 1:n_re) {
      if(is.factor(re_mf[,i])) {
        re_mf_n[,i] <- as.numeric(re_mf[,i])
        re_unique[[names_re[i]]] <- 1:length(levels(re_mf[,i]))
        ID_re[, i] <- sapply(1:n,
                             function(j) which(re_mf_n[j,i]==re_unique[[names_re[i]]]))
        re_unique_f[[names_re[i]]] <-levels(re_mf[,i])
      } else if(is.numeric(re_mf[,i])) {
        re_unique[[names_re[i]]] <- unique(re_mf[,i])
        ID_re[, i] <- sapply(1:n,
                             function(j) which(re_mf_n[j,i]==re_unique[[names_re[i]]]))
        re_unique_f[[names_re[i]]] <- re_unique[[names_re[i]]]
      }
    }
  } else {
    n_re <- 0
    re_unique <- NULL
    ID_re <- NULL
  }

  # Number of covariates
  p <- ncol(D)

  if(!is.null(model_fit)) {
    ind_beta <- 1:p
    par_hat <- coef(model_fit)

    beta <- par_hat[ind_beta]
    ind_sigma2 <- p+1
    sigma2 <- par_hat[ind_sigma2]
    ind_phi <- p+2
    phi <- par_hat[ind_phi]

    if(is.null(model_fit$fix_tau2)) {
      ind_tau2 <- p+3
      tau2 <- par_hat[ind_tau2]
      if(is.null(model_fit$fix_var_me)) {
        ind_sigma2_me <- p+4
        sigma2_me <- par_hat[ind_sigma2_me]
      } else {
        sigma2_me <- model_fit$fix_var_me
      }
      if(n_re>0) {
        ind_sigma2_re <- (p+5):(p+4+n_re)
        sigma2_re <- par_hat[ind_sigma2_re]
      }
    } else {
      tau2 <- model_fit$fix_tau2
      if(is.null(model_fit$fix_var_me)) {
        ind_sigma2_me <- p+3
        sigma2_me <- par_hat[ind_sigma2_me]
      } else {
        sigma2_me <- model_fit$fix_var_me
      }
      if(n_re>0) {
        ind_sigma2_re <- (p+4):(p+3+n_re)
        sigma2_re <- par_hat[ind_sigma2_re]
      }
    }
  } else {
    if(is.null(sim_pars$beta)) stop("'beta' is missing")
    beta <- sim_pars$beta
    if(length(beta)!=p) stop("the number of values provided for 'beta' does not match
    the number of covariates specified in the formula")
    if(is.null(sim_pars$sigma2)) stop("'sigma2' is missing")
    sigma2 <- sim_pars$sigma2
    if(is.null(sim_pars$phi)) stop("'phi' is missing")
    phi <- sim_pars$phi
    if(is.null(sim_pars$tau2)) stop("'tau2' is missing")
    tau2 <- sim_pars$tau2
    if(is.null(sim_pars$sigma2_me)) stop("'sigma2_me' is missing")
    sigma2_me <- sim_pars$sigma2_me
    if(n_re>0) {
      if(is.null(sim_pars$sigma2_re)) stop("'sigma2_re' is missing")
      if(length(sim_pars$sigma2_re)!=n_re) stop("the values passed to 'sigma2_re' in 'sim_pars'
      does not match the number of random effects specfied in re() in the formula")
      sigma2_re <- sim_pars$sigma2_re
    }
  }

  # Extract coordinates
  if(!is.null(convert_to_crs)) {
    if(!is.numeric(convert_to_crs)) stop("'convert_to_utm' must be a numeric object")
    data <- st_transform(data, crs = convert_to_crs)
    crs <- convert_to_crs
  }
  if(messages) message("The CRS used is", as.list(st_crs(data))$input, "\n")

  coords_o <- st_coordinates(data)
  coords <- unique(coords_o)

  m <- nrow(coords_o)
  ID_coords <- sapply(1:m, function(i)
    which(coords_o[i,1]==coords[,1] &
            coords_o[i,2]==coords[,2]))
  s_unique <- unique(ID_coords)




  if(all(table(ID_coords)==1) & (tau2!=0 & sigma2_me!=0)) {
    stop("When there is only one observation per location, both the nugget and measurement error cannot
         be estimate. Consider removing either one of them. ")
  }

  if(scale_to_km) {
    coords_o <- coords_o/1000
    coords <- coords/1000
    if(messages) message("Distances between locations are computed in kilometers \n")
  } else {
    if(messages) message("Distances between locations are computed in meters \n")
  }

  # Simulate S
  Sigma <- sigma2*matern_cor(dist(coords), phi = phi, kappa = kappa,
                             return_sym_matrix = TRUE)
  diag(Sigma) <- diag(Sigma) + tau2
  Sigma_sroot <- t(chol(Sigma))
  S_sim <- t(sapply(1:n_sim, function(i) Sigma_sroot%*%rnorm(nrow(coords))))

  # Simulate random effects
  if(n_re>0) {
    re_sim <- list()
    if(!is.null(model_fit)) {
      re_names <- names(model_fit$re)
    } else {
      re_names <- inter_f$re.spec$term
    }

    dim_re <- sapply(1:n_re, function(j) length(re_unique[[j]]))
    for(i in 1:n_sim) {
      re_sim[[i]] <- list()
      for(j in 1:n_re) {
        re_sim[[i]][[paste(re_names[j])]] <- rnorm(dim_re[j])*sqrt(sigma2_re[j])
      }
    }
  }

  # Linear predictor
  eta_sim <- t(sapply(1:n_sim, function(i) D%*%beta + S_sim[i,][ID_coords]))

  if(n_re > 0) {
    for(i in 1:n_sim) {
      for(j in 1:n_re) {
        eta_sim[i,] <- eta_sim[i,] + re_sim[[i]][[paste(re_names[j])]][ID_re[,j]]
      }
    }
  }

  y_sim <- matrix(NA, nrow=n_sim, ncol=n)
  if(family=="gaussian") {
    lin_pred <- eta_sim

    for(i in 1:n_sim) {
      y_sim[i,] <- lin_pred[i,] + sqrt(sigma2_me)*rnorm(n)
    }
  }

  if(!is.null(model_fit)) {
    data_sim <- model_fit$data_sf
  } else {
    data_sim <- data
  }

  for(i in 1:n_sim) {
    data_sim[[paste(inter_f$response,"_sim",i,sep="")]] <- y_sim[i,]
  }
  out <- list(data_sim = data_sim,
              S_sim = S_sim,
              lin_pred_sim = lin_pred,
              beta = beta,
              sigma2 = sigma2,
              tau2 = tau2,
              phi = phi)
  if(family=="gaussian") {
    out$sigma2_me <- sigma2_me
  }
  if(n_re>0) {
    out$sigma2_re <- sigma2_re
    out$re_sim <- re_sim
  }
  return(out)
}

##' Maximization of the Integrand for Generalized Linear Gaussian Process Models
##'
##' Maximizes the integrand function for Generalized Linear Gaussian Process Models (GLGPMs), which involves the evaluation of likelihood functions with spatially correlated random effects.
##'
##' @param y Response variable vector.
##' @param units_m Units of measurement for the response variable.
##' @param mu Mean vector of the response variable.
##' @param Sigma Covariance matrix of the spatial process.
##' @param ID_coords Indices mapping response to locations.
##' @param ID_re Indices mapping response to unstructured random effects.
##' @param family Distribution family for the response variable. Must be one of 'gaussian', 'binomial', or 'poisson'.
##' @param sigma2_re Variance of the unstructured random effects.
##' @param hessian Logical; if TRUE, compute the Hessian matrix.
##' @param gradient Logical; if TRUE, compute the gradient vector.
##'
##' @details
##' This function maximizes the integrand for GLGPMs using the Nelder-Mead optimization algorithm. It computes the likelihood function incorporating spatial covariance and unstructured random effects, if provided.
##'
##' The integrand includes terms for the spatial process (Sigma), unstructured random effects (sigma2_re), and the likelihood function (llik) based on the specified distribution family ('gaussian', 'binomial', or 'poisson').
##'
##' @return A list containing the mode estimate, and optionally, the Hessian matrix and gradient vector.
##' @export
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
maxim.integrand <- function(y,units_m,mu,Sigma,ID_coords, ID_re = NULL,family,
                            sigma2_re = NULL,
                            hessian=FALSE, gradient=FALSE) {

    Sigma.inv <- solve(Sigma)
    n_loc <- nrow(Sigma)
    n <- length(y)
    if((!is.null(ID_re) & is.null(sigma2_re)) | (is.null(ID_re) & !is.null(sigma2_re))) {
      stop("To introduce unstructured random effects both `ID_re` and `sigma2_re`
           must be provided.")
    }
    n_re <- length(sigma2_re)
    if(n_re > 0) {
      n_dim_re <- sapply(1:n_re, function(i) length(unique(ID_re[,i])))
      ind_re <- list()
      add_i <- 0
      for(i in 1:n_re) {
        ind_re[[i]] <- (add_i+n_loc+1):(add_i+n_loc+n_dim_re[i])
        if(i < n_re) add_i <- sum(n_dim_re[1:i])
      }
    }
    n_tot <- n_loc
    if(n_re > 0) n_tot <- n_tot + sum(n_dim_re)

    integrand <- function(S_tot) {
      S <- S_tot[1:n_loc]

      q.f_S <- as.numeric(t(S)%*%Sigma.inv%*%(S))

      q.f_re <- 0
      if(n_re > 0) {
        S_re <- NULL
        S_re_list <- list()
        for(i in 1:n_re) {
          S_re_list[[i]] <- S_tot[ind_re[[i]]]
          q.f_re <- q.f_re + sum(S_re_list[[i]]^2)/sigma2_re[i]
        }
      }

      eta <- mu + S[ID_coords]
      if(n_re > 0) {
        for(i in 1:n_re) {
          eta <- eta + S_re_list[[i]][ID_re[,i]]
        }
      }

      if(family=="poisson") {
        llik <- sum(y*eta-units_m*exp(eta))
      } else if(family=="binomial") {
        llik <- sum(y*eta-units_m*log(1+exp(eta)))
      }

    out <- -0.5*q.f_S-0.5*q.f_re+llik
    return(out)
  }


  C_S <- t(sapply(1:n_loc,function(i) ID_coords==i))

  if(n_re>0) {
    C_re <- list()
    C_S_re <- list()
    C_re_re <- list()
    for(j in 1:n_re){
      C_S_re[[j]] <- array(FALSE,dim = c(n_loc, n_dim_re[j], n))
      C_re[[j]] <- t(sapply(1:n_dim_re[j],function(i) ID_re[,j]==i))
      for(l in 1:n_dim_re[j]) {
        for(k in 1:n_loc) {
          ind_kl <- which(ID_coords==k & ID_re[,j]==l)
          if(length(ind_kl) > 0) {
            C_S_re[[j]][k,l,ind_kl] <- TRUE
          }
        }
      }

      if(j < n_re) {
        C_re_re[[j]] <- list()
        counter <- 0
        for(w in (j+1):n_re) {
          counter <- counter+1
          C_re_re[[j]][[counter]] <- array(FALSE,dim = c(n_dim_re[j], n_dim_re[w], n))
          for(l in 1:n_dim_re[j]) {
            for(k in 1:n_dim_re[w]) {
              ind_lk <- which(ID_re[,j]==l & ID_re[,w]==k)
              if(length(ind_kl) > 0) {
                C_re_re[[j]][[counter]][l,k,ind_lk] <- TRUE
              }
            }
          }
        }
      }
    }
  }

  grad.integrand <- function(S_tot) {
    S <- S_tot[1:n_loc]

    if(n_re > 0) {
      S_re_list <- list()
      for(i in 1:n_re) {
        S_re_list[[i]] <- S_tot[ind_re[[i]]]
      }
    }

    eta <- mu + S[ID_coords]
    if(n_re > 0) {
      for(i in 1:n_re) {
        eta <- eta + S_re_list[[i]][ID_re[,i]]
      }
    }

    if(family=="poisson") {
      h <- units_m*exp(eta)
    } else if(family=="binomial") {
      h <- units_m*exp(eta)/(1+exp(eta))
    }

    out <- rep(NA,n_tot)
    out[1:n_loc] <- as.numeric(-Sigma.inv%*%S+
                              sapply(1:n_loc,function(i) sum((y-h)[C_S[i,]])))
    if(n_re>0) {
      for(j in 1:n_re) {
        out[ind_re[[j]]] <- as.numeric(-S_re_list[[j]]/sigma2_re[[j]]+
                                         sapply(1:n_dim_re[[j]],
                                          function(x) sum((y-h)[C_re[[j]][x,]])))
      }
    }
    return(out)
  }

  hessian.integrand <- function(S_tot) {
    S <- S_tot[1:n_loc]


    if(n_re > 0) {
      S_re <- NULL
      S_re_list <- list()
      for(i in 1:n_re) {
        S_re_list[[i]] <- S_tot[ind_re[[i]]]
      }
    }

    eta <- mu + S[ID_coords]
    if(n_re > 0) {
      for(i in 1:n_re) {
        eta <- eta + S_re_list[[i]][ID_re[,i]]
      }
    }

    if(family=="poisson") {
      h <- units_m*exp(eta)
      h1 <- h
    } else if(family=="binomial") {
      h <- units_m*exp(eta)/(1+exp(eta))
      h1 <- h/(1+exp(eta))

    }

    out <- matrix(0,nrow = n_tot, ncol = n_tot)

    out[1:n_loc, 1:n_loc] <-  -Sigma.inv
    diag(out[1:n_loc, 1:n_loc]) <- diag(out[1:n_loc, 1:n_loc])+
                                  -sapply(1:n_loc,function(i) sum(h1[C_S[i,]]))
    if(n_re>0) {
      for(j in 1:n_re) {
        diag(out[ind_re[[j]], ind_re[[j]]]) <- -1/sigma2_re[j]
        diag(out[ind_re[[j]], ind_re[[j]]]) <- diag(out[ind_re[[j]], ind_re[[j]]])+
          -sapply(1:n_dim_re[j],function(i) sum(h1[C_re[[j]][i,]]))

        out[1:n_loc,ind_re[[j]]]

        for(k in 1:n_dim_re[[j]]) {
          out[1:n_loc, ind_re[[j]]][,k] <- -sapply(1:n_loc,function(i) sum(h1[C_S_re[[j]][i,k,]]))
          out[ind_re[[j]], 1:n_loc][k,] <- out[1:n_loc,ind_re[[j]]][,k]
        }

        if(j < n_re) {
          counter <- 0
          for(w in (j+1):n_re) {
            counter <- counter + 1
            for(k in 1:n_dim_re[[w]]) {
              out[ind_re[[j]], ind_re[[w]]][,k] <- -sapply(1:n_dim_re[j],function(i) sum(h1[C_re_re[[j]][[counter]][i,k,]]))
              out[ind_re[[w]], ind_re[[j]]][k,] <- out[ind_re[[j]], ind_re[[w]]][,k]
            }
          }
        }
      }
    }
    return(out)
  }


  estim <- nlminb(rep(0,n_tot),
                  function(x) -integrand(x),
                  function(x) -grad.integrand(x),
                  function(x) -hessian.integrand(x))


  out <- list()
  out$mode <- estim$par
  if(hessian) {
    out$hessian <- hessian.integrand(out$mode)
  } else {
    out$Sigma.tilde <- solve(-hessian.integrand(out$mode))
  }

  if(gradient) {
    out$gradient <- grad.integrand(out$mode)
  }

  return(out)
}

##' Laplace Sampling Markov Chain Monte Carlo (MCMC) for Generalized Linear Gaussian Process Models
##'
##' Performs MCMC sampling using Laplace approximation for Generalized Linear Gaussian Process Models (GLGPMs).
##'
##' @param y Response variable vector.
##' @param units_m Units of measurement for the response variable.
##' @param mu Mean vector of the response variable.
##' @param Sigma Covariance matrix of the spatial process.
##' @param ID_coords Indices mapping response to locations.
##' @param ID_re Indices mapping response to unstructured random effects.
##' @param sigma2_re Variance of the unstructured random effects.
##' @param family Distribution family for the response variable. Must be one of 'gaussian', 'binomial', or 'poisson'.
##' @param control_mcmc List with control parameters for the MCMC algorithm:
##'   \describe{
##'     \item{n_sim}{Number of MCMC iterations.}
##'     \item{burnin}{Number of burn-in iterations.}
##'     \item{thin}{Thinning parameter for saving samples.}
##'     \item{h}{Step size for proposal distribution. Defaults to 1.65/(n_tot^(1/6)).}
##'     \item{c1.h, c2.h}{Parameters for adaptive step size tuning.}
##'   }
##' @param Sigma_pd Precision matrix (optional) for Laplace approximation.
##' @param mean_pd Mean vector (optional) for Laplace approximation.
##' @param messages Logical; if TRUE, print progress messages.
##'
##' @details
##' This function implements a Laplace sampling MCMC approach for GLGPMs. It maximizes the integrand using `maxim.integrand` function for Laplace approximation if `Sigma_pd` and `mean_pd` are not provided.
##'
##' The MCMC procedure involves adaptive step size adjustment based on the acceptance probability (`acc_prob`) and uses a Gaussian proposal distribution centered on the current mean (`mean_curr`) with variance `h`.
##'
##' @return An object of class "mcmc.RiskMap" containing:
##'   \describe{
##'     \item{samples$S}{Samples of the spatial process.}
##'     \item{samples$<re_names[i]>}{Samples of each unstructured random effect, named according to columns of ID_re if provided.}
##'     \item{tuning_par}{Vector of step size (h) values used during MCMC iterations.}
##'     \item{acceptance_prob}{Vector of acceptance probabilities across MCMC iterations.}
##'   }
##' @export
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
Laplace_sampling_MCMC <- function(y, units_m, mu, Sigma,
                                  ID_coords, ID_re = NULL,
                                  sigma2_re = NULL,
                                  family, control_mcmc,
                                  Sigma_pd=NULL, mean_pd=NULL, messages = TRUE) {

  Sigma.inv <- solve(Sigma)
  n_loc <- nrow(Sigma)
  n <- length(y)
  if((!is.null(ID_re) & is.null(sigma2_re)) | (is.null(ID_re) & !is.null(sigma2_re))) {
    stop("To introduce unstructured random effects both `ID_re` and `sigma2_re`
           must be provided.")
  }
  n_re <- length(sigma2_re)
  if(n_re > 0) {
    n_dim_re <- sapply(1:n_re, function(i) length(unique(ID_re[,i])))
    ind_re <- list()
    add_i <- 0
    for(i in 1:n_re) {
      ind_re[[i]] <- (add_i+n_loc+1):(add_i+n_loc+n_dim_re[i])
      if(i < n_re) add_i <- sum(n_dim_re[1:i])
    }
  }
  n_tot <- n_loc
  if(n_re > 0) n_tot <- n_tot + sum(n_dim_re)

  C_S <- t(sapply(1:n_loc,function(i) ID_coords==i))

  if(n_re>0) {
    C_re <- list()
    C_S_re <- list()
    C_re_re <- list()
    for(j in 1:n_re){
      C_S_re[[j]] <- array(FALSE,dim = c(n_loc, n_dim_re[j], n))
      C_re[[j]] <- t(sapply(1:n_dim_re[j],function(i) ID_re[,j]==i))
      for(l in 1:n_dim_re[j]) {
        for(k in 1:n_loc) {
          ind_kl <- which(ID_coords==k & ID_re[,j]==l)
          if(length(ind_kl) > 0) {
            C_S_re[[j]][k,l,ind_kl] <- TRUE
          }
        }
      }

      if(j < n_re) {
        C_re_re[[j]] <- list()
        counter <- 0
        for(w in (j+1):n_re) {
          counter <- counter+1
          C_re_re[[j]][[counter]] <- array(FALSE,dim = c(n_dim_re[j], n_dim_re[w], n))
          for(l in 1:n_dim_re[j]) {
            for(k in 1:n_dim_re[w]) {
              ind_lk <- which(ID_re[,j]==l & ID_re[,w]==k)
              if(length(ind_kl) > 0) {
                C_re_re[[j]][[counter]][l,k,ind_lk] <- TRUE
              }
            }
          }
        }
      }
    }
  }

  if(is.null(Sigma_pd) | is.null(mean_pd)) {
    out_maxim <-
      maxim.integrand(y = y, units_m = units_m, Sigma = Sigma, mu = mu,
                      ID_coords = ID_coords, ID_re = ID_re,
                      sigma2_re = sigma2_re,
                      family = family,
                      hessian = FALSE, gradient = TRUE)

    if(is.null(Sigma_pd)) Sigma_pd <- out_maxim$Sigma.tilde
    if(is.null(mean_pd)) mean_pd <- out_maxim$mode
  }

  n_sim <- control_mcmc$n_sim
  n <- length(y)
  Sigma_pd_sroot <- t(chol(Sigma_pd))
  A <- solve(Sigma_pd_sroot)

  if(n_re == 0) {
    Sigma_tot <- Sigma
  } else {
    Sigma_tot <- matrix(0, n_tot, n_tot)
    Sigma_tot[1:n_loc, 1:n_loc] <- Sigma
    for(i in 1:n_re) {
      diag(Sigma_tot)[ind_re[[i]]] <- sigma2_re[i]
    }
  }
  Sigma_w_inv <- solve(A%*%Sigma_tot%*%t(A))
  mu_w <- -as.numeric(A%*%mean_pd)

  cond.dens.W <- function(W, S_tot) {
    S <- S_tot[1:n_loc]
    if(n_re > 0) {
      S_re <- NULL
      S_re_list <- list()
      for(i in 1:n_re) {
        S_re_list[[i]] <- S_tot[ind_re[[i]]]
      }
    }

    eta <- mu + S[ID_coords]
    if(n_re > 0) {
      for(i in 1:n_re) {
        eta <- eta + S_re_list[[i]][ID_re[,i]]
      }
    }

    if(family=="poisson") {
      llik <- sum(y*eta-units_m*exp(eta))
    } else if(family=="binomial") {
      llik <- sum(y*eta-units_m*log(1+exp(eta)))
    }
    diff_w <- W-mu_w
    -0.5*as.numeric(t(diff_w)%*%Sigma_w_inv%*%diff_w)+
    llik
  }

  lang.grad <- function(W, S_tot) {
    diff.w <- W-mu_w
    S <- S_tot[1:n_loc]
    if(n_re > 0) {
      S_re <- NULL
      S_re_list <- list()
      for(i in 1:n_re) {
        S_re_list[[i]] <- S_tot[ind_re[[i]]]
      }
    }

    eta <- mu + S[ID_coords]
    if(n_re > 0) {
      for(i in 1:n_re) {
        eta <- eta + S_re_list[[i]][ID_re[,i]]
      }
    }

    if(family=="poisson") {
      der <- units_m*exp(eta)
    } else if(family=="binomial") {
      der <- units_m*exp(eta)/(1+exp(eta))
    }

    grad_S_tot_r <- rep(NA,n_tot)
    grad_S_tot_r[1:n_loc] <- as.numeric(sapply(1:n_loc,function(i) sum((y-der)[C_S[i,]])))
    if(n_re>0) {
      for(j in 1:n_re) {
        grad_S_tot_r[ind_re[[j]]] <- as.numeric(sapply(1:n_dim_re[[j]],
                                                function(x) sum((y-der)[C_re[[j]][x,]])))
      }
    }

    out <- as.numeric(-Sigma_w_inv%*%(W-mu_w)+
                   t(Sigma_pd_sroot)%*%grad_S_tot_r)
  }

  h <- control_mcmc$h
  if(is.null(h)) h <- 1.65/(n_tot^(1/6))
  burnin <- control_mcmc$burnin
  thin <- control_mcmc$thin
  c1.h <- control_mcmc$c1.h
  c2.h <- control_mcmc$c2.h
  W_curr <- rep(0,n_tot)
  S_tot_curr <- as.numeric(Sigma_pd_sroot%*%W_curr+mean_pd)
  mean_curr <- as.numeric(W_curr + (h^2/2)*lang.grad(W_curr, S_tot_curr))
  lp_curr <- cond.dens.W(W_curr, S_tot_curr)
  acc <- 0
  n_samples <- (n_sim-burnin)/thin
  sim <- matrix(NA,nrow=n_samples, ncol=n_tot)

  if(messages) message("\n - Conditional simulation (burnin=",
                     control_mcmc$burnin,", thin=",control_mcmc$thin,"): \n \n",sep="")
  h.vec <- rep(NA,n_sim)
  acc_prob <- rep(NA,n_sim)
  for(i in 1:n_sim) {
    W_prop <- mean_curr+h*rnorm(n_tot)
    S_tot_prop <-  as.numeric(Sigma_pd_sroot%*%W_prop+mean_pd)
    mean_prop <- as.numeric(W_prop + (h^2/2)*lang.grad(W_prop, S_tot_prop))
    lp_prop <- cond.dens.W(W_prop, S_tot_prop)

    dprop_curr <- -sum((W_prop-mean_curr)^2)/(2*(h^2))
    dprop_prop <- -sum((W_curr-mean_prop)^2)/(2*(h^2))

    log_prob <- lp_prop+dprop_prop-lp_curr-dprop_curr

    if(log(runif(1)) < log_prob) {
      acc <- acc+1
      W_curr <- W_prop
      S_tot_curr <- S_tot_prop
      lp_curr <- lp_prop
      mean_curr <- mean_prop
    }

    if( i > burnin & (i-burnin)%%thin==0) {
      cnt <- (i-burnin)/thin
      sim[cnt,] <- S_tot_curr
    }

    acc_prob[i] <- acc/i
    h.vec[i] <- h <- max(10e-20,h + c1.h*i^(-c2.h)*(acc/i-0.57))

  }

  out_sim <- list()
  out_sim$samples <- list()
  out_sim$samples$S <- sim[,1:n_loc]
  if(n_re > 0) {
    re_names <- colnames(ID_re)
    for(i in 1:n_re) {
      out_sim$samples[[re_names[i]]] <- sim[,ind_re[[i]]]
    }
  }
  out_sim$tuning_par <- h.vec
  out_sim$acceptance_prob <- acc_prob
  class(out_sim) <- "mcmc.RiskMap"
  return(out_sim)
}

##' Set Control Parameters for Simulation
##'
##' This function sets control parameters for running simulations, particularly for MCMC methods.
##' It allows users to specify the number of simulations, burn-in period, thinning interval, and various
##' other parameters necessary for the simulation.
##'
##' @param n_sim Integer. The total number of simulations to run. Default is 12000.
##' @param burnin Integer. The number of initial simulations to discard (burn-in period, used for the MCMC algorithm). Default is 2000.
##' @param thin Integer. The interval at which simulations are recorded (thinning interval, used for the MCMC algorithm). Default is 10.
##' @param h Numeric. An optional parameter. Must be non-negative if specified.
##' @param c1.h Numeric. A control parameter for the simulation. Must be positive. Default is 0.01.
##' @param c2.h Numeric. Another control parameter for the simulation. Must be between 0 and 1. Default is 1e-04.
##' @param linear_model Logical. If TRUE, the function sets up parameters for a linear model and
##' only returns \code{n_sim}. Default is FALSE.
##'
##' @details
##' The function validates the input parameters and ensures they are appropriate for the simulation that is used
##' in the \code{\link{glgpm}} fitting function.
##' For non-linear models, it checks that \code{n_sim} is greater than \code{burnin}, that \code{thin} is positive
##' and a divisor of \code{(n_sim - burnin)}, and that \code{h}, \code{c1.h}, and \code{c2.h} are within their
##' respective valid ranges.
##'
##' If \code{linear_model} is TRUE, only \code{n_sim} and \code{linear_model} are required, and the function
##' returns a list containing these parameters.
##'
##' If \code{linear_model} is FALSE, the function returns a list containing \code{n_sim}, \code{burnin}, \code{thin},
##' \code{h}, \code{c1.h}, \code{c2.h}, and \code{linear_model}.
##'
##' @return A list of control parameters for the simulation with class attribute "mcmc.RiskMap".
##'
##' @examples
##' # Example with default parameters
##' control_params <- set_control_sim()
##'
##' # Example with custom parameters
##' control_params <- set_control_sim(n_sim = 15000, burnin = 3000, thin = 20)
##'
##' @seealso \code{\link[Matrix]{Matrix}}, \code{\link[Matrix]{forceSymmetric}}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##' @importFrom Matrix Matrix forceSymmetric
##' @export
set_control_sim <- function (n_sim  = 12000, burnin = 2000, thin = 10, h = NULL, c1.h = 0.01, c2.h = 1e-04,
                             linear_model = FALSE)
{

  if (!linear_model & n_sim < burnin)
    stop("n_sim cannot be smaller than burnin.")
  if (!linear_model & thin <= 0)
    stop("thin must be positive")
  if (!linear_model & (n_sim - burnin)%%thin != 0)
    stop("thin must be a divisor of (n_sim-burnin)")
  if (!linear_model & !is.null(h) && h < 0)
    stop("h must be positive.")
  if (!linear_model & c1.h < 0)
    stop("c1.h must be positive.")
  if (!linear_model & (c2.h < 0 | c2.h > 1))
    stop("c2.h must be between 0 and 1.")
  if(linear_model) {
    res <- list(n_sim = n_sim, linear_model = linear_model)
  } else {
    res <- list(n_sim = n_sim, burnin = burnin, thin = thin,
                h = h, c1.h = c1.h, c2.h = c2.h, linear_model = linear_model)
  }
  class(res) <- "mcmc.RiskMap"
  return(res)
}

##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @importFrom Matrix Matrix forceSymmetric
glgpm_nong <-
function(y, D, coords, units_m, kappa,
           par0, cov_offset,
           ID_coords, ID_re, s_unique, re_unique,
           fix_tau2, family,return_samples,
           start_beta,
           start_cov_pars,
           control_mcmc,
           messages = TRUE) {

  beta0 <- par0$beta
  mu0 <- D%*%beta0+cov_offset

  sigma2_0 <- par0$sigma2

  phi0 <- par0$phi

  tau2_0 <- par0$tau2

  if(is.null(tau2_0)) tau2_0 <- fix_tau2

  sigma2_re_0 <- par0$sigma2_re

  n_loc <- nrow(coords)
  n_re <- length(sigma2_re_0)
  n_samples <- (control_mcmc$n_sim-control_mcmc$burnin)/control_mcmc$thin

  u = dist(coords)

  Sigma0 <- sigma2_0*matern_cor(u = u, phi = phi0, kappa = kappa,
                                return_sym_matrix = TRUE)

  diag(Sigma0) <- diag(Sigma0) + tau2_0

  sigma2_re_0 <- par0$sigma2_re


  if(messages) message("\n - Obtaining covariance matrix and mean for the proposal distribution of the MCMC \n \n")
  out_maxim <-
    maxim.integrand(y = y, units_m = units_m, Sigma = Sigma0, mu = mu0,
                    ID_coords = ID_coords, ID_re = ID_re,
                    sigma2_re = sigma2_re_0,
                    family = family,
                    hessian = FALSE, gradient = TRUE)

  Sigma_pd <- out_maxim$Sigma.tilde
  mean_pd <- out_maxim$mode

  simulation <-
    Laplace_sampling_MCMC(y = y, units_m = units_m, mu = mu0, Sigma = Sigma0,
                          sigma2_re = sigma2_re_0,
                          ID_coords = ID_coords, ID_re = ID_re,
                          family = family, control_mcmc = control_mcmc,
                          Sigma_pd = Sigma_pd, mean_pd = mean_pd,
                          messages = messages)

  S_tot_samples <- simulation$samples$S

  p <- ncol(D)

  ind_beta <- 1:p

  ind_sigma2 <- p+1

  ind_phi <- p+2

  if(!is.null(fix_tau2)) {
    if(n_re>0) {
      ind_sigma2_re <- (p+2+1):(p+2+n_re)
      n_dim_re <- sapply(1:n_re, function(i) length(unique(ID_re[,i])))
    }
  } else {
    ind_nu2 <- p+3
    if(n_re>0) {
      ind_sigma2_re <- (p+3+1):(p+3+n_re)
      n_dim_re <- sapply(1:n_re, function(i) length(unique(ID_re[,i])))
    }
  }

  if(n_re> 0) {
    for(i in 1:n_re) {
      S_tot_samples <- cbind(S_tot_samples, simulation$samples[[i+1]])
    }
    ind_re <- list()
    add_i <- 0
    for(i in 1:n_re) {
      ind_re[[i]] <- (add_i+n_loc+1):(add_i+n_loc+n_dim_re[i])
      if(i < n_re) add_i <- sum(n_dim_re[1:i])
    }
  }


  log.integrand <- function(S_tot, val) {
    n <- length(y)
    S <- S_tot[1:n_loc]

    q.f_re <- 0
    if(n_re > 0) {
      S_re <- NULL
      S_re_list <- list()
      for(i in 1:n_re) {
        S_re_list[[i]] <- S_tot[ind_re[[i]]]
        q.f_re <- q.f_re + n_dim_re[i]*log(val$sigma2_re[i])+
          sum(S_re_list[[i]]^2)/val$sigma2_re[i]
      }
    } else {
      q.f_re <- 0
    }

    eta <- val$mu + S[ID_coords]
    if(n_re > 0) {
      for(i in 1:n_re) {
        eta <- eta + S_re_list[[i]][ID_re[,i]]
      }
    }
    if(family=="poisson") {
      llik <-  sum(y*eta-units_m*exp(eta))
    } else if(family=="binomial") {
      llik <- sum(y*eta-units_m*log(1+exp(eta)))
    }
    q.f_S <- n_loc*log(val$sigma2)+val$ldetR+t(S)%*%val$R.inv%*%S/val$sigma2
    out <- -0.5*(q.f_S+q.f_re)+llik
    return(out)
  }

  compute.log.f <- function(par,ldetR=NA,R.inv=NA) {
    beta <- par[ind_beta]
    sigma2 <- exp(par[ind_sigma2])
    if(length(fix_tau2)>0) {
      nu2 <- fix_tau2/sigma2
    } else {
      nu2 <- exp(par[ind_nu2])
    }
    phi <- exp(par[ind_phi])
    val <- list()
    val$sigma2 <- sigma2
    val$mu <- as.numeric(D%*%beta)+cov_offset
    if(n_re > 0) {
      val$sigma2_re <- exp(par[ind_sigma2_re])
    }
    if(is.na(ldetR) & is.na(as.numeric(R.inv)[1])) {
      R <- matern_cor(u, phi = phi, kappa=kappa,return_sym_matrix = TRUE)
      diag(R) <- diag(R)+nu2
      val$ldetR <- determinant(R)$modulus
      val$R.inv <- solve(R)
    } else {
      val$ldetR <- ldetR
      val$R.inv <- R.inv
    }
    sapply(1:n_samples,
           function(i) log.integrand(S_tot_samples[i,],val))
  }

  par0_vec <- c(par0$beta, log(c(par0$sigma2, par0$phi)))

  if(is.null(fix_tau2)) {
    par0_vec <- c(par0_vec, log(par0$tau2/par0$sigma2))
  }

  if(n_re > 0) {
    par0_vec <- c(par0_vec, log(par0$sigma2_re))
  }

  log.f.tilde <- compute.log.f(par0_vec)

  MC.log.lik <- function(par) {
    log(mean(exp(compute.log.f(par)-log.f.tilde)))
  }

  grad.MC.log.lik <- function(par) {
    beta <- par[ind_beta]; mu <- as.numeric(D%*%beta)+cov_offset
    sigma2 <- exp(par[ind_sigma2])
    if(length(fix_tau2)>0) {
      nu2 <- fix_tau2/sigma2
    } else {
      nu2 <- exp(par[ind_nu2])
    }
    phi <- exp(par[ind_phi])
    if(n_re > 0) {
      sigma2_re <- exp(par[ind_sigma2_re])
    }

    R <- matern_cor(u, phi = phi, kappa=kappa,return_sym_matrix = TRUE)
    diag(R) <- diag(R)+nu2

    R.inv <- solve(R)
    ldetR <- determinant(R)$modulus

    exp.fact <- exp(compute.log.f(par,ldetR,R.inv)-log.f.tilde)
    L.m <- sum(exp.fact)
    exp.fact <- exp.fact/L.m

    R1.phi <- matern.grad.phi(u,phi,kappa)
    m1.phi <- R.inv%*%R1.phi
    t1.phi <- -0.5*sum(diag(m1.phi))
    m2.phi <- m1.phi%*%R.inv; rm(m1.phi)

    if(is.null(fix_tau2)){
      t1.nu2 <- -0.5*sum(diag(R.inv))
      m2.nu2 <- R.inv%*%R.inv
    }

    gradient.S <- function(S_tot) {
      S <- S_tot[1:n_loc]

      if(n_re > 0) {
        S_re_list <- list()
        for(i in 1:n_re) {
          S_re_list[[i]] <- S_tot[ind_re[[i]]]
        }
      }

      eta <- mu + S[ID_coords]
      if(n_re > 0) {
        for(i in 1:n_re) {
          eta <- eta + S_re_list[[i]][ID_re[,i]]
        }
      }


      if(family=="poisson") {
        h <- units_m*exp(eta)
      } else if(family=="binomial") {
        h <- units_m*exp(eta)/(1+exp(eta))
      }

      q.f_S <- t(S)%*%R.inv%*%S

      grad.beta <-  t(D)%*%(y-h)

      grad.log.sigma2 <- (-n_loc/(2*sigma2)+0.5*q.f_S/(sigma2^2))*sigma2

      grad.log.phi <- (t1.phi+0.5*as.numeric(t(S)%*%m2.phi%*%(S))/sigma2)*phi

      out <- c(grad.beta,grad.log.sigma2,grad.log.phi)

      if(is.null(fix_tau2)) {
        grad.log.nu2 <-  (t1.nu2+0.5*as.numeric(t(S)%*%m2.nu2%*%(S))/sigma2)*nu2
        out <- c(out,grad.log.nu2)
      }

      if(n_re > 0) {
        grad.log.sigma2_re <- rep(NA, n_re)
        for(i in 1:n_re) {
          grad.log.sigma2_re[i] <- (-n_dim_re[i]/(2*sigma2_re[i])+0.5*sum(S_re_list[[i]]^2)/
                                      (sigma2_re[i]^2))*sigma2_re[i]
        }
        out <- c(out,grad.log.sigma2_re)
      }
      out
    }
    out <- rep(0,length(par))
    for(i in 1:n_samples) {
      out <- out + exp.fact[i]*gradient.S(S_tot_samples[i,])
    }
    out
  }

  hess.MC.log.lik <- function(par) {
    beta <- par[ind_beta]; mu <- as.numeric(D%*%beta)+cov_offset
    sigma2 <- exp(par[ind_sigma2])
    if(!is.null(fix_tau2)) {
      nu2 <- fix_tau2/sigma2
    } else {
      nu2 <- exp(par[ind_nu2])
    }
    phi <- exp(par[ind_phi])
    if(n_re > 0) {
      sigma2_re <- exp(par[ind_sigma2_re])
    }

    R <- matern_cor(u, phi = phi, kappa=kappa,return_sym_matrix = TRUE)
    diag(R) <- diag(R)+nu2

    R.inv <- solve(R)
    ldetR <- determinant(R)$modulus

    exp.fact <- exp(compute.log.f(par,ldetR,R.inv)-log.f.tilde)
    L.m <- sum(exp.fact)
    exp.fact <- exp.fact/L.m

    R1.phi <- matern.grad.phi(u,phi,kappa)
    m1.phi <- R.inv%*%R1.phi
    t1.phi <- -0.5*sum(diag(m1.phi))
    m2.phi <- m1.phi%*%R.inv; rm(m1.phi)

    if(is.null(fix_tau2)){
      t1.nu2 <- -0.5*sum(diag(R.inv))
      m2.nu2 <- R.inv%*%R.inv
      t2.nu2 <- 0.5*sum(diag(m2.nu2))
      n2.nu2 <- 2*R.inv%*%m2.nu2
      t2.nu2.phi <- 0.5*sum(diag(R.inv%*%R1.phi%*%R.inv))
      n2.nu2.phi <- R.inv%*%(R.inv%*%R1.phi+
                               R1.phi%*%R.inv)%*%R.inv
    }

    R2.phi <- matern.hessian.phi(u,phi,kappa)
    t2.phi <- -0.5*sum(diag(R.inv%*%R2.phi-R.inv%*%R1.phi%*%R.inv%*%R1.phi))
    n2.phi <- R.inv%*%(2*R1.phi%*%R.inv%*%R1.phi-R2.phi)%*%R.inv

    H <- matrix(0,nrow=length(par),ncol=length(par))

    hessian.S <- function(S_tot,ef) {
      S <- S_tot[1:n_loc]

      if(n_re > 0) {
        S_re_list <- list()
        for(i in 1:n_re) {
          S_re_list[[i]] <- S_tot[ind_re[[i]]]
        }
      }

      eta <- mu + S[ID_coords]
      if(n_re > 0) {
        for(i in 1:n_re) {
          eta <- eta + S_re_list[[i]][ID_re[,i]]
        }
      }

      if(family=="poisson") {
        h <- units_m*exp(eta)
        h1 <- h
      } else if(family=="binomial") {
        h <- units_m*exp(eta)/(1+exp(eta))
        h1 <- h/(1+exp(eta))
      }

      q.f_S <- t(S)%*%R.inv%*%S

      grad.beta <-  t(D)%*%(y-h)

      grad.log.sigma2 <- (-n_loc/(2*sigma2)+0.5*q.f_S/(sigma2^2))*sigma2

      grad.log.phi <- (t1.phi+0.5*as.numeric(t(S)%*%m2.phi%*%(S))/sigma2)*phi

      g <- c(grad.beta,grad.log.sigma2,grad.log.phi)
      if(is.null(fix_tau2)) {
        grad.log.nu2 <-  (t1.nu2+0.5*as.numeric(t(S)%*%m2.nu2%*%(S))/sigma2)*nu2
        g <- c(g,grad.log.nu2)
      }

      if(n_re > 0) {
        grad.log.sigma2_re <- rep(NA, n_re)
        for(i in 1:n_re) {
          grad.log.sigma2_re[i] <- (-n_dim_re[i]/(2*sigma2_re[i])+0.5*sum(S_re_list[[i]]^2)/
                                      (sigma2_re[i]^2))*sigma2_re[i]
        }
        g <- c(g,grad.log.sigma2_re)
      }

      grad2.log.lsigma2.lsigma2 <- (n_loc/(2*sigma2^2)-q.f_S/(sigma2^3))*sigma2^2+
        grad.log.sigma2

      grad2.log.lphi.lphi <-(t2.phi-0.5*t(S)%*%n2.phi%*%(S)/sigma2)*phi^2+
        grad.log.phi

      H[ind_beta, ind_beta] <- -t(D)%*%(D*h1)
      H[ind_sigma2, ind_sigma2] <-  grad2.log.lsigma2.lsigma2
      H[ind_sigma2,ind_phi] <-
        H[ind_phi, ind_sigma2] <- (grad.log.phi/phi-t1.phi)*(-phi)
      H[ind_phi,ind_phi] <- grad2.log.lphi.lphi

      if(is.null(fix_tau2)) {
        grad2.log.lnu2.lnu2 <- (t2.nu2-0.5*t(S)%*%n2.nu2%*%(S)/sigma2)*nu2^2+
          grad.log.nu2
        grad2.log.lnu2.lphi <- (t2.nu2.phi-0.5*t(S)%*%n2.nu2.phi%*%(S)/sigma2)*phi*nu2
        H[ind_sigma2,ind_nu2] <- H[ind_nu2,ind_sigma2] <- (grad.log.nu2/nu2-t1.nu2)*(-nu2)
        H[ind_nu2,ind_nu2] <- grad2.log.lnu2.lnu2
        H[ind_phi,ind_nu2] <- H[ind_nu2,ind_phi] <- grad2.log.lnu2.lphi
      }

      if(n_re > 0) {
        grad2.log.sigma2_re <- rep(NA, n_re)
        for(i in 1:n_re) {
          grad2.log.sigma2_re[i] <- (n_dim_re[i]/(2*sigma2_re[i]^2)-
                                       sum(S_re_list[[i]]^2)/(sigma2_re[i]^3))*
            sigma2_re[i]^2+grad.log.sigma2_re[i]
          H[ind_sigma2_re[i],ind_sigma2_re[i]] <- grad2.log.sigma2_re[i]

        }
      }
      out <- list()
      out$mat1<- ef*(g%*%t(g)+H)
      out$g <- g*ef
      out
    }

    a <- rep(0,length(par))
    A <- matrix(0,length(par),length(par))
    for(i in 1:n_samples) {
      out.i <- hessian.S(S_tot_samples[i,],exp.fact[i])
      a <- a+out.i$g
      A <- A+out.i$mat1
    }
    (A-a%*%t(a))

  }

  start_cov_pars[-(1:2)] <- start_cov_pars[-(1:2)]/start_cov_pars[1]
  start_par <- c(start_beta, log(start_cov_pars))

  out <- list()
  estim <- nlminb(start_par,
                  function(x) -MC.log.lik(x),
                  function(x) -grad.MC.log.lik(x),
                  function(x) -hess.MC.log.lik(x),
                  control=list(trace=1*messages))

  out$estimate <- estim$par
  out$grad.MLE <- grad.MC.log.lik(estim$par)
  hess.MLE <- hess.MC.log.lik(estim$par)
  out$covariance <- solve(-hess.MLE)
  out$log.lik <- -estim$objective
  if(return_samples) out$S_samples <- S_tot_samples
  class(out) <- "RiskMap"
  return(out)
}

##' Check MCMC Convergence for Spatial Random Effects
##'
##' This function checks the Markov Chain Monte Carlo (MCMC) convergence of spatial random effects
##' for either a \code{RiskMap} or \code{RiskMap.pred.re} object.
##' It plots the trace plot and autocorrelation function (ACF) for the MCMC chain
##' and calculates the effective sample size (ESS).
##'
##' @param object An object of class \code{RiskMap} or \code{RiskMap.pred.re}.
##'  \code{RiskMap} is the output from \code{\link{glgpm}} function, and
##'  \code{RiskMap.pred.re} is obtained from the \code{\link{pred_over_grid}} function.
##' @param check_mean Logical. If \code{TRUE}, checks the MCMC chain for the mean of the spatial random effects.
##'  If \code{FALSE}, checks the chain for a specific component of the random effects vector.
##' @param component Integer. The index of the spatial random effects component to check when \code{check_mean = FALSE}.
##'  Must be a positive integer corresponding to a location in the data. Ignored if \code{check_mean = TRUE}.
##' @param ... Additional arguments passed to the \code{\link[stats]{acf}} function for customizing the ACF plot.
##'
##' @details
##' The function first checks that the input object is either of class \code{RiskMap} or \code{RiskMap.pred.re}.
##' Depending on the value of \code{check_mean}, it either calculates the mean of the spatial random effects
##' across all locations for each iteration or uses the specified component.
##' It then generates two plots:
##' - A trace plot of the selected spatial random effect over iterations.
##' - An autocorrelation plot (ACF) with the effective sample size (ESS) displayed in the title.
##'
##' The ESS is computed using the \code{\link[sns]{ess}} function, which provides a measure of the effective number
##' of independent samples in the MCMC chain.
##'
##' If \code{check_mean = TRUE}, the \code{component} argument is ignored, and a warning is issued.
##' To specify a particular component of the random effects vector, set \code{check_mean = FALSE} and provide
##' a valid \code{component} value.
##'
##' @return
##' No return value, called for side effects (plots and warnings).
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @importFrom sns ess
##' @importFrom graphics par
##' @importFrom stats acf
##' @export
check_mcmc <- function(object, check_mean = TRUE,
                       component = NULL, ...) {
  if(!inherits(object,
               what = "RiskMap", which = FALSE) &
     !inherits(object,
               what = "RiskMap.pred.re", which = FALSE)) {
    stop("'object' must be either of one of these objects:
           a 'RiskMap' object obtained as an output from glgpm;
           a 'RiskMap.pred.re' object obtained as an output from 'pred_over_grid'")
  }

  if(inherits(object,
              what = "RiskMap", which = FALSE)) {
    S_samples <- object$S_samples
  } else if (inherits(object,
                      what = "RiskMap.pred.re", which = FALSE)) {
    S_samples <- t(object$S_samples)
  }

  if(check_mean & !is.null(component)) {
    warning("if check_mean = TRUE, the value passed to 'component' is ignored;
            set check_mean = FALSE when specifying a value for 'component'")
  }
  n_samples <- nrow(S_samples)
  n_loc <- ncol(S_samples)
  if(check_mean) {
    S_chain <- apply(S_samples, 1, mean)
  } else {
    if(is.null(component)) stop("When check_mean = FALSE a component of the
                                random effects vector must be specified through 'component'
                                by providing a positive integer")
    if(component < 0 | component > n_loc) stop("'component' must be a positive integer
                                              between 1 and the number of locations in the data")
    S_chain <- S_samples[,component]
  }

  par(mfrow = c(1,2))
  plot(S_chain, type = "l",
       ylab = "", xlab = "Iteration",
       main = "Spatial random effect")

  S_chain_ess <- round(ess(S_chain),3)
  acf(S_chain, main = paste("Effective sample size:",
                            S_chain_ess), ...)
  par(mfrow = c(1,1))
}
