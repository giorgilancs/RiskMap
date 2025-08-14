
compute_mda_effect <- function(survey_times_data, mda_times, intervention,
                               alpha, gamma, kappa) {
  n <- length(survey_times_data)
  effect <- rep(NA, n)

  mda_effect_f <- function(v, alpha, gamma, kappa) {
    alpha*exp(-(v/gamma)^kappa)
  }
  mda_effect_f <- Vectorize(mda_effect_f, "v")

  f <- function(t, mda_times, int, alpha, gamma, kappa) {
    ind_t <- which(t > mda_times)
    u_j <- mda_times[ind_t]
    if(length(u_j) > 0) {
      out <- prod((1-mda_effect_f(t-u_j, alpha, gamma, kappa))^int[ind_t])
    } else {
      out <- 1
    }
    return(out)
  }

  f <- Vectorize(f, "t")

  for(i in 1:n) {
    effect[i] <- f(survey_times_data[i], mda_times, intervention[i,],
                   alpha, gamma, kappa)
  }
  return(effect)
}


compute_mda_effect_derivatives <- function(survey_times_data, mda_times, intervention,
                                           alpha, gamma, kappa) {
  n <- length(survey_times_data)
  effect <- numeric(n)
  d_alpha <- numeric(n)
  d_gamma <- numeric(n)
  d2_alpha <- numeric(n)
  d2_gamma <- numeric(n)
  d2_alpha_gamma <- numeric(n)

  for (i in 1:n) {
    t <- survey_times_data[i]
    int <- intervention[i,]

    ind_t <- which(t > mda_times)
    u_j <- mda_times[ind_t]
    int_j <- int[ind_t]

    if (length(u_j) > 0) {
      v_j <- t - u_j
      z_j <- (v_j/gamma)^kappa
      e_j <- exp(-z_j)
      denom <- 1 - alpha * e_j
      effect[i] <- exp(sum(int_j * log(denom)))

      # First derivatives (CORRECT - matches numeric)
      term_alpha <- int_j * e_j / denom
      dlogP_dalpha <- -sum(term_alpha)

      gamma_term <- kappa * z_j / gamma
      term_gamma <- alpha * term_alpha * gamma_term
      dlogP_dgamma <- -sum(term_gamma)

      # Second derivatives
      # ∂²(log P)/∂α² (CORRECT - matches numeric)
      d2logP_dalpha2 <- -sum(int_j * e_j^2 / denom^2)

      #######################################################
      # CORRECTED ∂²(log P)/∂γ² (FIXED SIGN)
      # f = αe^{-z}, z = (v/γ)^κ
      # ∂f/∂γ = f * (κz/γ)
      # ∂²f/∂γ² = f * [κz(κ+1 - κz)/γ²]
      term1 <- (alpha * e_j * kappa * z_j * (kappa + 1 - kappa*z_j)) /
        (gamma^2 * denom)
      term2 <- (alpha * e_j * kappa * z_j / (gamma * denom))^2
      d2logP_dgamma2 <- sum(int_j * (term1 - term2))
      #######################################################

      # ∂²(log P)/∂α∂γ (CORRECT - matches numeric)
      term_cross <- e_j * gamma_term / denom
      d2logP_dalphadgamma <- sum(int_j * (-term_cross - alpha * e_j * term_cross / denom))

      # Convert to derivatives of P
      d_alpha[i] <- effect[i] * dlogP_dalpha
      d_gamma[i] <- effect[i] * dlogP_dgamma
      d2_alpha[i] <- effect[i] * (dlogP_dalpha^2 + d2logP_dalpha2)
      d2_gamma[i] <- effect[i] * (dlogP_dgamma^2 + d2logP_dgamma2)
      d2_alpha_gamma[i] <- effect[i] * (dlogP_dalpha * dlogP_dgamma + d2logP_dalphadgamma)

    } else {
      effect[i] <- 1
      d_alpha[i] <- 0
      d_gamma[i] <- 0
      d2_alpha[i] <- 0
      d2_gamma[i] <- 0
      d2_alpha_gamma[i] <- 0
    }
  }

  return(list(effect = effect,
              d_alpha = d_alpha,
              d_gamma = d_gamma,
              d2_alpha = d2_alpha,
              d2_gamma = d2_gamma,
              d2_alpha_gamma = d2_alpha_gamma))
}

dast_initial_value <- function(y, D, units_m, int_mat, survey_times_data,
                               penalty,
                               mda_times, fix_alpha, power_val) {

  p <- ncol(D)
  n <- nrow(D)

  llik <- function(par) {
    beta <- par[1:p]
    if(is.null(fix_alpha)) {
      alpha <- exp(par[p+1])/(1+exp(par[p+1]))
      gamma <- exp(par[p+2])
    } else {
      alpha <- fix_alpha
      gamma <- exp(par[p+1])
    }

    fact <- compute_mda_effect(survey_times_data, mda_times, int_mat,
                               alpha, gamma, kappa = power_val)
    eta <- as.numeric(D%*%beta)
    prob_star <- 1/(1+exp(-eta))
    prob <- fact*prob_star

    out <- -(sum(y*log(prob/(1-prob)) + units_m*log(1-prob))-penalty[[1]](alpha))
    return(out)
  }

  start <- c(runif(p+1,-1,1), mean(dist(survey_times_data)))
  est <- nlminb(start, llik)

  est$beta <- est$par[1:p]
  if(is.null(fix_alpha)) {
    est$alpha <- exp(est$par[p+1])/(1+exp(est$par[p+1]))
    est$gamma <- exp(est$par[p+2])
  } else {
    est$gamma <- exp(est$par[p+1])
  }

  return(est)
}


##' @title Fitting of decay-adjusted spatio-temporal (DAST) model
##'
##' @description
##' The function fits a decay-adjusted spatio-temporal (DAST) model using Monte Carlo maximum liklihood.
##' The DAST model allows for the incorporation of the temporal decay in diease prevalence due
##' to the impact of mass drug administration (MDA). The function requires the full MDA history as detailed in the arguments below.
##'
##' @param formula A model formula specifying the response variable and predictors.
##' @param data A `data.frame` or `sf` object containing the dataset.
##' @param den The denominator for binomial models.
##' @param survey_times The variable indicating the survey times.
##' @param mda_times A vector specifying the mass drug administration (MDA) times.
##' @param int_mat Intervention matrix specifying the timing and coverage of MDA; the dimension of the matrix
##' must be \code{n * n_mda}, where \code{n} is the number of rows of \code{data} and \code{n_mda} is the length of \code{mda_times}.
##' @param penalty Optional list specifying penalty functions for regularization, used in the estimation of the "drop" parameter \code{alpha}.
##' @param drop Optional value used for fixing the "drop" parameter of the MDA impact function.
##' @param power_val Value expressing the power of the MDA impact function.
##' @param crs Optional coordinate reference system (CRS) for spatial data.
##' @param convert_to_crs CRS to which spatial data should be converted.
##' @param scale_to_km Logical; whether to scale distances to kilometers (default: `TRUE`).
##' @param control_mcmc A list of MCMC control parameters, typically from `set_control_sim()`.
##' @param par0 Optional list of initial parameter values.
##' @param S_samples Number of posterior samples to retain.
##' @param return_samples Logical; whether to return posterior samples (default: `TRUE`).
##' @param messages Logical; whether to print messages (default: `TRUE`).
##' @param start_pars List of starting values for parameters.
##'
##' @return A list containing model estimates, posterior samples, and metadata.
##' @return A list containing:
##'
##' - `y`: Response variable values.
##' - `D`: Covariate matrix.
##' - `coords`: Unique spatial coordinates.
##' - `mda_times`: MDA time points.
##' - `survey_times_data`: Survey time data.
##' - `int_mat`: Intervention matrix.
##' - `ID_coords`: Indices of spatial locations.
##' - `re`: Random effects levels (if applicable).
##' - `ID_re`: Indices of random effects (if applicable).
##' - `power_val`: Power of the MDA impact function.
##' - `fix_tau2`: Fixed tau-squared value (if applicable).
##' - `fix_alpha`: Fixed alpha value (if applicable).
##' - `formula`: Model formula.
##' - `crs`: Coordinate reference system.
##' - `scale_to_km`: Indicator of distance scaling.
##' - `data_sf`: Processed spatial dataset.
##' - `family`: Model family (e.g., "binomial").
##' - `kappa`: Smoothness parameter.
##' - `units_m`: Denominator for binomial models.
##' - `cov_offset`: Offset for covariates.
##' - `call`: Function call.
##' - `penalty`: Penalty function details (if applicable).
##' - `posterior_samples`: Posterior samples if `return_samples = TRUE`.
##'
##' @seealso \code{\link{set_control_sim}}, \code{\link{summary.RiskMap}}, \code{\link{to_table}}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##' @export
dast <- function(formula,
                 data,
                 den = NULL, mda_times, int_mat,
                 penalty = NULL,
                 sst = FALSE,
                 drop = NULL, power_val,
                 crs = NULL, convert_to_crs = NULL,
                 scale_to_km = TRUE,
                 control_mcmc = set_control_sim(),
                 par0=NULL,
                 S_samples = NULL,
                 return_samples = TRUE,
                 messages = TRUE,
                 start_pars = list(beta = NULL,
                                   sigma2 = NULL,
                                   tau2 = NULL,
                                   phi = NULL,
                                   psi = NULL,
                                   sigma2_re = NULL,
                                   gamma = NULL,
                                   alpha = NULL)) {

  nong <- TRUE


  if(!inherits(formula,
               what = "formula", which = FALSE)) {
    stop("'formula' must be a 'formula'
         object indicating the variables of the
         model to be fitted")
  }

  inter_f <- interpret.formula(formula)

  if(inter_f$gp.spec$dim != 3) {
    stop("Three elements must be specified through 'gp' in the formula:
         - the x-coordinates;
         - the y-coordinates;
         - the times of obeservation")
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
    if(length(inter_f$gp.spec$term)>1) {
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

  if(is.null(penalty)) {
    no_penalty <- TRUE
    penalty <- list(pn = function(x) return(0),
                    pn_d1 = function(x) return(0),
                    pn_d2 = function(x) return(0))
  } else {
    no_penalty <- FALSE
  }

  kappa <- inter_f$gp.spec$kappa
  if(kappa < 0) stop("kappa must be positive.")


  mf <- model.frame(inter_f$pf,data=data, na.action = na.fail)

  # Extract outcome data
  y <- as.numeric(model.response(mf))
  n <- length(y)

  # Extract covariates matrix
  D <- as.matrix(model.matrix(attr(mf,"terms"),data=data))

  if(is.null(inter_f$offset)) {
    cov_offset <- rep(0, nrow(data))
  } else {
    cov_offset <- data[[inter_f$offset]]
  }

  # Define denominators for Binomial and Poisson distributions

  # units_m <- data$m
  if(nong) {
    do_name <- deparse(substitute(den))
    if(do_name=="NULL") {
      units_m <- rep(1, nrow(data))
      warning("'den' is assumed to be 1 for all observations \n")
    } else {
      units_m <- data[[do_name]]
    }
    if(is.integer(units_m)) units_m <- as.numeric(units_m)
    if(!is.numeric(units_m)) stop("the variable passed to `den` must be numeric")
    if(any(y > units_m)) stop("The counts identified by the outcome variable cannot be larger
                              than `den` in the case of a Binomial distribution")
    if(!inherits(control_mcmc,
                 what = "mcmc.RiskMap", which = FALSE)) {
      stop ("the argument passed to 'control_mcmc' must be an output
                                                  from the function set_control_sim; see ?set_control_sim
                                                  for more details")

    }

  }

  survey_times_data <- data[[inter_f$gp.spec$term[3]]]

  if(length(inter_f$re.spec) > 0) {
    hr_re <- inter_f$re.spec$term
    re_names <- inter_f$re.spec$term
  } else {
    hr_re <- NULL
  }
  if(!is.null(drop)) {
    fix_alpha <- drop
  } else {
    fix_alpha <- NULL
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

  if(messages) message("The CRS used is ", as.list(st_crs(data))$input, "\n")

  coords_o <- st_coordinates(data)
  if(sst) {
    coords_time <- unique(cbind(coords_o, survey_times_data))
    coords <- coords_time[,1:2]
    time <- coords_time[,3]
  } else {
    coords <- unique(coords_o)
  }

  m <- nrow(coords_o)
  if(sst) {
    ID_coords <- sapply(1:m, function(i)
      which(coords_o[i,1]==coords[,1] &
            coords_o[i,2]==coords[,2] &
            survey_times_data[i]==time))
  } else {
    ID_coords <- sapply(1:m, function(i)
      which(coords_o[i,1]==coords[,1] &
              coords_o[i,2]==coords[,2]))
  }
  s_unique <- unique(ID_coords)

  fix_tau2 <- inter_f$gp.spec$nugget


  if(scale_to_km) {
    coords_o <- coords_o/1000
    coords <- coords/1000
    if(messages) message("Distances between locations are computed in kilometers ")
  } else {
    if(messages) message("Distances between locations are computed in meters ")
  }

  if(is.null(start_pars) & !is.null(par0)) {
    start_pars <- par0
    if(length(par0$beta)!=ncol(D)) stop("the values passed to `beta` in par0 do not match the
                                        variables specified in the formula")
  }

  if(is.null(start_pars$beta)) {
    aux_data <- data.frame(y=y, units_m = units_m, D[,-1])
    if(length(cov_offset)==1) cov_offset_aux <- rep(cov_offset, n)
    glm_fitted <- glm(cbind(y, units_m - y) ~ ., offset = cov_offset,
                      data = aux_data, family = binomial)
    start_pars$beta <- stats::coef(glm_fitted)
  }

  if(is.null(start_pars$gamma) | (is.null(drop) & is.null(start_pars$alpha)) ) {
    dast_i <- dast_initial_value(y, D, units_m, int_mat = int_mat, survey_times_data,
                                 fix_alpha = fix_alpha,penalty=penalty,
                      mda_times, power_val = power_val)
    start_pars$beta <- dast_i$beta
    if(is.null(drop)) {
      start_pars$alpha <- dast_i$alpha
    }
    start_pars$gamma <- dast_i$gamma
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
  if(sst) {
    if(is.null(start_pars$psi)) {
      start_pars$psi <- quantile(dist(unique(survey_times_data)),0.1)
    } else {
      if(start_pars$psi<0) stop("the starting value for psi must be positive")
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



  if(is.null(par0)) {
    par0 <- start_pars
  }
  res <- dast_fit(y = y, D, coords, time, units_m = units_m,
                  mda_times = mda_times, survey_times_data = survey_times_data,
                  sst = sst,
                  int_mat = int_mat,
                  kappa = inter_f$gp.spec$kappa,
                  penalty = penalty,
                      ID_coords, ID_re, s_unique, re_unique,
                      fix_tau2, fix_alpha,
                      return_samples = return_samples,
                      par0 = par0, cov_offset = cov_offset,
                      power_val = power_val,
                      start_beta = start_pars$beta,
                      start_alpha = start_pars$alpha,
                      start_gamma = start_pars$gamma,
                      start_cov_pars = c(start_pars$sigma2,
                                         start_pars$phi,
                                         start_pars$tau2,
                                         start_pars$sigma2_re),
                      start_psi = start_pars$psi,
                      control_mcmc = control_mcmc,
                      messages = messages)



  res$y <- y
  res$D <- D
  res$coords <- coords
  res$sst <- sst
  if(sst) res$time <- time
  res$mda_times <- mda_times
  res$survey_times_data <- survey_times_data
  res$int_mat <- int_mat
  res$ID_coords <- ID_coords
  if(n_re>0) {
    res$re <- re_unique_f
    res$ID_re <- as.data.frame(ID_re)
    colnames(res$ID_re) <- names_re
  }
  res$power_val <- power_val
  res$fix_tau2 <- fix_tau2
  res$fix_alpha <- fix_alpha
  res$formula <- formula
  if(!is.null(convert_to_crs)) {
    crs <- convert_to_crs
  } else {
    crs <- sf::st_crs(data)$input
  }
  if(no_penalty) {
    res$penalty <- NULL
  } else {
    res$penalty <- penalty
  }
  res$crs <- crs
  res$scale_to_km <- scale_to_km
  res$data_sf <- data
  res$family <- "binomial"
  res$kappa <- kappa
  if(nong) res$units_m <- units_m
  res$cov_offset <- cov_offset
  res$call <- match.call()


  invlink <- NULL
  if (is.null(invlink)) {
    inv <- function(eta) stats::plogis(eta)
    d1  <- function(eta) { p <- inv(eta); p * (1 - p) }
    d2  <- function(eta) { p <- inv(eta); d <- p * (1 - p); d * (1 - 2 * p) }
    invlink <- list(inv = inv, d1 = d1, d2 = d2, name = "canonical")
  }
  res$linkf <- invlink

  return(res)
}


dast_fit <-
  function(y, D, coords, time, units_m, kappa, penalty,
           mda_times, survey_times_data,
           sst,
           int_mat,
           par0, cov_offset, power_val,
           ID_coords, ID_re, s_unique, re_unique,
           fix_tau2, fix_alpha, return_samples,
           start_beta,
           start_alpha, start_gamma,
           start_cov_pars,
           start_psi,
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
    if(sst) v = dist(time)

    if(sst) {
      psi0 <- par0$psi
      Sigma0 <- sigma2_0*matern_cor(u = u, phi = phi0, kappa = kappa,
                                    return_sym_matrix = TRUE) *
                         matern_cor(u = v, phi = psi0, kappa = 0.5,
                                    return_sym_matrix = TRUE)
    } else {
      Sigma0 <- sigma2_0*matern_cor(u = u, phi = phi0, kappa = kappa,
                                    return_sym_matrix = TRUE)
    }

    diag(Sigma0) <- diag(Sigma0) + tau2_0

    sigma2_re_0 <- par0$sigma2_re

    alpha0 <- par0$alpha
    if(is.null(alpha0)) alpha0<- fix_alpha

    gamma0 <- par0$gamma

    mda_effect0 <- compute_mda_effect(survey_times_data, mda_times, int_mat,
                       alpha0, gamma0, kappa = power_val)
    if(messages) message("\n - Obtaining covariance matrix and mean for the proposal distribution of the MCMC \n")
    out_maxim <-
      maxim.integrand.dast(y = y, units_m = units_m, Sigma = Sigma0, mu = mu0,
                      mda_effect = mda_effect0,
                      ID_coords = ID_coords, ID_re = ID_re,
                      sigma2_re = sigma2_re_0,
                      hessian = FALSE, gradient = TRUE)

    Sigma_pd <- out_maxim$Sigma.tilde
    mean_pd <- out_maxim$mode

    simulation <-
      Laplace_sampling_MCMC_dast(y = y, units_m = units_m, mu = mu0,
                            mda_effect = mda_effect0, Sigma = Sigma0,
                            sigma2_re = sigma2_re_0,
                            ID_coords = ID_coords, ID_re = ID_re,
                            control_mcmc = control_mcmc,
                            Sigma_pd = Sigma_pd, mean_pd = mean_pd,
                            messages = messages)

    S_tot_samples <- simulation$samples$S

    # Get the number of columns in D
    p <- ncol(D)

    # Define index for beta, sigma2, and phi
    ind_beta <- 1:p
    ind_sigma2 <- p + 1
    ind_phi <- p + 2

    # Conditional indexing based on fix_tau2 and n_re
    if (!is.null(fix_tau2)) {
      if (n_re > 0) {
        ind_sigma2_re <- (p + 2 + 1):(p + 2 + n_re)
        n_dim_re <- sapply(1:n_re, function(i) length(unique(ID_re[, i])))
      } else {
        n_dim_re <- NULL
        ind_re <- NULL
      }

      if (is.null(fix_alpha)) {
        ind_alpha <- p + n_re + 3
        ind_gamma <- p + n_re + 4
        ind_psi <- p + n_re + 5  # psi when alpha is not fixed
      } else {
        ind_gamma <- p + n_re + 3
        ind_psi <- p + n_re + 4  # psi when alpha is fixed
      }
    } else {
      ind_nu2 <- p + 3
      if (n_re > 0) {
        ind_sigma2_re <- (p + 3 + 1):(p + 3 + n_re)
        n_dim_re <- sapply(1:n_re, function(i) length(unique(ID_re[, i])))
      } else {
        n_dim_re <- NULL
        ind_re <- NULL
      }

      if (is.null(fix_alpha)) {
        ind_alpha <- p + n_re + 4
        ind_gamma <- p + n_re + 5
        ind_psi <- p + n_re + 6  # psi when alpha is not fixed
      } else {
        ind_gamma <- p + n_re + 4
        ind_psi <- p + n_re + 5  # psi when alpha is fixed
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

      prob_star <- 1/(1+exp(-eta))
      prob <- val$mda_effect*prob_star

      llik <- sum(y*log(prob)+(units_m-y)*log(1-prob))

      q.f_S <- n_loc*log(val$sigma2)+val$ldetR+t(S)%*%val$R.inv%*%S/val$sigma2
      out <- -0.5*(q.f_S+q.f_re)+llik - val$pen_alpha
      return(out)
    }

    compute.log.f <- function(par,ldetR=NA,R.inv=NA) {
      beta <- par[ind_beta]
      sigma2 <- exp(par[ind_sigma2])

      if(is.null(fix_alpha)) {
        alpha <- exp(par[ind_alpha])/(1+exp(par[ind_alpha]))
      } else {
        alpha <- fix_alpha
      }

      gamma <- exp(par[ind_gamma])

      if(length(fix_tau2)>0) {
        nu2 <- fix_tau2/sigma2
      } else {
        nu2 <- exp(par[ind_nu2])
      }
      phi <- exp(par[ind_phi])

      if(sst) psi <- exp(par[ind_psi])

      val <- list()
      val$sigma2 <- sigma2
      val$mu <- as.numeric(D%*%beta)+cov_offset
      val$mda_effect <- compute_mda_effect(survey_times_data, mda_times,
                                           intervention = int_mat,
                                           alpha, gamma, kappa = power_val)
      val$pen_alpha <- penalty[[1]](alpha)
      if(n_re > 0) {
        val$sigma2_re <- exp(par[ind_sigma2_re])
      }
      if(is.na(ldetR) & is.na(as.numeric(R.inv)[1])) {
        if(sst) {
          R <- matern_cor(u, phi = phi, kappa=kappa,return_sym_matrix = TRUE)*
               matern_cor(v, phi = psi, kappa=kappa,return_sym_matrix = TRUE)
        } else {
          R <- matern_cor(u, phi = phi, kappa=kappa,return_sym_matrix = TRUE)
        }
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


    if(!is.null(fix_tau2)) {
      par0_vec <- c(par0$beta, log(c(par0$sigma2, par0$phi)))
    } else {
      par0_vec <- c(par0$beta, log(c(par0$sigma2, par0$phi, par0$tau2/par0$sigma2)))
    }

    if(n_re > 0) {
      if(is.null(fix_alpha)) {
        par0_vec <- c(par0_vec,log(par0$sigma2_re),
                      log(par0$alpha/(1-par0$alpha)), log(par0$gamma))
      } else {
        par0_vec <- c(par0_vec, log(par0$sigma2_re),
                      log(par0$gamma))
      }
    } else {
      if(is.null(fix_alpha)) {
        par0_vec <- c(par0_vec,
                      log(par0$alpha/(1-par0$alpha)), log(par0$gamma))
      } else {
        par0_vec <- c(par0_vec,
                      log(par0$gamma))
      }
    }
    if(sst) par0_vec <- c(par0_vec, log(par0$psi))

    log.f.tilde <- compute.log.f(par0_vec)

    MC.log.lik <- function(par) {
      log(mean(exp(compute.log.f(par)-log.f.tilde)))
    }

    grad.MC.log.lik <- function(par) {
      beta <- par[ind_beta]
      mu <- as.numeric(D %*% beta) + cov_offset
      sigma2 <- exp(par[ind_sigma2])
      if (is.null(fix_alpha)) {
        alpha <- exp(par[ind_alpha]) / (1 + exp(par[ind_alpha]))
      } else {
        alpha <- fix_alpha
      }
      gamma <- exp(par[ind_gamma])
      mda_effect_all <- compute_mda_effect_derivatives(survey_times_data, mda_times,
                                                       intervention = int_mat,
                                                       alpha, gamma, kappa = power_val)
      mda_effect <- mda_effect_all$effect
      mda_der_alpha <- mda_effect_all$d_alpha
      mda_der_gamma <- mda_effect_all$d_gamma

      if (length(fix_tau2) > 0) {
        nu2 <- fix_tau2 / sigma2
      } else {
        nu2 <- exp(par[ind_nu2])
      }
      phi <- exp(par[ind_phi])

      if(sst) psi <- exp(par[ind_psi])
      if (n_re > 0) {
        sigma2_re <- exp(par[ind_sigma2_re])
      }

      if (sst) {
        R_u <- matern_cor(u, phi = phi, kappa = kappa, return_sym_matrix = TRUE)
        R_v <- matern_cor(v, phi = psi, kappa = 0.5, return_sym_matrix = TRUE)
        R <- R_u * R_v
      } else {
        R <- matern_cor(u, phi = phi, kappa=kappa,return_sym_matrix = TRUE)
      }
      diag(R) <- diag(R) + nu2
      R.inv <- solve(R)
      ldetR <- determinant(R)$modulus

      if(sst) {
        psi <- exp(par[ind_psi])
        R1.psi <- R_u * matern.grad.phi(v, psi, 0.5)
        m1.psi <- R.inv %*% R1.psi
        t1.psi <- -0.5 * sum(diag(m1.psi))
        m2.psi <- m1.psi %*% R.inv
      }

      exp.fact <- exp(compute.log.f(par, ldetR, R.inv) - log.f.tilde)
      L.m <- sum(exp.fact)
      exp.fact <- exp.fact / L.m

      if (sst) {
        R1.phi <- matern.grad.phi(u, phi, kappa) * R_v
      } else {
        R1.phi <- matern.grad.phi(u, phi, kappa)
      }
      m1.phi <- R.inv %*% R1.phi
      t1.phi <- -0.5 * sum(diag(m1.phi))
      m2.phi <- m1.phi %*% R.inv

      if (is.null(fix_tau2)) {
        t1.nu2 <- -0.5 * sum(diag(R.inv))
        m2.nu2 <- R.inv %*% R.inv
      }

      gradient.S <- function(S_tot) {
        S <- S_tot[1:n_loc]
        if (n_re > 0) {
          S_re_list <- list()
          for (i in 1:n_re) {
            S_re_list[[i]] <- S_tot[ind_re[[i]]]
          }
        }

        eta <- mu + S[ID_coords]
        if (n_re > 0) {
          for (i in 1:n_re) {
            eta <- eta + S_re_list[[i]][ID_re[, i]]
          }
        }

        prob_star <- 1 / (1 + exp(-eta))
        prob <- mda_effect * prob_star
        d_S <- (y / prob - (units_m - y) / (1 - prob)) * mda_effect * prob_star / (1 + exp(eta))
        d_S <- as.numeric(d_S)

        q.f_S <- t(S) %*% R.inv %*% S
        grad.beta <- t(D) %*% d_S
        grad.log.sigma2 <- (-n_loc / (2 * sigma2) + 0.5 * q.f_S / (sigma2^2)) * sigma2
        grad.log.phi <- (t1.phi + 0.5 * as.numeric(t(S) %*% m2.phi %*% S) / sigma2) * phi

        if (is.null(fix_alpha)) {
          der.alpha <- exp(par[ind_alpha]) / ((1 + exp(par[ind_alpha]))^2)
          grad.alpha.t <- der.alpha * (
            sum((y / mda_effect - (units_m - y) * prob_star / (1 - prob)) * mda_der_alpha) -
              penalty[[2]](alpha)
          )
        }

        grad.log.gamma <- gamma * sum((y / mda_effect - (units_m - y) * prob_star / (1 - prob)) * mda_der_gamma)

        out <- c(grad.beta, grad.log.sigma2, grad.log.phi)

        if (is.null(fix_tau2)) {
          grad.log.nu2 <- (t1.nu2 + 0.5 * as.numeric(t(S) %*% m2.nu2 %*% S) / sigma2) * nu2
          out <- c(out, grad.log.nu2)
        }

        if (n_re > 0) {
          grad.log.sigma2_re <- rep(NA, n_re)
          for (i in 1:n_re) {
            grad.log.sigma2_re[i] <- (-n_dim_re[i] / (2 * sigma2_re[i]) +
                                        0.5 * sum(S_re_list[[i]]^2) / (sigma2_re[i]^2)) * sigma2_re[i]
          }
          out <- c(out, grad.log.sigma2_re)
        }

        if (is.null(fix_alpha)) {
          out <- c(out, grad.alpha.t, grad.log.gamma)
        } else {
          out <- c(out, grad.log.gamma)
        }

        if (sst) {
          grad.log.psi <- (t1.psi + 0.5 * as.numeric(t(S) %*% m2.psi %*% S) / sigma2) * psi
          out <- c(out, grad.log.psi)
        }

        return(out)
      }

      out <- rep(0, length(par))
      for (i in 1:n_samples) {
        out <- out + exp.fact[i] * gradient.S(S_tot_samples[i, ])
      }
      out
    }

    hess.MC.log.lik <- function(par) {
      # Unpack parameters
      beta   <- par[ind_beta]
      mu     <- as.numeric(D %*% beta) + cov_offset
      sigma2 <- exp(par[ind_sigma2])
      if (is.null(fix_alpha)) {
        alpha <- exp(par[ind_alpha]) / (1 + exp(par[ind_alpha]))
      } else {
        alpha <- fix_alpha
      }
      gamma <- exp(par[ind_gamma])

      # MDA effect + derivatives
      mda_all           <- compute_mda_effect_derivatives(
        survey_times_data, mda_times,
        intervention = int_mat,
        alpha, gamma, kappa = power_val
      )
      mda_eff        <- mda_all$effect
      mda_der_alpha  <- mda_all$d_alpha
      mda_der_gamma  <- mda_all$d_gamma
      mda_der2_alpha <- mda_all$d2_alpha
      mda_der2_gamma <- mda_all$d2_gamma
      mda_der2_alpha_gamma <- mda_all$d2_alpha_gamma

      # nu2
      if (!is.null(fix_tau2)) {
        nu2 <- fix_tau2 / sigma2
      } else {
        nu2 <- exp(par[ind_nu2])
      }

      # Build R and its derivatives
      phi <- exp(par[ind_phi])
      if (sst) {
        psi    <- exp(par[ind_psi])
        R_u    <- matern_cor(u,   phi = phi, kappa = kappa, return_sym_matrix = TRUE)
        R_v    <- matern_cor(v,   phi = psi, kappa = 0.5,    return_sym_matrix = TRUE)
        R      <- R_u * R_v
        R1.phi <- matern.grad.phi(u, phi, kappa) * R_v
        R2.phi <- matern.hessian.phi(u, phi, kappa) * R_v
        R1.psi <- R_u * matern.grad.phi(v, psi, 0.5)
        R2.psi <- R_u * matern.hessian.phi(v, psi, 0.5)
      } else {
        R      <- matern_cor(u, phi = phi, kappa = kappa, return_sym_matrix = TRUE)
        R1.phi <- matern.grad.phi(u, phi, kappa)
        R2.phi <- matern.hessian.phi(u, phi, kappa)
      }
      diag(R) <- diag(R) + nu2
      R.inv   <- solve(R)

      # Precompute phi derivatives
      m1.phi <- R.inv %*% R1.phi
      t1.phi <- -0.5 * sum(diag(m1.phi))
      m2.phi <- m1.phi %*% R.inv
      t2.phi <- -0.5 * sum(diag(R.inv %*% R2.phi - R.inv %*% R1.phi %*% R.inv %*% R1.phi))
      n2.phi <- R.inv %*% (2 * R1.phi %*% R.inv %*% R1.phi - R2.phi) %*% R.inv

      if (sst) {
        R2.psi.phi <- matern.grad.phi(u, phi, 0.5)*matern.grad.phi(v, psi, 0.5)
        m1.psi   <- R.inv %*% R1.psi
        t1.psi   <- -0.5 * sum(diag(m1.psi))
        m2.psi   <- m1.psi %*% R.inv
        t2.psi   <- -0.5 * sum(diag(R.inv %*% R2.psi - R.inv %*% R1.psi %*% R.inv %*% R1.psi))
        n2.psi   <- R.inv %*% (2 * R1.psi %*% R.inv %*% R1.psi - R2.psi) %*% R.inv

        t2.psi.phi <- -0.5*(sum(R.inv*R2.psi.phi)-
                              sum(m1.phi*t(m1.psi)))
        n2.psi.phi <- R.inv%*%(R1.phi%*%m1.psi+
                                 R1.psi%*%m1.phi-
                                 R2.psi.phi)%*%R.inv

      }

      # Precompute nu2 derivatives
      if (is.null(fix_tau2)) {
        m2.nu2     <- R.inv %*% R.inv
        t1.nu2     <- -0.5 * sum(diag(R.inv))
        t2.nu2     <-  0.5 * sum(diag(m2.nu2))
        n2.nu2     <-  2 * R.inv %*% m2.nu2
        t2.nu2.phi <-  0.5 * sum(diag(R.inv %*% R1.phi %*% R.inv))
        n2.nu2.phi <-  R.inv %*% (R.inv %*% R1.phi + R1.phi %*% R.inv) %*% R.inv
        if (sst) {
          t2.nu2.psi <-  0.5 * sum(diag(R.inv %*% R1.psi %*% R.inv))
          n2.nu2.psi <-  R.inv %*% (R.inv %*% R1.psi + R1.psi %*% R.inv) %*% R.inv
        }
      }

      # Monte Carlo weights
      ldetR <- determinant(R)$modulus
      expf  <- exp(compute.log.f(par, ldetR, R.inv) - log.f.tilde)
      w     <- expf / sum(expf)

      # Accumulators
      H_acc <- matrix(0, nrow = length(par), ncol = length(par))
      g_acc <- rep(0, length(par))

      hessian.S <- function(S_tot, ef) {
        S <- S_tot[1:n_loc]
        if (n_re > 0) {
          S_re_list <- lapply(seq_len(n_re), function(i) S_tot[ind_re[[i]]])
        }

        # Linear predictor
        eta <- mu + S[ID_coords]
        if (n_re > 0) for (i in seq_len(n_re)) eta <- eta + S_re_list[[i]][ID_re[, i]]
        prob_star <- 1 / (1 + exp(-eta))
        prob      <- mda_eff * prob_star

        # Derivatives wrt S
        d_S  <- (y/prob - (units_m-y)/(1-prob)) * mda_eff * prob_star / (1+exp(eta))
        d2_S <- (-y/prob^2 - (units_m-y)/(1-prob)^2) * (mda_eff*prob_star/(1+exp(eta)))^2 +
          (y/prob - (units_m-y)/(1-prob)) * mda_eff * exp(eta)*(1-exp(eta))/(1+exp(eta))^3
        d_S  <- as.numeric(d_S)

        # d_S_alpha
        if (is.null(fix_alpha)) {
          d_S_alpha <- (-y/prob^2 - (units_m-y)/(1-prob)^2) *
            mda_der_alpha * prob_star * mda_eff * prob_star/(1+exp(eta)) +
            (y/prob - (units_m-y)/(1-prob)) * mda_der_alpha * prob_star/(1+exp(eta))
        }

        # d_S_gamma
        d_S_gamma <- (-y/prob^2 - (units_m-y)/(1-prob)^2) *
          mda_der_gamma * prob_star * mda_eff * prob_star/(1+exp(eta)) +
          (y/prob - (units_m-y)/(1-prob)) * mda_der_gamma * prob_star/(1+exp(eta))

        qSS <- as.numeric(t(S) %*% R.inv %*% S)

        # Gradient components
        grad_beta  <- t(D) %*% d_S
        grad_ls2   <- (-n_loc/(2*sigma2) + 0.5*qSS/sigma2^2) * sigma2
        grad_lphi  <- (t1.phi + 0.5*as.numeric(t(S)%*%m2.phi%*%S)/sigma2) * phi

        g_vec <- c(grad_beta, grad_ls2, grad_lphi)

        if (is.null(fix_tau2)) {
          grad_lnu2 <- (t1.nu2 + 0.5*as.numeric(t(S)%*%m2.nu2%*%S)/sigma2)*nu2
          g_vec     <- c(g_vec, grad_lnu2)
        }

        if (is.null(fix_alpha)) {
          der_alpha  <- exp(par[ind_alpha])/(1+exp(par[ind_alpha]))^2
          grad_alpha <- der_alpha * (
            sum((y/mda_eff - (units_m-y)*prob_star/(1-prob))*mda_der_alpha) -
              penalty[[2]](alpha)
          )
          g_vec <- c(g_vec, grad_alpha)
        }

        grad_lgamma <- gamma * sum((y/mda_eff - (units_m-y)*prob_star/(1-prob))*mda_der_gamma)
        g_vec       <- c(g_vec, grad_lgamma)

        if (n_re > 0) {
          grad_ls2re <- sapply(seq_len(n_re), function(i)
            (-n_dim_re[i]/(2*sigma2_re[i]) + 0.5*sum(S_re_list[[i]]^2)/sigma2_re[i]^2)*sigma2_re[i]
          )
          g_vec <- c(g_vec, grad_ls2re)
        }

        if (sst) {
          grad_lpsi <- (t1.psi + 0.5*as.numeric(t(S)%*%m2.psi%*%S)/sigma2)*psi
          g_vec     <- c(g_vec, grad_lpsi)
        }

        # Build H_loc
        H_loc <- matrix(0, nrow=length(par), ncol=length(par))

        # beta-beta
        H_loc[ind_beta, ind_beta] <- t(D) %*% (D * d2_S)

        # beta-alpha
        if (is.null(fix_alpha)) {
          H_loc[ind_beta, ind_alpha] <- t(D) %*% d_S_alpha * der_alpha
          H_loc[ind_alpha, ind_beta] <- H_loc[ind_beta, ind_alpha]
        }

        # beta-gamma
        H_loc[ind_beta, ind_gamma] <- t(D) %*% d_S_gamma * gamma
        H_loc[ind_gamma, ind_beta] <- H_loc[ind_beta, ind_gamma]

        # alpha-alpha
        if (is.null(fix_alpha)) {
          der2_alpha <- exp(par[ind_alpha])*(1-exp(par[ind_alpha]))/(1+exp(par[ind_alpha]))^3
          H_loc[ind_alpha, ind_alpha] <- der2_alpha * (
            sum((y/mda_eff - (units_m-y)*prob_star/(1-prob))*mda_der_alpha) -
              penalty[[2]](alpha)
          ) + der_alpha^2 * (
            sum((y/mda_eff - (units_m-y)*prob_star/(1-prob))*mda_der2_alpha +
                  (-y/mda_eff^2 - (units_m-y)*(prob_star^2)/(1-prob)^2)*mda_der_alpha^2
            ) - penalty[[3]](alpha)
          )
        }

        # alpha-gamma
        if (is.null(fix_alpha)) {
          H_loc[ind_alpha, ind_gamma] <- (gamma*der_alpha)*sum(
            (y/mda_eff - (units_m-y)*prob_star/(1-prob))*mda_der2_alpha_gamma +
              (-y/mda_eff^2 - (units_m-y)*(prob_star^2)/(1-prob)^2)*mda_der_alpha*mda_der_gamma
          )
          H_loc[ind_gamma, ind_alpha] <- H_loc[ind_alpha, ind_gamma]
        }

        # gamma-gamma
        H_loc[ind_gamma, ind_gamma] <- gamma*sum(
          (y/mda_eff - (units_m-y)*prob_star/(1-prob))*mda_der_gamma
        ) + gamma^2*sum(
          (y/mda_eff - (units_m-y)*prob_star/(1-prob))*mda_der2_gamma +
            (-y/mda_eff^2 - (units_m-y)*(prob_star^2)/(1-prob)^2)*mda_der_gamma^2
        )

        # sigma2-sigma2
        H_loc[ind_sigma2, ind_sigma2] <-
          (n_loc/(2*sigma2^2) - qSS/sigma2^3)*sigma2^2 + grad_ls2

        # sigma2-phi
        H_loc[ind_sigma2, ind_phi] <- (grad_lphi/phi - t1.phi)*(-phi)
        H_loc[ind_phi, ind_sigma2] <- H_loc[ind_sigma2, ind_phi]

        # phi-phi
        H_loc[ind_phi, ind_phi] <- (t2.phi - 0.5*as.numeric(t(S)%*%n2.phi%*%S)/sigma2)*phi^2 + grad_lphi

        # nu2-block
        if (is.null(fix_tau2)) {
          H_loc[ind_sigma2, ind_nu2] <- (grad_lnu2/nu2 - t1.nu2)*(-nu2)
          H_loc[ind_nu2, ind_sigma2] <- H_loc[ind_sigma2, ind_nu2]
          H_loc[ind_nu2, ind_nu2]   <- (t2.nu2 - 0.5*as.numeric(t(S)%*%n2.nu2%*%S)/sigma2)*nu2^2 + grad_lnu2
          H_loc[ind_phi, ind_nu2]   <- (t2.nu2.phi - 0.5*as.numeric(t(S)%*%n2.nu2.phi%*%S)/sigma2)*phi*nu2
          H_loc[ind_nu2, ind_phi]   <- H_loc[ind_phi, ind_nu2]
        }

        # psi-block
        if (sst) {
          H_loc[ind_psi, ind_psi]     <- (t2.psi - 0.5*as.numeric(t(S)%*%n2.psi%*%S)/sigma2)*psi^2 + grad_lpsi
          H_loc[ind_sigma2, ind_psi]    <- (grad_lpsi/psi - t1.psi)*(-psi)
          H_loc[ind_psi, ind_sigma2]    <- H_loc[ind_sigma2, ind_psi]
          H_loc[ind_phi, ind_psi]       <- (t2.psi.phi-0.5*t(S)%*%n2.psi.phi%*%S/sigma2)*phi*psi
          H_loc[ind_psi, ind_phi]       <- H_loc[ind_phi, ind_psi]
          if (is.null(fix_tau2)) {
            val <- (t2.nu2.psi - 0.5*as.numeric(t(S)%*%n2.nu2.psi%*%S)/sigma2)*nu2*psi
            H_loc[ind_nu2, ind_psi] <- val
            H_loc[ind_psi, ind_nu2] <- val
          }
        }

        list(
          mat1 = ef * (g_vec %*% t(g_vec) + H_loc),
          g    = ef * g_vec
        )
      }

      # Monte Carlo aggregation
      for (i in seq_len(n_samples)) {
        tmp   <- hessian.S(S_tot_samples[i, ], w[i])
        g_acc <- g_acc + tmp$g
        H_acc <- H_acc + tmp$mat1
      }

      H_acc - g_acc %*% t(g_acc)
    }

    start_cov_pars[-(1:2)] <- start_cov_pars[-(1:2)]/start_cov_pars[1]
    if(is.null(fix_alpha)) {
      start_par <- c(start_beta, log(start_cov_pars),
                     log(start_alpha/(1-start_alpha)), log(start_gamma))
    } else {
      start_par <- c(start_beta, log(start_cov_pars), log(start_gamma))
    }
    if(sst) start_par <- c(start_par, log(start_psi))


    out <- list()
    estim <- nlminb(
      start     = start_par,
      objective = function(x)   -MC.log.lik(x),
      gradient  = function(x)   -grad.MC.log.lik(x),
      hessian   = function(x)   -hess.MC.log.lik(x),
      control   = list(trace=1*messages)
    )

    out$estimate <- estim$par
    out$grad.MLE <- grad.MC.log.lik(estim$par)
    hess.MLE <- hess.MC.log.lik(estim$par)
    out$covariance <- solve(-hess.MLE)
    out$log.lik <- -estim$objective
    if(return_samples) out$S_samples <- S_tot_samples
    class(out) <- "RiskMap"
    return(out)
}


maxim.integrand.dast <- function(y,units_m,mu,Sigma,ID_coords, ID_re = NULL,
                                 mda_effect,
                            sigma2_re = NULL,
                            hessian=FALSE, gradient=FALSE) {
  # Sigma <- Sigma0
  # mda_effect <- mda_effect0
  # mu <- mu0
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
    prob_star <- 1/(1+exp(-eta))
    prob <- mda_effect*prob_star

    llik <- sum(y*log(prob)+(units_m-y)*log(1-prob))


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

    prob_star <- 1 / (1 + exp(-eta))
    prob <- mda_effect * prob_star

    # Compute derivative of log-likelihood with respect to eta
    d_S <- (y/prob - (units_m-y)/(1-prob))*mda_effect*prob_star/(1+exp(eta))
    d_S <- as.numeric(d_S)

    out <- rep(NA,n_tot)
    out[1:n_loc] <- as.numeric(-Sigma.inv%*%S+
                                 sapply(1:n_loc,function(i) sum(d_S[C_S[i,]])))
    if(n_re>0) {
      for(j in 1:n_re) {
        out[ind_re[[j]]] <- as.numeric(-S_re_list[[j]]/sigma2_re[[j]]+
                                         sapply(1:n_dim_re[[j]],
                                                function(x) sum(d_S[C_re[[j]][x,]])))
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

    prob_star <- 1 / (1 + exp(-eta))
    prob <- mda_effect * prob_star

    # Compute derivative of log-likelihood with respect to eta
    d2_S <- (-y/prob^2 - (units_m-y)/((1-prob)^2))*(mda_effect*prob_star/(1+exp(eta)))^2+
            (y/prob - (units_m-y)/(1-prob))*mda_effect*exp(eta)*(1-exp(2*eta))/((1+exp(eta))^4)
    d2_S <- as.numeric(d2_S)
    out <- matrix(0,nrow = n_tot, ncol = n_tot)

    out[1:n_loc, 1:n_loc] <-  -Sigma.inv
    diag(out[1:n_loc, 1:n_loc]) <- diag(out[1:n_loc, 1:n_loc])+
      sapply(1:n_loc,function(i) sum(d2_S[C_S[i,]]))
    if(n_re>0) {
      for(j in 1:n_re) {
        diag(out[ind_re[[j]], ind_re[[j]]]) <- -1/sigma2_re[j]
        diag(out[ind_re[[j]], ind_re[[j]]]) <- diag(out[ind_re[[j]], ind_re[[j]]])+
          sapply(1:n_dim_re[j],function(i) sum(d2_S[C_re[[j]][i,]]))

        out[1:n_loc,ind_re[[j]]]

        for(k in 1:n_dim_re[[j]]) {
          out[1:n_loc, ind_re[[j]]][,k] <- -sapply(1:n_loc,function(i) sum(d2_S[C_S_re[[j]][i,k,]]))
          out[ind_re[[j]], 1:n_loc][k,] <- out[1:n_loc,ind_re[[j]]][,k]
        }

        if(j < n_re) {
          counter <- 0
          for(w in (j+1):n_re) {
            counter <- counter + 1
            for(k in 1:n_dim_re[[w]]) {
              out[ind_re[[j]], ind_re[[w]]][,k] <- -sapply(1:n_dim_re[j],function(i) sum(d2_S[C_re_re[[j]][[counter]][i,k,]]))
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
Laplace_sampling_MCMC_dast <- function(y, units_m, mu, mda_effect, Sigma, ID_coords, ID_re = NULL,
                                       sigma2_re = NULL, control_mcmc,
                                       Sigma_pd = NULL, mean_pd = NULL, messages = TRUE) {

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # Precompute values that will not change
  Sigma.inv <- solve(Sigma)
  n_loc <- nrow(Sigma)
  n <- length(y)

  if ((!is.null(ID_re) & is.null(sigma2_re)) | (is.null(ID_re) & !is.null(sigma2_re))) {
    stop("To introduce unstructured random effects both `ID_re` and `sigma2_re` must be provided.")
  }

  n_re <- length(sigma2_re)
  n_tot <- if (n_re > 0) n_loc + sum(sapply(1:n_re, function(i) length(unique(ID_re[, i])))) else n_loc

  # Precompute the matrices
  C_S <- t(sapply(1:n_loc, function(i) ID_coords == i))

  if (n_re > 0) {
    ind_re <- vector("list", n_re)
    add_i <- 0
    n_dim_re <- sapply(1:n_re, function(i) length(unique(ID_re[, i])))

    for (i in 1:n_re) {
      ind_re[[i]] <- (add_i + n_loc + 1):(add_i + n_loc + n_dim_re[i])
      if (i < n_re) add_i <- sum(n_dim_re[1:i])
    }

    C_re <- vector("list", n_re)
    C_S_re <- vector("list", n_re)
    C_re_re <- vector("list", n_re)

    for (j in 1:n_re) {
      C_S_re[[j]] <- array(FALSE, dim = c(n_loc, n_dim_re[j], n))
      C_re[[j]] <- t(sapply(1:n_dim_re[j], function(i) ID_re[, j] == i))
      for (l in 1:n_dim_re[j]) {
        for (k in 1:n_loc) {
          ind_kl <- which(ID_coords == k & ID_re[, j] == l)
          if (length(ind_kl) > 0) C_S_re[[j]][k, l, ind_kl] <- TRUE
        }
      }

      if (j < n_re) {
        C_re_re[[j]] <- vector("list", n_re - j)
        counter <- 0
        for (w in (j + 1):n_re) {
          counter <- counter + 1
          C_re_re[[j]][[counter]] <- array(FALSE, dim = c(n_dim_re[j], n_dim_re[w], n))
          for (l in 1:n_dim_re[j]) {
            for (k in 1:n_dim_re[w]) {
              ind_lk <- which(ID_re[, j] == l & ID_re[, w] == k)
              if (length(ind_lk) > 0) {
                C_re_re[[j]][[counter]][l, k, ind_lk] <- TRUE
              }
            }
          }
        }
      }
    }
  }

  if (is.null(Sigma_pd) | is.null(mean_pd)) {
    out_maxim <- maxim.integrand.dast(y = y, units_m = units_m, Sigma = Sigma, mu = mu,
                                      mda_effect = mda_effect, ID_coords = ID_coords,
                                      ID_re = ID_re, sigma2_re = sigma2_re,
                                      hessian = FALSE, gradient = TRUE)

    if (is.null(Sigma_pd)) Sigma_pd <- out_maxim$Sigma.tilde
    if (is.null(mean_pd)) mean_pd <- out_maxim$mode
  }

  n_sim  <- control_mcmc$n_sim
  burnin <- control_mcmc$burnin
  thin   <- control_mcmc$thin
  h      <- control_mcmc$h %||% (1.65 / (n_tot^(1 / 6)))  # Default if not provided
  c1.h   <- control_mcmc$c1.h
  c2.h   <- control_mcmc$c2.h

  Sigma_pd_sroot <- t(chol(Sigma_pd))
  A <- solve(Sigma_pd_sroot)

  # Construct total covariance matrix (Sigma_tot)
  Sigma_tot <- if (n_re == 0) Sigma else {
    Sigma_tot <- matrix(0, n_tot, n_tot)
    Sigma_tot[1:n_loc, 1:n_loc] <- Sigma
    for (i in 1:n_re) diag(Sigma_tot)[ind_re[[i]]] <- sigma2_re[i]
    Sigma_tot
  }

  Sigma_w_inv <- solve(A %*% Sigma_tot %*% t(A))
  mu_w <- -as.numeric(A %*% mean_pd)

  cond.dens.W <- function(W, S_tot) {
    S <- S_tot[1:n_loc]

    eta <- mu + S[ID_coords]
    if (n_re > 0) {
      for (i in 1:n_re) {
        eta <- eta + S_tot[ind_re[[i]]][ID_re[, i]]
      }
    }

    prob_star <- stats::plogis(eta)             # in (0,1)
    prob <- mda_effect * prob_star              # expected in (0,1)
    eps <- .Machine$double.eps
    prob <- pmin(pmax(prob, eps), 1 - eps)      # clamp for numerical stability

    llik <- sum(y * log(prob) + (units_m - y) * log(1 - prob))
    diff_w <- W - mu_w

    as.numeric(-0.5 * crossprod(diff_w, Sigma_w_inv %*% diff_w) + llik)
  }

  lang.grad <- function(W, S_tot) {
    S <- S_tot[1:n_loc]

    eta <- mu + S[ID_coords]
    if (n_re > 0) {
      for (i in 1:n_re) {
        eta <- eta + S_tot[ind_re[[i]]][ID_re[, i]]
      }
    }

    prob_star <- stats::plogis(eta)
    prob <- mda_effect * prob_star
    eps <- .Machine$double.eps
    prob <- pmin(pmax(prob, eps), 1 - eps)

    # derivative wrt eta
    # d/deta prob = mda_effect * prob_star * (1 - prob_star)
    dprob_deta <- mda_effect * prob_star * (1 - prob_star)
    d_llik_deta <- (y / prob - (units_m - y) / (1 - prob)) * dprob_deta

    grad_S_tot_r <- rep(NA_real_, n_tot)
    # contributions to S (size n_loc)
    grad_S_tot_r[1:n_loc] <- as.numeric(sapply(1:n_loc, function(i) sum(d_llik_deta[C_S[i, ]])))

    if (n_re > 0) {
      for (j in 1:n_re) {
        # random effect prior term and likelihood term
        grad_like_re <- sapply(1:n_dim_re[[j]], function(x) sum(d_llik_deta[C_re[[j]][x, ]]))
        grad_S_tot_r[ind_re[[j]]] <- as.numeric(-S_tot[ind_re[[j]]] / sigma2_re[[j]] + grad_like_re)
      }
    }

    as.numeric(-Sigma_w_inv %*% (W - mu_w) + t(Sigma_pd_sroot) %*% grad_S_tot_r)
  }

  # MCMC loop
  W_curr <- rep(0, n_tot)
  S_tot_curr <- as.numeric(Sigma_pd_sroot %*% W_curr + mean_pd)
  mean_curr <- as.numeric(W_curr + (h^2 / 2) * lang.grad(W_curr, S_tot_curr))
  lp_curr <- cond.dens.W(W_curr, S_tot_curr)

  acc <- 0L
  n_samples <- floor((n_sim - burnin) / thin)
  sim <- matrix(NA_real_, nrow = n_samples, ncol = n_tot)

  if (messages) message("\n - Conditional simulation (burnin=", burnin, ", thin=", thin, "):")
  pb <- NULL
  if (messages && interactive()) {
    pb <- utils::txtProgressBar(min = 0, max = n_sim, style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }

  h.vec <- rep(NA_real_, n_sim)
  acc_prob <- rep(NA_real_, n_sim)

  for (i in 1:n_sim) {
    W_prop <- mean_curr + h * rnorm(n_tot)
    S_tot_prop <- as.numeric(Sigma_pd_sroot %*% W_prop + mean_pd)
    mean_prop <- as.numeric(W_prop + (h^2 / 2) * lang.grad(W_prop, S_tot_prop))
    lp_prop <- cond.dens.W(W_prop, S_tot_prop)

    dprop_curr <- -sum((W_prop - mean_curr)^2) / (2 * (h^2))
    dprop_prop <- -sum((W_curr - mean_prop)^2) / (2 * (h^2))

    log_prob <- lp_prop + dprop_prop - lp_curr - dprop_curr

    if (log(runif(1)) < log_prob) {
      acc <- acc + 1L
      W_curr <- W_prop
      S_tot_curr <- S_tot_prop
      lp_curr <- lp_prop
      mean_curr <- mean_prop
    }

    if (i > burnin && (i - burnin) %% thin == 0) {
      cnt <- (i - burnin) %/% thin
      sim[cnt, ] <- S_tot_curr
    }

    acc_prob[i] <- acc / i
    h <- max(1e-19, h + c1.h * i^(-c2.h) * (acc / i - 0.57))  # actually update h
    h.vec[i] <- h

    if (messages) {
      if (!is.null(pb)) {
        utils::setTxtProgressBar(pb, i)
      } else if (i %% max(1, floor(n_sim/20)) == 0) {
        message(sprintf("   %3d%%", round(100 * i / n_sim)))
      }
    }
  }

  if (!is.null(pb)) close(pb)
  if (messages) message(" done.\n")

  out_sim <- list(samples = list(S = sim[, 1:n_loc, drop = FALSE]))
  if (n_re > 0) {
    re_names <- if (!is.null(colnames(ID_re))) colnames(ID_re) else paste0("re", seq_len(n_re))
    for (i in 1:n_re) {
      out_sim$samples[[re_names[i]]] <- sim[, ind_re[[i]], drop = FALSE]
    }
  }

  out_sim$tuning_par <- h.vec
  out_sim$acceptance_prob <- acc_prob
  class(out_sim) <- "mcmc.RiskMap"

  return(out_sim)
}


