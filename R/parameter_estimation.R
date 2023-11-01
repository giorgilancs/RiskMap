##' @title Estimation of Generalized Linear Gaussian Process Models
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @importFrom sf st_crs st_as_sf st_drop_geometry
##' @export
glgpm <- function(formula,
                 data,
                 family,
                 distr_offset = NULL,
                 cov_offset = NULL,
                 crs = NULL, convert_to_crs = NULL,
                 scale_to_km = TRUE,
                 control_MCMC = NULL,
                 S_samples = NULL,
                 save_samples = F,
                 messages = TRUE,
                 fix_sigma2_me = NULL,
                 start_pars = list(beta = NULL,
                                   sigma2 = NULL,
                                   tau2 = NULL,
                                   phi = NULL,
                                   sigma2_me = NULL,
                                   sigma2_re = NULL)) {

  if(family=="binomial" | family=="poisson") {
    if(is.null(control_MCMC)) stop("if family='binomial' or family='poisson'
                                   'control_MCMC' must be provided.")
    if(class(control_MCMC)!="mcmc.RiskMap") stop("'control_MCMC' must be of class 'mcmc.PrevMap'")
  }

  if(class(formula)!="formula") stop("'formula' must be a 'formula'
                                     object indicating the variables of the
                                     model to be fitted")

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
      st_crs(data) <- crs
    }
  }


  if(family=="binomial") {
    if(is.null(distr_offset)) stop("if family='binomial', the argument 'm_offset'
                               must be provided")
    m_offset_ch <- as.charaster(as.name(subsitute(distr_offset)))
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


  # Extract coordinates
  if(!is.null(convert_to_crs)) {
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

  fix_tau2 <- inter_f$gp.spec$nugget

  if(all(table(ID_coords)==1) &
    is.null(fix_tau2) & is.null(fix_sigma2_me)) {
    stop("When there is only one observation per location, both the nugget and measurement error cannot
         be estimate. Consider removing either one of them. ")
  }

  if(scale_to_km) {
    coords_o <- coords_o/1000
    coords <- coords/1000
    if(messages) cat("Distances between locations are computed in kilometers \n")
  } else {
    if(messages) cat("Distances between locations are computed in meters \n")
  }


  if(is.null(start_pars$beta)) {
    start_pars$beta <- as.numeric(solve(t(D)%*%D)%*%t(D)%*%y)
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
    start_pars$phi <- quantile(dist(coords),0.25)
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
      if(length(sigma2_re)!=n_re) stop("starting values for 'sigma2_re' do not
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



  if(family=="gaussian") {
    if(is.null(fix_sigma2_me)) {
      if(is.null(start_pars$sigma2_me)) {
        start_pars$sigma2_me <- 1
      } else {
        if(start_pars$sigma2_me<0) stop("the starting value for sigma2_me must be positive")
      }
    }
    res <- glgpm_lm(y, D, coords, kappa = inter_f$gp.spec$kappa,
            ID_coords, ID_re, s_unique, re_unique,
            fix_sigma2_me, fix_tau2,
            start_beta = start_pars$beta,
            start_cov_pars = c(start_pars$sigma2,
                               start_pars$phi,
                               start_pars$tau2,
                               start_pars$sigma2_re,
                               start_pars$sigma2_me),
            messages = messages)
  } else if(family=="binomial" | family=="poisson") {

  }

  res$y <- y
  res$D <- D
  res$coords <- coords
  res$ID_coords <- ID_coords
  res$re <- re_unique_f
  if(n_re>1) {
    res$ID_re <- as.data.frame(ID_re)
    colnames(res$ID_re) <- names_re
  }
  res$fix_tau2 <- fix_tau2
  res$fix_sigma2_me <- fix_sigma2_me
  res$formula <- formula
  res$family <- family
  res$crs <- crs
  res$convert_to_crs <- convert_to_crs
  res$scale_to_km <- scale_to_km
  res$data_sf <- data
  res$kappa <- kappa
  res$call <- match.call()
  return(res)
}


##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @importFrom maxLik maxBFGS
##' @importFrom Matrix Matrix forceSymmetric
glgpm_lm <- function(y, D, coords, kappa, ID_coords, ID_re, s_unique, re_unique,
                    fix_sigma2_me, fix_tau2, start_beta, start_cov_pars, messages) {

  m <- length(y)
  n <- nrow(coords)
  p <- ncol(D)
  U <- dist(coords)
  if(is.null(ID_re)) {
    n_re <- 0
  } else {
    n_re <- ncol(ID_re)
  }

  if(!is.null(fix_sigma2_me)) {
    if(fix_sigma2_me==0) {
      fix_sigma2_me <- 10e-10
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

    if(!is.null(fix_sigma2_me)) {
      omega2 <- fix_sigma2_me
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
    if(!is.null(fix_sigma2_me)) {
      omega2 <- fix_sigma2_me
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


    if(is.null(fix_sigma2_me)) {
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
    if(!is.null(fix_sigma2_me)) {
      omega2 <- fix_sigma2_me
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


    if(is.null(fix_sigma2_me)) {
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

    if(is.null(fix_sigma2_me)) {
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

    if(is.null(fix_sigma2_me)) {
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

      if(is.null(fix_sigma2_me)) {
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

    if(is.null(fix_sigma2_me)) {
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
    if(is.null(object$fix_sigma2_me)) {
      ind_sigma2_me <- p+4
    } else {
      ind_sigma2_me <- NULL
    }
    if(n_re>0) {
      ind_sigma2_re <- (p+5):(p+4+n_re)
    }
  } else {
    ind_tau2 <- NULL
    if(is.null(object$fix_sigma2_me)) {
      ind_sigma2_me <- p+3
      names(object$estimate)[ind_sigma2_me] <- "sigma2_me"
    } else {
      ind_sigma2_me <- NULL
    }
    if(n_re>0) {
      ind_sigma2_re <- (p+4):(p+3+n_re)
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

  res <- object$estimate
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
    if(is.null(object$fix_sigma2_me)) {
      ind_sigma2_me <- p+4
    } else {
      ind_sigma2_me <- NULL
    }
    if(n_re>0) {
      ind_sigma2_re <- (p+5):(p+4+n_re)
    }
  } else {
    ind_tau2 <- NULL
    if(is.null(object$fix_sigma2_me)) {
      ind_sigma2_me <- p+3
      names(object$estimate)[ind_sigma2_me] <- "Measuremment error var."
    } else {
      ind_sigma2_me <- NULL
    }
    if(n_re>0) {
      ind_sigma2_re <- (p+4):(p+3+n_re)
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
    if(is.null(object$fix_sigma2_me)) {
      res$me <- cbind(Estimate = object$estimate[ind_sigma2_me],
                      'Lower limit' = exp(log(object$estimate[ind_sigma2_me])-qnorm(1-alpha/2)*se_par[ind_sigma2_me]),
                      'Upper limit' = exp(log(object$estimate[ind_sigma2_me])+qnorm(1-alpha/2)*se_par[ind_sigma2_me]))
    } else {
      res$me <- object$fix_sigma2_me
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
  res$aic <- 2*length(res$estimate)-2*res$log.lik

  class(res) <- "summary.RiskMap"
  return(res)
}

##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @method print summary.RiskMap
##' @export
print.summary.RiskMap <- function(x) {
  if(x$family=="gaussian") {
    cat("Geostatistical linear model \n")
  }

  cat("'Lower limit' and 'Upper limit' are the ",
      x$conf_level*100,
      "% confidence level intervals \n", sep="")


  cat("\n Regression coefficients \n")
  printCoefmat(x$reg_coef,P.values=TRUE,has.Pvalue=TRUE)

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
  cat("\n AIC: ",x$aic,"\n \n",sep="")

}

##' @title Simulation from Generalized Linear Gaussian Process Models
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @importFrom sf st_crs st_as_sf st_drop_geometry
##' @export
glgpm_sim <- function(n_sim,
                      model_fit = NULL,
                      formula = NULL,
                      data = NULL,
                      family = NULL,
                      distr_offset = NULL,
                      cov_offset = NULL,
                      crs = NULL, convert_to_crs = NULL,
                      scale_to_km = TRUE,
                      control_MCMC = NULL,
                      sim_pars = list(beta = NULL,
                                      sigma2 = NULL,
                                      tau2 = NULL,
                                      phi = NULL,
                                      sigma2_me = NULL,
                                      sigma2_re = NULL),
                      messages = TRUE) {

  if(!is.null(model_fit)) {
    if(class(model_fit)!="RiskMap") stop("'model_fit' must be of class 'RiskMap'")
    formula <- as.formula(model_fit$formula)
    data <- model_fit$data_sf
    family = model_fit$family
    crs <- model_fit$crs
    convert_to_crs <- model_fit$convert_to_crs
    scale_to_km <- model_fit$scale_to_km
  }
  inter_f <- interpret.formula(formula)

  if(family=="binomial" | family=="poisson") {
    if(is.null(control_MCMC)) stop("if family='binomial' or family='poisson'
                                   'control_MCMC' must be provided.")
    if(class(control_MCMC)!="mcmc.RiskMap") stop("'control_MCMC' must be of class 'mcmc.PrevMap'")
  }

  if(class(formula)!="formula") stop("'formula' must be a 'formula'
                                     object indicating the variables of the
                                     model to be fitted")

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
      st_crs(data) <- crs
    }
  }


  if(family=="binomial") {
    if(is.null(distr_offset)) stop("if family='binomial', the argument 'm_offset'
                               must be provided")
    m_offset_ch <- as.charaster(as.name(subsitute(distr_offset)))
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
      if(is.null(model_fit$fix_sigma2_me)) {
        ind_sigma2_me <- p+4
        sigma2_me <- par_hat[ind_sigma2_me]
      } else {
        sigma2_me <- model_fit$fix_sigma2_me
      }
      if(n_re>0) {
        ind_sigma2_re <- (p+5):(p+4+n_re)
        sigma2_re <- par_hat[ind_sigma2_re]
      }
    } else {
      tau2 <- model_fit$fix_tau2
      if(is.null(model_fit$fix_sigma2_me)) {
        ind_sigma2_me <- p+3
        sigma2_me <- par_hat[ind_sigma2_me]
      } else {
        sigma2_me <- model_fit$fix_sigma2_me
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
      if(is.null(sim_pars$sigma2_re)) steop("'sigma2_re' is missing")
      if(length(sim_pars$sigma2_re)!=n_re) stop("the values passed to 'sigma2_re' in 'sim_pars'
      does not match the number of random effects specfied in re() in the formula")
      sigma2_re <- sim_pars$sigma2_re
    }
  }

  # Extract coordinates
  if(!is.null(convert_to_crs)) {
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




  if(all(table(ID_coords)==1) & (tau2!=0 & sigma2_me!=0)) {
    stop("When there is only one observation per location, both the nugget and measurement error cannot
         be estimate. Consider removing either one of them. ")
  }

  if(scale_to_km) {
    coords_o <- coords_o/1000
    coords <- coords/1000
    if(messages) cat("Distances between locations are computed in kilometers \n")
  } else {
    if(messages) cat("Distances between locations are computed in meters \n")
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
