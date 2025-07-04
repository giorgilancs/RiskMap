##' @importFrom stats setNames
##' @importFrom utils head

##' @title Prediction of the random effects components and covariates effects over a spatial grid using a fitted generalized linear Gaussian process model
##'
##' @description This function computes predictions over a spatial grid using a fitted model
##' obtained from the \code{\link{glgpm}} function. It provides point predictions and uncertainty
##' estimates for the specified locations for each component of the model separately: the spatial random effects;
##' the unstructured random effects (if included); and the covariates effects.
##'
##' @param object A RiskMap object obtained from the `glgpm` function.
##' @param grid_pred An object of class 'sfc', representing the spatial grid over which predictions
##'                 are to be made. Must be in the same coordinate reference system (CRS) as the
##'                 object passed to 'object'.
##' @param predictors Optional. A data frame containing predictor variables used for prediction.
##' @param re_predictors Optional. A data frame containing predictors for unstructured random effects,
##'                      if applicable.
##' @param pred_cov_offset Optional. A numeric vector specifying covariate offsets at prediction locations.
##' @param control_sim Control parameters for MCMC sampling. Must be an object of class "mcmc.RiskMap" as returned by \code{\link{set_control_sim}}.
##' @param type Type of prediction. "marginal" for marginal predictions, "joint" for joint predictions.
##' @param messages Logical. If TRUE, display progress messages. Default is TRUE.
##'
##' @return An object of class 'RiskMap.pred.re' containing predicted values, uncertainty estimates,
##'         and additional information.
##'
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##' @importFrom Matrix solve
##' @export
##'

pred_over_grid <- function(object,
                   grid_pred,
                   predictors = NULL,
                   re_predictors = NULL,
                   pred_cov_offset = NULL,
                   control_sim = set_control_sim(),
                   type = "marginal",
                   messages = TRUE) {
  if(!inherits(object,
               what = "RiskMap", which = FALSE)) {
    stop("The object passed to 'object' must be an output of
         the function 'glgpm'")
  }

  if(!inherits(grid_pred,
               what = "sfc", which = FALSE)) {
    stop("The object passed to 'grid_pred' must be an object
         of class 'sfc'")
  }

  if(!inherits(control_sim,
                 what = "mcmc.RiskMap", which = FALSE)) {
      stop ("the argument passed to 'control_sim' must be an output
                                                  from the function set_control_sim; see ?set_control_sim
                                                  for more details")

  }

  grid_pred <- st_transform(grid_pred,crs = object$crs)

  grp <- st_coordinates(grid_pred)
  n_pred <- nrow(grp)
  par_hat <- coef(object)
  object$D <- as.matrix(object$D)
  p <- ncol(object$D)

  if(!all(length(object$cov_offset==0))) {
    if(is.null(pred_cov_offset)) {
      stop("The covariate offset must be specified at each of the prediciton
           locations.")
    } else{
      if(!inherits(pred_cov_offset,
                   what = "numeric", which = FALSE)) {
        stop("'pred_cov_offset' must be a numeric vector")
      }
      if(length(pred_cov_offset) != n_pred) {
        stop("The length of 'pred_cov_offset' does not match the number of
             provided prediction locations")
      }
    }
  } else {
    pred_cov_offset <- 0
  }

  if(type!="marginal" & type!="joint") {
    stop("the argument 'type' must be set to 'marginal' or 'joint'")
  }

  inter_f <- interpret.formula(object$formula)
  inter_lt_f <- inter_f
  inter_lt_f$pf <- update(inter_lt_f$pf, NULL ~.)

  if(p==1) {
    intercept_only <- prod(object$D[,1]==1)==1
  } else {
    intercept_only <- FALSE
  }

  n_re <- length(object$re)
  if(n_re > 0 & type=="marginal" & !is.null(re_predictors)) {
    stop("Predictions for the unstructured random effects will not be perfomed if type = 'marginal';
            if you wish to also include the predictions for the unstructured random
            effects then set type = 'joint'")
  }

  if(!is.null(predictors)) {

    if(!is.data.frame(predictors)) stop("'predictors' must be an object of class 'data.frame'")
    if(nrow(predictors)!=n_pred) stop("the values provided for 'predictors' do not match the prediction grid passed to 'grid_pred'")

    if(any(is.na(predictors))) {
      warning("There are missing values in 'predictors'; these values have been removed
              alongside the corresponding prediction locations")
      if(any(is.na(re_predictors))) {
        warning("There are missing values in 're_predictors'; these values have been removed
              alongside the corresponding prediction locations")
      }
      if(!is.null(re_predictors)) {
        if(nrow(re_predictors)!=n_pred) stop("the values provided for 're_predictors' do not match the prediction grid passed to 'grid_pred'")

        comb_pred <- data.frame(predictors,
                                re_predictors)
        ind_c <- complete.cases(as.data.frame(st_drop_geometry(comb_pred)))
        predictors_aux <- data.frame(na.omit(comb_pred)[,1:ncol(predictors)])
        colnames(predictors_aux) <- colnames(predictors)
        predictors <- predictors_aux

        re_predictors_aux <- data.frame(na.omit(comb_pred)[,(ncol(predictors)+1):
                                                            (ncol(predictors)+
                                                             ncol(re_predictors))])
        colnames(re_predictors_aux) <- colnames(re_predictors_aux)
        re_predictors <- re_predictors_aux
      } else {
        comb_pred <- predictors
        ind_c <- complete.cases(as.data.frame(st_drop_geometry(comb_pred)))
        predictors_aux <- data.frame(na.omit(comb_pred))
        colnames(predictors_aux) <- colnames(predictors)
        predictors <- predictors_aux
      }
      grid_pred_aux <- st_coordinates(grid_pred)[ind_c,]
      grid_pred <- st_as_sf(data.frame(grid_pred_aux),
                            coords = c("X","Y"),
                            crs = st_crs(grid_pred)$input)
      grp <- st_coordinates(grid_pred)
      n_pred <- nrow(grp)

    }

    mf_pred <- model.frame(inter_lt_f$pf,data=predictors, na.action = na.fail)
    D_pred <- as.matrix(model.matrix(attr(mf_pred,"terms"),data=predictors))
    if(ncol(D_pred)!=ncol(object$D)) stop("the provided variables in 'predictors' do not match the number of explanatory variables used to fit the model.")
    mu_pred <- as.numeric(D_pred%*%par_hat$beta)

  } else if(intercept_only)  {
    mu_pred <- par_hat$beta
  } else {
    mu_pred <- 0
  }

  if(n_re>0) {

    ID_g <- as.matrix(cbind(object$ID_coords, object$ID_re))
    re_unique <- object$re
    n_dim_re_tot <- sapply(1:(n_re+1), function(i) length(unique(ID_g[,i])))
    if(!is.null(re_predictors)) {
      if(any(is.na(re_predictors))) {
        warning("There are missing values in 're_predictors'; these values have been removed
              alongside the corresponding prediction locations")
        comb_pred <- data.frame(re_predictors)
        ind_c <- complete.cases(comb_pred)
        re_predictors_aux <- data.frame(na.omit(comb_pred))
        colnames(re_predictors_aux) <- colnames(re_predictors_aux)
        re_predictors <- re_predictors_aux
        grid_pred_aux <- st_coordinates(grid_pred)[ind_c,]
        grid_pred <- st_as_sf(data.frame(grid_pred_aux),
                              coords = c("X","Y"),
                              crs = st_crs(grid_pred)$input)
        grp <- st_coordinates(grid_pred)
        n_pred <- nrow(grp)
      }
      if(!is.data.frame(re_predictors)) stop("'re_predictors' must be an object of class 'data.frame'")

      if(nrow(re_predictors)!=n_pred) stop("the values provided for 're_predictors' do not match the prediction grid passed to 'grid_pred'")
      if(ncol(re_predictors)!=n_re) stop("the number of unstructured random effects provided in 're_predictors' does not match that of the fitted model")

      D_re_pred <- list()
      n_dim_re <- sapply(1:n_re, function(i) length(unique(object$ID_re[,i])))
      for(i in 1:n_re) {
        re_val_i <- re_predictors[,i]
        D_re_pred[[i]] <- matrix(0,nrow = n_pred, ncol = n_dim_re[i])
        for(j in 1:length(object$re[[i]])) {
          ind_j <- which(re_val_i==object$re[[i]][j])
          D_re_pred[[i]][ind_j, j] <- 1
        }
      }
    } else {
      D_re_pred <- NULL
    }
  } else  {
    D_re_pred <- NULL
  }

  out <- list()

  out$mu_pred <- mu_pred
  out$grid_pred <- grid_pred
  out$par_hat <- object$par_hat

  if(object$scale_to_km) grp <- grp/1000
  if(object$family=="gaussian") {
    U_pred <- t(sapply(1:n_pred,
                       function(i) sqrt((object$coords[object$ID_coords,1]-grp[i,1])^2+
                                        (object$coords[object$ID_coords,2]-grp[i,2])^2)))
  } else {
    U_pred <- t(sapply(1:n_pred,
                       function(i) sqrt((object$coords[,1]-grp[i,1])^2+
                                          (object$coords[,2]-grp[i,2])^2)))
  }
  U <- dist(object$coords)
  C <- par_hat$sigma2*matern_cor(U_pred, phi = par_hat$phi,
                                 kappa = object$kappa)
  mu <- as.numeric(object$D%*%par_hat$beta)

  if(control_sim$linear_model) {
    n_samples <- control_sim$n_sim
  } else {
    n_samples <- (control_sim$n_sim-control_sim$burnin)/control_sim$thin
  }

  R <- matern_cor(U,phi = par_hat$phi, kappa=object$kappa,return_sym_matrix = TRUE)
  diff.y <- object$y-mu
  if(!is.null(object$fix_tau2)) {
    nu2 <- object$fix_tau2/par_hat$sigma2
  } else {
    nu2 <- par_hat$tau2/par_hat$sigma2
  }
  if(nu2==0) nu2 <- 10e-10

  diag(R) <- diag(R)+nu2

  dast_model <- !is.null(object$power_val)
  if(object$family!="gaussian") {
    Sigma <- par_hat$sigma2*R
    Sigma_inv <- solve(Sigma)
    A <- C%*%Sigma_inv

    if(dast_model) {
      alpha <- par_hat$alpha
      if(is.null(alpha)) alpha <- object$fix_alpha

      gamma <- par_hat$gamma

      mda_effect <- compute_mda_effect(object$survey_times_data, object$mda_times,
                                        object$int_mat,
                                        alpha, gamma, kappa = object$power_val)
      simulation <-
        Laplace_sampling_MCMC_dast(y = object$y, units_m = object$units_m, mu = mu, Sigma = Sigma,
                              sigma2_re = par_hat$sigma2_re,
                              mda_effect = mda_effect,
                              ID_coords = object$ID_coords, ID_re = object$ID_re,
                              control_mcmc = control_sim,
                              messages = messages)
    } else {
      simulation <-
        Laplace_sampling_MCMC(y = object$y, units_m = object$units_m, mu = mu, Sigma = Sigma,
                              sigma2_re = par_hat$sigma2_re,
                              ID_coords = object$ID_coords, ID_re = object$ID_re,
                              family = object$family, control_mcmc = control_sim,
                              messages = messages)
    }
    mu_cond_S <- A%*%t(simulation$samples$S)
    if(type=="marginal") {
      sd_cond_S <- sqrt(par_hat$sigma2-diag(A%*%t(C)))
      out$S_samples <- sapply(1:n_samples,
                              function(i)
                                mu_cond_S[,i]+
                                sd_cond_S*rnorm(n_pred))
    } else {
      U_pred_o <- dist(grp)
      Sigma_pred <-  par_hat$sigma2*matern_cor(U_pred_o, phi = par_hat$phi,
                                               kappa = object$kappa,
                                               return_sym_matrix = TRUE)
      Sigma_cond <- Sigma_pred - A%*%t(C)
      Sigma_cond_sroot <- t(chol(Sigma_cond))
      out$S_samples <- sapply(1:n_samples,
                              function(i)
                                mu_cond_S[,i]+
                                Sigma_cond_sroot%*%rnorm(n_pred))
    }
  } else {

    if(!is.null(object$fix_var_me) && object$fix_var_me>0 ||
       is.null(object$fix_var_me)) {
      m <- length(object$y)
      s_unique <- unique(object$ID_coords)

      ID_g <- as.matrix(cbind(object$ID_coords, object$ID_re))

      n_dim_re_tot <- sapply(1:(n_re+1), function(i) length(unique(ID_g[,i])))
      C_g <- matrix(0, nrow = m, ncol = sum(n_dim_re_tot))

      for(i in 1:m) {
        ind_s_i <- which(s_unique==ID_g[i,1])
        C_g[i,1:n_dim_re_tot[1]][ind_s_i] <- 1
      }

      if(n_re>0) {
        for(j in 1:n_re) {
          select_col <- sum(n_dim_re_tot[1:j])

          for(i in 1:m) {
            ind_re_j_i <- which(re_unique[[j]]==ID_g[i,j+1])
            C_g[i,select_col+1:n_dim_re_tot[j+1]][ind_re_j_i] <- 1
          }
        }
      }
      C_g <- Matrix(C_g, sparse = TRUE, doDiag = FALSE)
      C_g_m <- Matrix::t(C_g)%*%C_g
      C_g_m <- forceSymmetric(C_g_m)


      Sigma_g <- matrix(0, nrow = sum(n_dim_re_tot), ncol = sum(n_dim_re_tot))
      Sigma_g_inv <- matrix(0, nrow = sum(n_dim_re_tot), ncol = sum(n_dim_re_tot))
      Sigma_g[1:n_dim_re_tot[1], 1:n_dim_re_tot[1]] <- par_hat$sigma2*R
      Sigma_g_inv[1:n_dim_re_tot[1], 1:n_dim_re_tot[1]] <-
        solve(R)/par_hat$sigma2
      if(n_re > 0) {
        for(j in 1:n_re) {
          select_col <- sum(n_dim_re_tot[1:j])

          diag(Sigma_g[select_col+1:n_dim_re_tot[j+1], select_col+1:n_dim_re_tot[j+1]]) <-
            par_hat$sigma2_re[j]

          diag(Sigma_g_inv[select_col+1:n_dim_re_tot[j+1], select_col+1:n_dim_re_tot[j+1]]) <-
            1/par_hat$sigma2_re[j]

        }
      }

      Sigma_star <- Sigma_g_inv+C_g_m/par_hat$sigma2_me
      Sigma_star_inv <- forceSymmetric(Matrix::solve(Sigma_star))

      B <- -C_g%*%Sigma_star_inv%*%Matrix::t(C_g)/(par_hat$sigma2_me^2)
      diag(B) <- Matrix::diag(B) + 1/par_hat$sigma2_me
      A <- C%*%B
    } else {
      Sigma <- par_hat$sigma2*R
      Sigma_inv <- solve(Sigma)
      A <- C%*%Sigma_inv
    }

    mu_cond_S <- as.numeric(A%*%diff.y)

    if(type=="marginal") {
      sd_cond_S <- sqrt(par_hat$sigma2-Matrix::diag(A%*%t(C)))
      out$S_samples <- sapply(1:n_samples,
                              function(i)
                                mu_cond_S+
                                sd_cond_S*rnorm(n_pred))
    } else {
      U_pred_o <- dist(grp)
      Sigma_pred <-  par_hat$sigma2*matern_cor(U_pred_o, phi = par_hat$phi,
                                               kappa = object$kappa,
                                               return_sym_matrix = TRUE)

      Sigma_cond <- Sigma_pred - A%*%t(C)
      Sigma_cond_sroot <- t(chol(Sigma_cond))
      out$S_samples <- sapply(1:n_samples,
                              function(i)
                                mu_cond_S+
                                Sigma_cond_sroot%*%rnorm(n_pred))
    }
  }

  if(n_re>0 & !is.null(re_predictors)) {
    out$re <- list()
    out$re$D_pred <- D_re_pred
    out$re$samples <- list()
    re_names <- colnames(object$ID_re)
    if(object$family=="gaussian") {
      Sigma_cond_inv <- solve(Sigma_cond)
      C_Z <- C_g[,-(1:n_dim_re_tot[1])]
      add <- 0
      for(i in 1:length(n_dim_re_tot[-1])) {
        C_Z[,add+1:n_dim_re_tot[i+1]] <- par_hat$sigma2_re[i]*C_Z[,add+1:n_dim_re_tot[i+1]]
        add <- n_dim_re_tot[i+1]
      }

      A_Z <- Matrix::t(C_Z)%*%B%*%t(C)%*%Sigma_cond_inv

      Sigma_Z_cond <- diag(rep(par_hat$sigma2_re,n_dim_re_tot[-1]))-
        Matrix::t(C_Z)%*%B%*%C_Z-
        A_Z%*%C%*%Matrix::t(B)%*%C_Z
      Sigma_Z_cond_sroot <- t(chol(Sigma_Z_cond))

      mu_Z_cond <- sapply(1:n_samples, function(i) as.matrix(A_Z%*%(out$S_samples[,i]-mu_cond_S)))
      add <- 0
      for(i in 1:n_re) {
        for(j in 1:n_dim_re_tot[1+i]) {
          add <- add + 1
          C_re_ij <- matrix(0, ncol= m)
          ind_ij <- which(object$ID_re[[i]]==re_unique[[i]][j])
          C_re_ij[,ind_ij] <- par_hat$sigma2_re[i]
          A_re <- C_re_ij%*%B

          mu_Z_cond[add,] <- as.numeric(A_re%*%diff.y)+mu_Z_cond[add,]
        }
      }
      re_samples <- sapply(1:n_samples,
                           function(i)
                             as.numeric(
                               mu_Z_cond[,i]+
                                 Sigma_Z_cond_sroot%*%rnorm(sum(n_dim_re_tot[-1]))))
    } else {
      re_samples <- matrix(0, nrow = sum(n_dim_re_tot[-1]), ncol = n_samples)
      add <- 0
      for(i in 1:n_re) {
        re_samples[1:n_dim_re_tot[i+1],] <- t(simulation$samples[[i+1]])
        add <- add+n_dim_re_tot[i+1]
      }
    }

    add <- 0
    for(i in 1:n_re) {
      for(j in 1:n_dim_re_tot[1+i]) {
        add <- add + 1
        out$re$samples[[paste(re_names[[i]])]][[paste(re_unique[[i]][[j]])]] <-
          re_samples[add,]
      }
    }

  } else {
    out$re <- list(D_pred = NULL, samples = NULL)
  }
  if(dast_model) {
    out$mda_times <- object$mda_times
    if(is.null(par_hat$alpha)) out$fix_alpha <- object$fix_alpha
    out$power_val <- object$power_val
  }
  out$inter_f <- inter_f
  out$family <- object$family
  out$par_hat <- par_hat
  out$cov_offset <- pred_cov_offset
  out$type <- type
  class(out) <- "RiskMap.pred.re"
  return(out)
}

##' @title Predictive Target Over a Regular Spatial Grid
##'
##' @description Computes predictions over a regular spatial grid using outputs from the
##' \code{\link{pred_over_grid}} function.
##' This function allows for incorporating covariates, offsets, MDA effects, and optional
##' unstructured random effects into the predictive target.
##'
##' @param object Output from `pred_over_grid`, a RiskMap.pred.re object.
##' @param include_covariates Logical. Include covariates in the predictive target.
##' @param include_nugget Logical. Include the nugget effect in the predictive target.
##' @param include_cov_offset Logical. Include the covariate offset in the predictive target.
##' @param include_mda_effect Logical. Include the MDA effect in the predictive target using a DAST model; see \code{\link{dast}}.
##' @param mda_grid Optional. Grid of MDA coverage values required for predictions using a DAST model; see \code{\link{dast}}.
##' @param time_pred Optional. Time point for prediction required for predictions using a DAST model; see \code{\link{dast}}.
##' @param include_re Logical. Include unstructured random effects in the predictive target.
##' @param f_target Optional. List of functions to apply on the linear predictor samples.
##' @param pd_summary Optional. List of summary functions to apply on the predicted values.
##'
##' @return An object of class 'RiskMap_pred_target_grid' containing predicted values
##'         and summaries over the regular spatial grid.
##' @seealso \code{\link{pred_over_grid}}
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##' @importFrom Matrix solve
##' @export
pred_target_grid <- function(object,
                        include_covariates = TRUE,
                        include_nugget = FALSE,
                        include_cov_offset = FALSE,
                        include_mda_effect = TRUE,
                        mda_grid = NULL,
                        time_pred = NULL,
                        include_re = FALSE,
                        f_target = NULL,
                        pd_summary = NULL) {
  if(!inherits(object,
               what = "RiskMap.pred.re", which = FALSE)) {
    stop("The object passed to 'object' must be an output of
         the function 'pred_over_grid'")
  }

  dast_model <- !is.null(object$par_hat$gamma)

  if(dast_model) {
    if(include_mda_effect & is.null(mda_grid)) {
      stop("The MDA coverage must be specified for each point on the grid through the argument 'mda_grid'")
    }
    if(is.null(time_pred)) {
      stop("For a DAST model, the time of prediction must be specified through the argument 'time_pred'")
    }
  }


  if(is.null(object$par_hat$tau2) &
     include_nugget) {
    stop("The nugget cannot be included in the predictive target
             because it was not included when fitting the model")
  }

  if(is.null(f_target)) {
    f_target <- list(linear_target = function(x) x)
  }

  if(is.null(pd_summary)) {
    pd_summary <- list(mean = mean, sd = sd)
  }

  n_summaries <- length(pd_summary)
  n_f <- length(f_target)
  n_samples <- ncol(object$S_samples)
  n_pred <- nrow(object$S_samples)

  n_re <- length(object$re$samples)
  re_names <- names(object$re$samples)

  out <- list()
  if(length(object$mu_pred)==1 && object$mu_pred==0 &&
     include_covariates) {
    stop("Covariates have not been provided; re-run pred_over_grid
         and provide the covariates through the argument 'predictors'")
  }

  if(n_re==0 &&
     include_re) {
    stop("The categories of the randome effects variables have not been provided;
         re-run pred_over_grid and provide the covariates through the argument 're_predictors'")
  }

  if(!include_covariates) {
    mu_target <- 0
  } else {

    if(is.null(object$mu_pred)) stop("the output obtained from 'pred_S' does not
                                     contain any covariates; if including covariates
                                     in the predictive target these shuold be included
                                     when running 'pred_S'")
    mu_target <- object$mu_pred
  }

  if(!include_cov_offset) {
    cov_offset <- 0
  } else {
    if(length(object$cov_offset)==1) {
      stop("No covariate offset was included in the model;
           set include_cov_offset = FALSE, or refit the model and include
           the covariate offset")
    }
    cov_offset <- object$cov_offset
  }

  if(include_nugget) {
    Z_sim <- matrix(rnorm(n_samples*n_pred,
                          sd = sqrt(object$par_hat$tau2)),
                    ncol = n_samples)
    object$S_samples <- object$S_samples+Z_sim
  }

  if(is.matrix(mu_target)) {
    out$lp_samples <- sapply(1:n_samples,
                             function(i)
                               mu_target[,i] + cov_offset +
                               object$S_samples[,i])
  } else {
    out$lp_samples <- sapply(1:n_samples,
                             function(i)
                               mu_target + cov_offset +
                               object$S_samples[,i])
  }


  if(include_re) {
    n_dim_re <- sapply(1:n_re, function(i) length(object$re$samples[[i]]))

    for(i in 1:n_re) {
      for(j in 1:n_dim_re[i]) {
        for(h in 1:n_samples) {
          out$lp_samples[,h] <- out$lp_samples[,h] +
          object$re$D_pred[[i]][,j]*object$re$samples[[i]][[j]][h]
        }
      }
    }
  }

  names_f <- names(f_target)
  if(is.null(names_f)) names_f <- paste("f_target_",1:length(f_target), sep = "")

  names_s <- names(pd_summary)
  if(is.null(pd_summary)) names_s <- paste("pd_summary_",1:length(f_target), sep = "")

  out$target <- list()

  for(i in 1:n_f) {
    target_samples_i <-
      f_target[[i]](out$lp_samples)
    if(dast_model && include_mda_effect) {
      alpha <- object$par_hat$alpha
      if(is.null(alpha)) alpha <- object$fix_alpha
      gamma <- object$par_hat$gamma
      mda_effect_time_pred <- compute_mda_effect(rep(time_pred, n_pred),
                                                 mda_times = object$mda_times,
                              mda_grid, alpha = alpha,
                              gamma = gamma, kappa = object$power_val)
      target_samples_i <- target_samples_i*mda_effect_time_pred

    }
    out$target[[paste(names_f[i])]] <- list()
    for(j in 1:n_summaries) {
      out$target[[paste(names_f[i])]][[paste(names_s[j])]] <-
        apply(target_samples_i, 1, pd_summary[[j]])
    }
  }
  out$grid_pred <- object$grid_pred
  out$f_target <- names(f_target)
  out$pd_summary <- names(pd_summary)
  class(out) <- "RiskMap_pred_target_grid"
  return(out)
}

##' Plot Method for RiskMap_pred_target_grid Objects
##'
##' Generates a plot of the predicted values or summaries over the regular spatial grid
##' from an object of class 'RiskMap_pred_target_grid'.
##'
##' @param x An object of class 'RiskMap_pred_target_grid'.
##' @param which_target Character string specifying which target prediction to plot.
##' @param which_summary Character string specifying which summary statistic to plot (e.g., "mean", "sd").
##' @param ... Additional arguments passed to the \code{\link[terra]{plot}} function of the \code{terra} package.
##' @return A \code{ggplot} object representing the specified prediction target or summary statistic over the spatial grid.
##' @details
##' This function requires the 'terra' package for spatial data manipulation and plotting.
##' It plots the values or summaries over a regular spatial grid, allowing for visual examination of spatial patterns.
##'
##' @seealso \code{\link{pred_target_grid}}
##'
##' @importFrom terra as.data.frame rast plot
##' @method plot RiskMap_pred_target_grid
##' @export
##'
##'
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
plot.RiskMap_pred_target_grid <- function(x, which_target = "linear_target", which_summary = "mean", ...) {
  t_data.frame <-
    terra::as.data.frame(cbind(st_coordinates(x$grid_pred),
                               x$target[[which_target]][[which_summary]]),
                         xy = TRUE)
  raster_out <- terra::rast(t_data.frame, crs = st_crs(x$grid_pred)$input)

  terra::plot(raster_out, ...)
}


##' Predictive Target over a Shapefile
##'
##' Computes predictions over a shapefile using outputs from the
##' \code{\link{pred_over_grid}} function.
##' This function allows for incorporating covariates, offsets, and optional
##' unstructured random effects into the predictive target.
##'
##' @param object Output from `pred_over_grid`, a RiskMap.pred.re object.
##' @param shp Spatial dataset (sf or data.frame) representing the shapefile over which predictions are computed.
##' @param shp_target Function defining the aggregation method for shapefile targets (default is mean).
##' @param weights Optional numeric vector of weights for spatial predictions.
##' @param standardize_weights Logical indicating whether to standardize weights (default is FALSE).
##' @param col_names Column name or index in 'shp' containing region names.
##' @param include_covariates Logical indicating whether to include covariates in predictions (default is TRUE).
##' @param include_nugget Logical indicating whether to include the nugget effect (default is FALSE).
##' @param include_cov_offset Logical indicating whether to include covariate offset in predictions (default is FALSE).
##' @param include_mda_effect Logical indicating whether to include the mass drug administration (MDA) effect as defined by the fitted DAST model (default is TRUE).
##' @param return_shp Logical indicating whether to return the shape file with the added predictive distribution summaries as defined through \code{pd_summary}.
##' @param include_re Logical indicating whether to include random effects in predictions (default is FALSE).
##' @param f_target List of target functions to apply to the linear predictor samples.
##' @param pd_summary List of summary functions (e.g., mean, sd) to summarize target samples.
##'
##' @importFrom terra rast as.data.frame
##' @export
##' @details
##' This function computes predictive targets or summaries over a spatial shapefile
##' using outputs from 'pred_S'. It requires the 'terra' package for spatial data
##' manipulation and should be used with 'sf' or 'data.frame' objects representing
##' the shapefile.
##'
##' @return An object of class 'RiskMap_pred_target_shp' containing computed targets,
##' summaries, and associated spatial data.
##'
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##'
##' @seealso
##' \code{\link{pred_target_grid}}
##'
pred_target_shp <- function(object, shp, shp_target=mean,
                            weights = NULL,standardize_weights = FALSE,
                            col_names=NULL,
                            include_covariates = TRUE,
                            include_nugget = FALSE,
                            include_cov_offset = FALSE,
                            include_mda_effect=TRUE,
                            return_shp=TRUE,
                            time_pred = NULL,
                            mda_grid = NULL,
                            include_re = FALSE,
                            f_target = NULL,
                            pd_summary = NULL,
                            messages = TRUE) {
  if(!inherits(object,
               what = "RiskMap.pred.re", which = FALSE)) {
    stop("The object passed to 'object' must be an output of
         the function 'pred_S'")
  }

  if(object$type!="joint") {
    stop("To run predictions with a shape file, joint predictions must be used;
         rerun 'pred_over_grid' and set type='joint'")
  }

  dast_model <- !is.null(object$par_hat$gamma)

  if(dast_model) {
    if(include_mda_effect & is.null(mda_grid)) {
      stop("The MDA coverage must be specified for each point on the grid through the argument 'mda_grid'")
    }
    if(is.null(time_pred)) {
      stop("For a DAST model, the time of prediction must be specified through the argument 'time_pred'")
    }
  }


  n_pred <- nrow(object$S_samples)

  n_re <- length(object$re$samples)
  re_names <- names(object$re$samples)

  if(n_re==0 &&
     include_re) {
    stop("The categories of the randome effects variables have not been provided;
         re-run pred_over_grid and provide the covariates through the argument 're_predictors'")
  }

  if(!is.null(weights)) {
    if(!inherits(weights,
                 what = c("numeric"), which = FALSE)) {
      stop("The object passed to 'weights' must be a numeric
         vector of weights with length equal to the number of
           prediction locations")
    }
    if(length(weights)!=n_pred) stop("The value passed to 'weights' do not match
                                   the prediction grid")
    no_weights <- FALSE
  } else {
    weights <- rep(1,n_pred)
    no_weights <- TRUE
  }

  if(is.null(object$par_hat$tau2) &
     include_nugget) {
    stop("The nugget cannot be included in the predictive target
             because it was not included when fitting the model")
  }

  if(is.null(f_target)) {
    f_target <- list(linear_target = function(x) x)
  }

  if(is.null(pd_summary)) {
    pd_summary <- list(mean = mean, sd = sd)
  }

  n_summaries <- length(pd_summary)
  n_f <- length(f_target)
  n_samples <- ncol(object$S_samples)

  if(length(object$mu_pred)==1 && object$mu_pred==0 &&
     include_covariates) {
    stop("Covariates have not been provided; re-run pred_over_grid
         and provide the covariates through the argument 'predictors'")
  }
  if(!include_covariates) {
    mu_target <- 0
  } else {
    if(is.null(object$mu_pred)) stop("the output obtained from 'pred_S' does not
                                     contain any covariates; if including covariates
                                     in the predictive target these shuold be included
                                     when running 'pred_S'")
    mu_target <- object$mu_pred
  }

  if(!include_cov_offset) {
    cov_offset <- 0
  } else {
    if(length(object$cov_offset)==1) {
      stop("No covariate offset was included in the model;
           set include_cov_offset = FALSE, or refit the model and include
           the covariate offset")
    }
    cov_offset <- object$cov_offset
  }

  if(include_nugget) {
    Z_sim <- matrix(rnorm(n_samples*n_pred,
                          sd = sqrt(object$par_hat$tau2)),
                    ncol = n_samples)
    object$S_samples <- object$S_samples+Z_sim
  }

  out <- list()
  if(is.matrix(mu_target)) {
    out$lp_samples <- sapply(1:n_samples,
                             function(i)
                               mu_target[,i] + cov_offset +
                               object$S_samples[,i])
  } else {
    out$lp_samples <- sapply(1:n_samples,
                             function(i)
                               mu_target + cov_offset +
                               object$S_samples[,i])
  }
  if(include_re) {
    n_dim_re <- sapply(1:n_re, function(i) length(object$re$samples[[i]]))

    for(i in 1:n_re) {
      for(j in 1:n_dim_re[i]) {
        for(h in 1:n_samples) {
          out$lp_samples[,h] <- out$lp_samples[,h] +
            object$re$D_pred[[i]][,j]*object$re$samples[[i]][[j]][h]
        }
      }
    }
  }

  names_f <- names(f_target)
  names_s <- names(pd_summary)
  out$target <- list()

  n_reg <- nrow(shp)
  if(is.null(col_names)) {
    shp$region <- paste("reg",1:n_reg, sep="")
    col_names <- "region"
    names_reg <- shp$region
  } else {
    names_reg <- shp[[col_names]]
    if(n_reg != length(names_reg)) {
      stop("The names in the column identified by 'col_names' do not
         provide a unique set of names, but there are duplicates")
    }
  }

  shp <- st_transform(shp, crs = st_crs(object$grid_pred)$input)
  inter <- st_intersects(shp, object$grid_pred)

  if(any(is.na(weights))) {
    warning("Missing values found in 'weights' are set to 0 \n")
    weights[is.na(weights)] <- 0
  }
  no_comp <- NULL
  for(h in 1:n_reg) {
    if(messages) message("Computing predictive target for:",shp[[col_names]][h])
    if(length(inter[[h]])==0) {
      warning(paste("No points on the grid fall within", shp[[col_names]][h],
                    "and no predictions are carried out for this area"))
      no_comp <- c(no_comp, h)
    } else {
      ind_grid_h <- inter[[h]]
      if(standardize_weights & !no_weights) {
        weights_h <- weights[ind_grid_h]/sum(weights[ind_grid_h])
      } else {
        weights_h <- weights[ind_grid_h]
      }
      for(i in 1:n_f) {
        target_grid_samples_i <- as.matrix(f_target[[i]](out$lp_samples[ind_grid_h,]))
        if(dast_model && include_mda_effect) {
          alpha <- object$par_hat$alpha
          if(is.null(alpha)) alpha <- object$fix_alpha
          gamma <- object$par_hat$gamma
          mda_effect_time_pred <- compute_mda_effect(rep(time_pred, length(ind_grid_h)),
                                                     mda_times = object$mda_times,
                                                     mda_grid[ind_grid_h,], alpha = alpha,
                                                     gamma = gamma, kappa = object$power_val)
          target_grid_samples_i <- target_grid_samples_i*mda_effect_time_pred

        }

        target_samples_i <-
          apply(target_grid_samples_i,
                2, function(x) shp_target(weights_h*x))

        out$target[[paste(names_reg[h])]][[paste(names_f[i])]] <- list()
        for(j in 1:n_summaries) {
          out$target[[paste(names_reg[h])]][[paste(names_f[i])]][[paste(names_s[j])]] <-
            pd_summary[[j]](target_samples_i)
        }
      }
    }
    if(messages) message(" \n")
  }

  if(return_shp) {
    if(length(no_comp) > 0) {
      ind_reg <- (1:n_reg)[-no_comp]
    } else {
      ind_reg <- 1:n_reg
    }
    for(i in 1:n_f) {
      for(j in 1:n_summaries) {
        name_ij <- paste(names_f[i],"_",paste(names_s[j]),sep="")
        shp[[name_ij]] <- rep(NA, n_reg)
        for(h in ind_reg) {
          which_reg <- which(shp[[col_names]]==names_reg[h])
          shp[which_reg,][[name_ij]] <-
            out$target[[paste(names_reg[h])]][[paste(names_f[i])]][[paste(names_s[j])]]
        }
      }
    }
  }
  out$shp <- shp
  out$f_target <- names(f_target)
  out$pd_summary <- names(pd_summary)
  out$grid_pred <- object$grid_pred
  class(out) <- "RiskMap_pred_target_shp"
  return(out)
}


##' Plot Method for RiskMap_pred_target_shp Objects
##'
##' Generates a plot of predictive target values or summaries over a shapefile.
##'
##' @param x An object of class 'RiskMap_pred_target_shp' containing computed targets,
##' summaries, and associated spatial data.
##' @param which_target Character indicating the target type to plot (e.g., "linear_target").
##' @param which_summary Character indicating the summary type to plot (e.g., "mean", "sd").
##' @param ... Additional arguments passed to 'scale_fill_distiller' in 'ggplot2'.
##' @return A \code{ggplot} object showing the plot of the specified predictive target or summary.
##' @details
##' This function plots the predictive target values or summaries over a shapefile.
##' It requires the 'ggplot2' package for plotting and 'sf' objects for spatial data.
##'
##' @seealso
##' \code{\link{pred_target_shp}}, \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{geom_sf}},
##' \code{\link[ggplot2]{aes}}, \code{\link[ggplot2]{scale_fill_distiller}}
##'
##' @importFrom ggplot2 ggplot geom_sf aes scale_fill_distiller
##' @method plot RiskMap_pred_target_shp
##' @export
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
plot.RiskMap_pred_target_shp <- function(x, which_target = "linear_target",
                                         which_summary = "mean", ...) {
  col_shp_name <- paste(which_target,"_",which_summary,sep="")

  out <- ggplot(x$shp) +
    geom_sf(aes(fill = x$shp[[col_shp_name]])) +
    scale_fill_distiller(...)
  return(out)
}

##' @title Update Predictors for a RiskMap Prediction Object
##'
##' @description
##' This function updates the predictors of a given RiskMap prediction object. It ensures that the new predictors match the original prediction grid and updates the relevant components of the object accordingly.
##'
##' @param object A `RiskMap.pred.re` object, which is the output of the \code{\link{pred_over_grid}} function.
##' @param predictors A data frame containing the new predictor values. The number of rows must match the prediction grid in the `object`.
##'
##' @details
##' The function performs several checks and updates:
##' \itemize{
##'   \item Ensures that `object` is of class `RiskMap.pred.re`.
##'   \item Ensures that the number of rows in `predictors` matches the prediction grid in `object`.
##'   \item Removes any rows with missing values in `predictors` and updates the corresponding components of the `object`.
##'   \item Updates the prediction locations, the predictive samples for the random effects, and the linear predictor.
##' }
##'
##' @return The updated `RiskMap.pred.re` object.
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##' @export
update_predictors <- function(object,
                              predictors) {
  if(!inherits(object,
               what = "RiskMap.pred.re", which = FALSE)) {
    stop("The object passed to 'object' must be an output of
         the function 'glgpm'")
  }
  grp <- st_coordinates(object$grid_pred)
  n_pred <- nrow(grp)
  par_hat <- object$par_hat
  p <- length(par_hat$beta)

  inter_f <- object$inter_f
  inter_lt_f <- inter_f
  inter_lt_f$pf <- update(inter_lt_f$pf, NULL ~.)

  if(p==1) {
    stop("No update of the prectors can be done for an intercept only model")
  }

  if(!is.data.frame(predictors)) stop("'predictors' must be an object of class 'data.frame'")
  if(nrow(predictors)!=n_pred) stop("the values provided for 'predictors' do not match the prediction grid passed to 'grid_pred'")

  if(any(is.na(predictors))) {
    warning("There are missing values in 'predictors'; these values have been removed
              alongside the corresponding prediction locations, and predictive samples
            for the random effects")
  }

  if(!is.null(object$re_predictors)) {

    comb_pred <- data.frame(predictors,
                            object$re_predictors)
    ind_c <- complete.cases(comb_pred)
    predictors_aux <- data.frame(na.omit(comb_pred)[,1:ncol(predictors)])
    colnames(predictors_aux) <- colnames(predictors)
    predictors <- predictors_aux

    re_predictors_aux <- data.frame(na.omit(comb_pred)[,(ncol(predictors)+1):
                                                         (ncol(predictors)+
                                                            ncol(object$re_predictors))])
    colnames(re_predictors_aux) <- colnames(re_predictors_aux)
    object$re_predictors <- re_predictors_aux

    n_re <- length(object$re$samples)
    for(i in 1:n_re) {
      object$re$D_pred[[i]] <- object$re$D_pred[[i]][ind_c,]
    }

  } else {
    comb_pred <- predictors
    ind_c <- complete.cases(comb_pred)
    predictors_aux <- data.frame(na.omit(comb_pred))
    colnames(predictors_aux) <- colnames(predictors)
    predictors <- predictors_aux
  }

  grid_pred_aux <- st_coordinates(object$grid_pred)[ind_c,]
  object$grid_pred <- st_as_sf(data.frame(grid_pred_aux),
                               coords = c("X","Y"),
                               crs = st_crs(object$grid_pred)$input)
  grp <- st_coordinates(object$grid_pred)
  n_pred <- nrow(grp)
  object$S_samples <- object$S_samples[ind_c,]

  mf_pred <- model.frame(inter_lt_f$pf,data=predictors, na.action = na.fail)
  D_pred <- as.matrix(model.matrix(attr(mf_pred,"terms"),data=predictors))
  if(ncol(D_pred)!=p) stop("the provided variables in 'predictors' do not match the number of explanatory variables used to fit the model.")
  mu_pred <- as.numeric(D_pred%*%par_hat$beta)

  object$mu_pred <- mu_pred
  out <- object
  class(out) <- "RiskMap.pred.re"
  return(out)
}

##' @title Assess Predictive Performance via Spatial Cross-Validation
##'
##' @description
##' This function evaluates the predictive performance of spatial models fitted to `RiskMap` objects using cross-validation. It supports two classes of diagnostic tools:
##'
##' - **Scoring rules**, including the Continuous Ranked Probability Score (CRPS) and its scaled version (SCRPS), which quantify the sharpness and calibration of probabilistic forecasts;
##' - **Calibration diagnostics**, based on the Probability Integral Transform (PIT) for Gaussian outcomes and Aggregated nonparametric PIT (AnPIT) curves for discrete outcomes (e.g., Poisson or Binomial).
##'
##' Cross-validation can be performed using either spatial clustering or regularized subsampling with a minimum inter-point distance. For each fold or subset, models can be refitted or evaluated with fixed parameters, offering flexibility in model validation. The function also provides visualizations of the spatial distribution of test folds.
##'
##' @param object A list of `RiskMap` objects, each representing a model fitted with `glgpm`.
##' @param keep_par_fixed Logical; if `TRUE`, parameters are kept fixed across folds, otherwise the model is re-estimated for each fold.
##' @param iter Integer; number of times to repeat the cross-validation.
##' @param fold Integer; number of folds for cross-validation (required if `method = "cluster"`).
##' @param n_size Optional; the size of the test set, required if `method = "regularized"`.
##' @param control_sim Control settings for simulation, an output from `set_control_sim`.
##' @param method Character; either `"cluster"` or `"regularized"` for the cross-validation method. The `"cluster"` method uses
##' spatial clustering as implemented by the \code{\link{spatial_clustering_cv}} function from the `spatialEco` package, while the `"regularized"` method
##' selects a subsample of the dataset by imposing a minimum distance, set by the `min_dist` argument, for a randomly selected
##' subset of locations.
##' @param min_dist Optional; minimum distance for regularized subsampling (required if `method = "regularized"`).
##' @param plot_fold Logical; if `TRUE`, plots each fold's test set.
##' @param messages Logical; if `TRUE`, displays progress messages.
##' @param which_metric Character vector; one or more of `"CRPS"`, `"SCRPS"`, or `"AnPIT"`, to specify the predictive performance metrics to compute.
##' @param ... Additional arguments passed to clustering or subsampling functions.
##'
##' @return A list of class `RiskMap.spatial.cv`, containing:
##' \describe{
##'   \item{test_set}{A list of test sets used for validation, each of class `'sf'`.}
##'   \item{model}{A named list, one per model, each containing:
##'     \describe{
##'       \item{score}{A list with CRPS and/or SCRPS scores for each fold if requested.}
##'       \item{PIT}{(if `family = "gaussian"` and `which_metric` includes `"AnPIT"`) A list of PIT values for test data.}
##'       \item{AnPIT}{(if `family` is discrete and `which_metric` includes `"AnPIT"`) A list of AnPIT curves for test data.}
##'     }
##'   }
##' }
##'
##' @seealso \code{\link{spatial_clustering_cv}}, \code{\link{subsample.distance}}, \code{\link{plot_AnPIT}}
##'
##' @references
##' Bolin, D., & Wallin, J. (2023). Local scale invariance and robustness of proper scoring rules. *Statistical Science*, 38(1), 140–159. \doi{10.1214/22-STS864}.
##'
##' @importFrom terra match
##' @importFrom ggplot2 ggplot geom_sf theme_minimal ggtitle
##' @importFrom gridExtra grid.arrange
##' @importFrom stats ecdf integrate rbinom rpois
##' @importFrom spatialEco subsample.distance
##' @importFrom spatialsample spatial_clustering_cv autoplot
##' @importFrom sf st_as_sfc
##' @export
##' @author Emanuele Giorgi
assess_pp <- function(object,
                      keep_par_fixed = TRUE,
                      iter = 1,
                      fold = NULL, n_size = NULL,
                      control_sim = set_control_sim(),
                      method, min_dist = NULL,
                      plot_fold = TRUE,
                      messages = TRUE,
                      which_metric = c("AnPIT", "CRPS", "SCRPS"),
                      ...) {

  ## ───────────────────────────────────── helpers ────────────────────────────── ##
  is_list_of_riskmap <- function(object) {
    is.list(object) &&
      all(vapply(object, inherits, logical(1), what = "RiskMap"))
  }

  crps_gaussian <- function(y, mu, sigma) {
    if (sigma == 0) return(0)
    z   <- (y - mu) / sigma
    2*dnorm(z) + z * (2 * pnorm(z) - 1) - 1 / sqrt(pi)
  }

  ##  Discrete CRPS  (Binomial / Poisson)  ------------------------------------ ##
  crps_discrete <- function(y, Fk) {
    k      <- seq_along(Fk) - 1
    sum((Fk - as.numeric(k >= y))^2)
  }

  ## expectation of CRPS under the predictive pmf  (needed for SCRPS) ---------- ##
  exp_crps_discrete <- function(Fk, pk) {
    k <- seq_along(pk) - 1
    sum(pk * vapply(k, crps_discrete, numeric(1), Fk = Fk))
  }

  ## ───────────────────────────── sanity checks (unchanged) ─────────────────── ##
  if (!is_list_of_riskmap(object))
    stop("`object` must be a list of model fits obtained with `glgpm()`")

  if (!all(which_metric %in% c("CRPS", "SCRPS", "AnPIT")))
    stop("`which_metric` must only contain 'CRPS', 'SCRPS' or 'AnPIT'")

  if (!method %in% c("cluster", "regularized"))
    stop("`method` must be either 'cluster' or 'regularized'")

  if (method == "regularized") {
    if (is.null(min_dist))
      stop("for 'regularized', supply `min_dist`")
    if (is.null(n_size))
      stop("for 'regularized', supply `n_size`")
  }

  if (method == "cluster" && is.null(fold))
    stop("for 'cluster', supply `fold`")

  if (!inherits(control_sim, "mcmc.RiskMap"))
    stop("`control_sim` must come from `set_control_sim()`")

  get_CRPS  <- "CRPS"  %in% which_metric
  get_SCRPS <- "SCRPS" %in% which_metric
  get_AnPIT <- "AnPIT" %in% which_metric      # will hold PIT for Gaussian

  ## ─────────────────────────────── test–train splitting ────────────────────── ##
  object1 <- object[[1]]
  data_sf <- object1$data_sf

  if (method == "cluster") {
    data_split <- spatial_clustering_cv(data = data_sf,
                                        v = fold, repeats = iter, ...)
  } else {                                            # regularised distance split
    data_split <- list(splits = vector("list", iter))
    for (i in seq_len(iter)) {
      data_split$splits[[i]]            <- list()
      ## --------- sample at the *location* level (duplicates removed) -------- ##
      locations_sf                       <- data_sf[!duplicated(st_as_text(data_sf$geometry)), ]
      data_split$splits[[i]]$data_test   <- subsample.distance(
        locations_sf, size = n_size, d = min_dist * 1000, ...
      )
      test_geom                          <- st_as_text(data_split$splits[[i]]$data_test$geometry)
      in_test                            <- st_as_text(data_sf$geometry) %in% test_geom
      data_split$splits[[i]]$out_id      <- which(in_test)
      data_split$splits[[i]]$in_id       <- which(!in_test)
      data_split$splits[[i]]$data        <- data_sf[!in_test, ]
    }
  }

  n_iter <- if (method == "cluster") iter * fold else iter
  if (plot_fold) {
    if (method == "cluster") {
      print(autoplot(data_split))
    } else {
      library(gridExtra)
      plots <- lapply(seq_len(n_iter), function(i)
        ggplot(data_split$splits[[i]]$data_test) +
          geom_sf() + theme_minimal() +
          ggtitle(paste("Subset", i))
      )
      if (n_iter > 1) grid.arrange(grobs = plots, ncol = 2) else print(plots[[1]])
    }
  }

  ## ───────────────────────── initialise output containers ──────────────────── ##
  n_models   <- length(object)
  model_names <- names(object)
  out <- list(test_set = vector("list", n_iter), model = list())

  ## ─────────────────────────── iterate over models ─────────────────────────── ##
  for (h in seq_len(n_models)) {

    par_hat   <- coef(object[[h]])
    den_name  <- as.character(object[[h]]$call$den)
    fam       <- object[[h]]$family             # "gaussian", "binomial", "poisson"
    linkfun   <- switch(fam,
                        gaussian = identity,
                        binomial = function(x) plogis(x),
                        poisson  = function(x) exp(x))

    ## containers for this model --------------------------------------------- ##
    if (get_CRPS)   CRPS  <- vector("list", n_iter)
    if (get_SCRPS) { y_CRPS <- vector("list", n_iter)
    SCRPS  <- vector("list", n_iter) }
    if (get_AnPIT) {           # will store PIT if Gaussian
      if (fam == "gaussian")  PIT  <- vector("list", n_iter)
      else                    AnPIT <- vector("list", n_iter)
    }

    ## ───────────── iterate over CV splits (refit where requested) ─────────── ##
    for (i in seq_len(n_iter)) {

      ## ------------- optionally refit -------------------------------------- ##
      if (!keep_par_fixed) {
        message("\nRe-estimating model for subset ", i)
        refit_i <- eval(bquote(
          glgpm(.(
            formula      = object[[h]]$formula,
            data         = data_sf[data_split$splits[[i]]$in_id, ],
            cov_offset   = .(object[[h]]$cov_offset),
            family       = .(fam),
            crs          = .(object[[h]]$crs),
            scale_to_km  = .(object[[h]]$scale_to_km),
            control_mcmc = control_sim,
            fix_var_me   = .(object[[h]]$fix_var_me),
            den          = .(as.name(den_name)),
            messages     = FALSE,
            start_pars   = par_hat
          ))
        ))
      } else {                                  # quick copy without re-fitting
        refit_i <- object[[h]]
        keep    <- data_split$splits[[i]]$in_id
        drop    <- data_split$splits[[i]]$out_id
        refit_i$data_sf  <- refit_i$data_sf [keep, ]
        refit_i$units_m  <- refit_i$units_m[keep]

        keep_coord <- unique(refit_i$ID_coords[keep])   # rows in coords to keep
        refit_i$coords <- refit_i$coords[keep_coord, , drop = FALSE]

        refit_i$y        <- refit_i$y      [keep]
        refit_i$D        <- refit_i$D      [keep, ]
        if (!is.null(refit_i$cov_offset))
          refit_i$cov_offset <- refit_i$cov_offset[keep]
        refit_i$ID_coords <- compute_ID_coords(refit_i$data_sf)$ID_coords
      }

      ## ------------- posterior predictive samples for *held-out* points ---- ##
      data_test_i <- data_sf[data_split$splits[[i]]$out_id, ]
      data_test_i <- data_test_i[complete.cases(st_drop_geometry(data_test_i)), ]

      out$test_set[[i]] <- data_test_i
      pred_coff_i <- if(is.null(object[[h]]$cov_offset)) {
        rep(0, nrow(data_test_i))
      } else {
        object[[h]]$cov_offset[data_split$splits[[i]]$out_id]
      }

      message("\nModel: ", model_names[h],
              "\nSpatial prediction for subset ", i, "\n")

      pred_S <- pred_over_grid(object       = refit_i,
                               grid_pred     = st_as_sfc(data_test_i),
                               control_sim   = control_sim,
                               predictors    = data_test_i,
                               pred_cov_offset = pred_coff_i,
                               type          = "marginal",
                               messages      = FALSE)

      pred_lp  <- pred_target_grid(pred_S,
                                   include_nugget    = is.null(refit_i$fix_tau2) ||
                                     refit_i$fix_tau2 != 0,
                                   include_cov_offset = !all(refit_i$cov_offset==0))

      eta_samp <- pred_lp$lp_samples                  # linear predictor samples
      if(is.null(refit_i$fix_var_me)) {
        sigma2_me <- coef(refit_i)$sigma2_me
      } else {
        sigma2_me <- refit_i$fix_var_me
      }
      eta_samp <- eta_samp + sqrt(sigma2_me)*rnorm(length(eta_samp))


      mu_samp  <- linkfun(eta_samp)                   # on response scale

      n_pred   <- nrow(eta_samp)
      n_draw   <- ncol(eta_samp)

      if (get_CRPS)   CRPS [[i]] <- numeric(n_pred)
      if (get_SCRPS) { y_CRPS[[i]] <- numeric(n_pred)
      SCRPS [[i]] <- numeric(n_pred) }

      if (get_AnPIT) {
        if (fam == "gaussian")  {
          PIT_i  <- numeric(n_pred)
        } else {
          u_val      <- seq(0, 1, length.out = 1000)
          AnPIT_i    <- matrix(NA_real_, nrow = length(u_val), ncol = n_pred)
          npit_fun   <- function(y, u, Fk) {
            F_y1 <- if (y == 0) 0 else Fk[y]
            F_y  <- Fk[y + 1]
            ifelse(u <= F_y1, 0,
                   ifelse(u <= F_y, (u - F_y1) / (F_y - F_y1), 1))
          }
        }
      }

      ## --------------------------- loop over locations --------------------- ##
      units_m_i <- object[[h]]$units_m[data_split$splits[[i]]$out_id]
      y_i       <- object[[h]]$y      [data_split$splits[[i]]$out_id]

      for (j in seq_len(n_pred)) {
        ## ---------- family-specific predictive pmf/cdf --------------------- ##
        if (fam == "gaussian") {

          mu_j <- mean(mu_samp[j, ])
          sd_j <- sd  (mu_samp[j, ])

          ## CRPS closed form
          if (get_CRPS)
            CRPS[[i]][j] <- sd_j * crps_gaussian(y_i[j], mu_j, sd_j)

          ## SCRPS (scaled)
          if (get_SCRPS) {
            y_CRPS[[i]][j] <- sd_j / sqrt(pi)
            SCRPS [[i]][j] <- -0.5 * (1 + CRPS[[i]][j] / y_CRPS[[i]][j] +
                                        log(2 * abs(y_CRPS[[i]][j])))
          }

          ## PIT
          if (get_AnPIT) {
            PIT_i[j] <- pnorm(y_i[j], mean = mu_j, sd = sd_j)
          }
        } else {                                # ----- Binomial / Poisson ---- #
          if (fam == "binomial") {
            y_samp <- rbinom(n_draw, size = units_m_i[j], prob = mu_samp[j, ])
            support <- 0:units_m_i[j]
          } else {                              # Poisson
            lambda  <- units_m_i[j] * mu_samp[j, ]
            y_samp <- rpois(n_draw, lambda)
            support <- 0:max(max(y_samp), y_i[j], qpois(0.999, mean(lambda)))
          }

          ## empirical pmf / cdf from the draws
          pk           <- tabulate(y_samp + 1, nbins = length(support)) / n_draw
          Fk           <- cumsum(pk)

          ## CRPS discrete
          if (get_CRPS)
            CRPS[[i]][j] <- crps_discrete(y_i[j], Fk)

          ## SCRPS
          if (get_SCRPS) {
            y_CRPS[[i]][j] <- exp_crps_discrete(Fk, pk)
            SCRPS [[i]][j] <- -0.5 * (1 + CRPS[[i]][j] / y_CRPS[[i]][j] +
                                        log(2 * abs(y_CRPS[[i]][j])))
          }

          ## AnPIT (aggregate NPIT over u grid)
          if (get_AnPIT) {
            npit_u <- vapply(u_val,
                             npit_fun, numeric(1),
                             y = y_i[j], Fk = Fk)
            AnPIT_i[, j] <- npit_u
          }
        } # end family branch
      }   # end j loop

      ## store metrics for this subset --------------------------------------- ##
      if (get_AnPIT) {
        if (fam == "gaussian") {
          PIT[[i]]   <- PIT_i
        } else  {
          AnPIT[[i]] <- rowMeans(AnPIT_i)
        } # end j loop

        ## store metrics for this subset --------------------------------------- ##
        if (get_AnPIT) {
          if (fam == "gaussian") {
            PIT[[i]]   <- PIT_i
          } else  {
            AnPIT[[i]] <- rowMeans(AnPIT_i)
          }
        }
      }   # end i (fold) loop

      ## ─────────────── finalise output for this model ─────────────────────── ##
      out$model[[model_names[h]]] <- list(score = list())
      if (get_CRPS)  out$model[[model_names[h]]]$score$CRPS  <- CRPS
      if (get_SCRPS) out$model[[model_names[h]]]$score$SCRPS <- SCRPS

      if (get_AnPIT) {
        if (fam == "gaussian") {
          out$model[[model_names[h]]]$PIT   <- PIT
        } else {
          out$model[[model_names[h]]]$AnPIT <- AnPIT
        }
      }
    }
  }

  class(out) <- "RiskMap.spatial.cv"
  return(out)
}


##' Simulate surface data based on a spatial model
##'
##' This function simulates surface data based on a user-defined formula and other parameters. It allows for simulation of spatial data with various model families (Gaussian, Binomial, or Poisson). The simulation involves creating spatially correlated random fields and generating outcomes for data points in a given prediction grid.
##'
##' @param n_sim The number of simulations to run.
##' @param pred_grid A spatial object (either `sf` or `data.frame`) representing the prediction grid where the simulation will take place.
##' @param formula A formula object specifying the model to be fitted. It should include both fixed effects and random effects if applicable.
##' @param sampling_f A function that returns a sampled dataset (of class `sf` or `data.frame`) to simulate data from.
##' @param family A character string specifying the family of the model. Must be one of "gaussian", "binomial", or "poisson".
##' @param scale_to_km A logical indicating whether the coordinates should be scaled to kilometers. Defaults to `TRUE`.
##' @param control_mcmc A list of control parameters for MCMC (not used in this implementation but can be expanded later).
##' @param par0 A list containing initial parameter values for the simulation, including `beta`, `sigma2`, `phi`, `tau2`, and `sigma2_me`.
##' @param include_covariates A logical indicateing if the covariates (or the intercept if no covariates are used) should be included in the linear
##' predictor. By default \code{include_covariates = TRUE}
##' @param nugget_over_grid A logical indicating whether to include a nugget effect over the entire prediction grid.
##' @param fix_var_me A parameter to fix the variance of the random effects for the measurement error. Defaults to `NULL`.
##' @param messages A logical value indicating whether to print messages during the simulation. Defaults to `TRUE`.
##'
##' @return A list containing the simulated data (\code{data_sim}), the linear predictors (\code{lp_grid_sim}),
##' a logical value indicating if covariates have been included in the linear predictor (\code{include_covariates}),
##' a logical value indicating if the nugget has been included into the simulations of the linear predictor over the grid
##' (\code{nugget_over_grid}), a logical  indicating if a covariate offset has been included in the linear predictor (\code{include_cov_offset}),
##' the model parameters set for the simulation (\code{par0}) and the family used in the model (\code{family}).
##'
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @export
surf_sim <- function(n_sim, pred_grid,
                     formula, sampling_f,
                     family,
                     scale_to_km = TRUE,
                     control_mcmc = set_control_sim(),
                     par0, nugget_over_grid = FALSE,
                     include_covariates = TRUE,
                     fix_var_me = NULL,
                     messages = TRUE) {


  if(!inherits(formula,
               what = "formula", which = FALSE)) {
    stop("'formula' must be a 'formula'
                                     object indicating the variables of the
                                     model to be fitted")
  }
  inter_f <- interpret.formula(formula)
  include_cov_offset <- !is.null(inter_f$offset)
  if(!inherits(pred_grid,
               what = c("sf", "data.frame"), which = FALSE)) {
    stop("'pred_grid' must be an 'sf'
          object indicating the variables of the
          model to be fitted")
  }

  if(!inherits(sampling_f,
               what = c("function"), which = FALSE)) {
    stop("'sampling_f' must be an object of class 'function'")
  }

  data_test <- sampling_f()
  if(!inherits(data_test,
               what = c("sf", "data.frame"), which = FALSE)) {
    stop("The object return by 'sampling_f' must be of an 'sf' object")
  }

  if (!"units_m" %in% colnames(data_test)) {
    stop("The object returned by 'sampling_f' must contain a column named 'units_m'.")
  }


  data_sim <- list()
  coords_sim <- list()
  for(i in 1:n_sim) {
    data_sim[[i]] <- sampling_f()
    coords_sim[[i]] <- st_coordinates(data_sim[[i]])
    if(scale_to_km) coords_sim[[i]] <- coords_sim[[i]]/1000
    if (st_crs(data_sim[[i]]) != st_crs(pred_grid)) {
      pred_grid <- st_transform(pred_grid, st_crs(data_sim[[i]]))
    }

    # Find nearest neighbor in 'pred_grid' for each feature in 'data'
    nearest_indices <- st_nearest_feature(data_sim[[i]], pred_grid)

    # Extract variables from the nearest features in 'pred_grid'
    pred_grid_vars <- st_drop_geometry(pred_grid[nearest_indices, ])

    # Bind the extracted variables to 'data'
    data_sim[[i]] <- cbind(data_sim[[i]], pred_grid_vars)
  }


  kappa <- inter_f$gp.spec$kappa
  if(kappa < 0) stop("kappa must be positive.")

  if(family != "gaussian" & family != "binomial" &
     family != "poisson") stop("'family' must be either 'gaussian', 'binomial'
                               or 'poisson'")

  inter_lt_f <- inter_f
  inter_lt_f$pf <- update(inter_lt_f$pf, NULL ~.)
  D <- list()
  n <- list()
  for(i in 1:n_sim) {
    mf <- model.frame(inter_lt_f$pf,data=data_sim[[i]], na.action = na.fail)


    # Extract covariates matrix
    D[[i]] <- as.matrix(model.matrix(attr(mf,"terms"),data=data_sim[[i]]))
    n[[i]] <- nrow(D[[i]])
  }


  if(length(inter_f$re.spec) > 0) {
    stop("In the current impletementation of 'surf_sim' the addition of random effects
        with re() is not supported")
  }
  mf_pred <- model.frame(inter_lt_f$pf,data=pred_grid, na.action = na.fail)
  D_pred <- as.matrix(model.matrix(attr(mf_pred,"terms"),data=pred_grid))
  if(ncol(D_pred)!=ncol(D[[1]])) stop("the provided variables in 'grid_pred' do not match the number of explanatory variables used to fit the model.")

  # Number of covariates
  p <- ncol(D[[1]])

  beta <- par0$beta

  if(length(beta)!=p) stop("The values passed to 'beta' do not match the variables in 'grid_pred'")

  sigma2 <- par0$sigma2
  phi <- par0$phi

  if(is.null(par0$tau2)) {
    tau2 <- 0
  } else {
    tau2 <- par0$tau2
  }

  if(is.null(fix_var_me)) {
    sigma2_me <- 0
  } else {
    sigma2_me <- par0$sigma2_me
  }

  ID_coords <- sapply(1:n_sim, function(i) 1:n[[i]])

  grid_pred <- st_coordinates(pred_grid)
  if(scale_to_km) grid_pred <- grid_pred/1000
  n_pred <- nrow(grid_pred)

  # Simulate on the grid
  # Simulate S
  Sigma <- sigma2*matern_cor(dist(grid_pred), phi = phi, kappa = kappa,
                             return_sym_matrix = TRUE)
  if(nugget_over_grid) {
    diag(Sigma) <- diag(Sigma) + tau2
  }
  Sigma_sroot <- t(chol(Sigma))

  S_sim <- sapply(1:n_sim, function(i) Sigma_sroot%*%rnorm(n_pred))

  sim_columns <- paste0("sim_", 1:n_sim)
  S_sim_df <- as.data.frame(S_sim)
  colnames(S_sim_df) <- sim_columns

  # Combine simulations with grid_pred
  S_sim <- cbind(grid_pred, S_sim_df)
  S_sim <- st_sf(S_sim, geometry = st_geometry(pred_grid))

  coords_tot <- list()
  out <- list()

  for(i in 1:n_sim) {
    coords_tot[[i]] <- rbind(coords_sim[[i]], grid_pred)
    ind_data <- 1:n[[i]]
    ind_pred <- (n[[i]]+1):(n[[i]]+n_pred)

    D_tot <- rbind(D[[i]], D_pred)

    # Find the nearest features in grid_pred_with_sim for each location in data_sim[[i]]
    nearest_indices <- st_nearest_feature(data_sim[[i]],
                                          S_sim)

    # Extract the i-th simulation values
    sim_column <- paste0("sim_", i)  # Column name for the i-th simulation
    S_sim_data <- S_sim[nearest_indices, sim_column, drop = TRUE]

    if(tau2>0 & !nugget_over_grid) {
      S_sim_data <- S_sim_data + sqrt(tau2)*rnorm(n[[i]])
    }
    # Linear predictor
    if(include_covariates) {
      eta_sim_tot <- D_tot%*%beta +
        c(S_sim_data, S_sim[[sim_column]])
    } else {
      eta_sim_tot <- c(S_sim_data, S_sim[[sim_column]])
    }

    data_sim[[i]]$y <- NA

    if(family=="gaussian") {

      data_sim[[i]]$y <- eta_sim_tot[1:n[[i]]] + sqrt(sigma2_me)*rnorm(n)

    } else if(family=="binomial") {
      prob_i <- exp(eta_sim_tot[1:n[[i]]])/(1+exp(eta_sim_tot[1:n[[i]]]))
      data_sim[[i]]$y <- rbinom(n[[i]], size = data_sim[[i]]$units_m,
                                prob = prob_i)
    } else if(family=="poisson") {
      mean_i <- data_sim[[i]]$units_m*exp(eta_sim_tot[1:n[[i]]])
      data_sim[[i]]$y <- rpois(n[[i]], lambda = mean_i)
    }

    data_sim[[i]]$lp_data <- S_sim_data

    pred_grid[[paste0("lp_",sim_column)]] <- eta_sim_tot[-(1:n[[i]])]

  }

  out$data_sim <- data_sim
  out$lp_grid_sim <- pred_grid
  out$include_covariates <- include_covariates
  out$nugget_over_grid <- nugget_over_grid
  out$include_cov_offset <- include_cov_offset
  out$par0 <- par0
  out$family <- family
  class(out) <- "RiskMap.sim"
  return(out)
}

##' Plot simulated surface data for a given simulation
##'
##' This function plots the simulated surface data for a specific simulation from the result of `surf_sim`. It visualizes the linear predictor values on a raster grid along with the actual data points.
##'
##' @param surf_obj The output object from `surf_sim`, containing both simulated data (`data_sim`) and predicted grid simulations (`lp_grid_sim`).
##' @param sim The simulation index to plot.
##' @param ... Additional graphical parameters to be passed to the plotting function of the `terra` package.
##'
##' @return A plot of the simulation results.
##'
##' @importFrom stars st_rasterize
##'
##' @export
plot_sim_surf <-  function(surf_obj, sim, ...) {

  sf_object <- surf_obj$lp_grid_sim
  value_column <- paste0("lp_sim_",sim)
  r <- rast(st_rasterize(sf_object[,c("x",value_column)]))
  r[r == 0] <- NA

  plot(r, main = paste("Simulation no.", sim), ...)
  points(st_coordinates(surf_obj$data_sim[[sim]]), pch = 20)

}

##' @title Assess Simulations
##'
##' @description This function evaluates the performance of models based on simulation results from the `surf_sim` function.
##'
##' @param obj_sim An object of class `RiskMap.sim`, obtained as an output from the `surf_sim` function.
##' @param models A named list of models to be evaluated.
##' @param control_mcmc A control object for MCMC sampling, created with `set_control_sim()`. Default is `set_control_sim()`.
##' @param spatial_scale The scale at which predictions are assessed, either `"grid"` or `"area"`.
##' @param messages Logical, if `TRUE` messages will be displayed during processing. Default is `TRUE`.
##' @param f_grid_target A function for processing grid-level predictions.
##' @param f_area_target A function for processing area-level predictions.
##' @param shp A shapefile of class `sf` or `data.frame` for area-level analysis, required if `spatial_scale = "area"`.
##' @param col_names Column name in `shp` containing unique region names. If `NULL`, defaults to `"region"`.
##' @param pred_objective A character vector specifying objectives, either `"mse"`, `"classify"`, or both.
##' @param categories A numeric vector of thresholds defining categories for classification. Required if `pred_objective = "classify"`.
##'
##' @return A list of class `RiskMap.sim.res` containing model evaluation results.
##'
##' @export
assess_sim <- function(obj_sim,
                       models,
                       control_mcmc = set_control_sim(),
                       spatial_scale,
                       messages = TRUE,
                       f_grid_target = NULL,
                       f_area_target = NULL,
                       shp = NULL, col_names = NULL,
                       pred_objective = c("mse","classify"),
                       categories= NULL) {

  if (!inherits(obj_sim, "RiskMap.sim")) {
    stop("'obj_sim' must be an object of class 'RiskMap.sim' obtained as an output from the 'surf_sim' function")
  }
  if (length(setdiff(pred_objective, c("mse","classify")))>0) {
    stop(paste("Invalid value for pred_objective. Allowed values are:", paste(c("mse","classify"), collapse = ", ")))
  }
  if(spatial_scale != "grid" & spatial_scale != "area") {
    stop("'spatial_scale' must be set to 'grid' or 'area'")
  }
  if(spatial_scale=="area" & is.null(shp)) {
    stop("if spatial_scale='area' then a shape file of the area(s) must be passed to
         'shp'")
  }
  units_m <- NULL
  if(any(pred_objective=="classify")) {
    if(is.null(categories)) stop("if pred_objective='class', a value for 'categories' must be specified")
    if (length(categories) < 3) {
      stop("Categories vector must contain at least three unique, strictly increasing values.")
    }
  }
  n_sim <- length(obj_sim$data_sim)
  n_models <- length(models)

  if(spatial_scale == "area" & is.null(f_area_target)) {
    stop("If 'spatial_scale' is set to 'area', then 'f_area_target' must be provided")
  }
  model_names <- names(models)

  include_covariates <- obj_sim$include_covariates
  include_cov_offset <- obj_sim$include_cov_offset
  include_nugget <- obj_sim$nugget_over_grid

  fits <- list()
  preds <- list()
  if(spatial_scale=="grid") {
    type <- "marginal"
  } else if(spatial_scale=="area") {
    type <- "joint"
    n_reg <- nrow(shp)
    if(is.null(shp)) stop("If spatial_scale='area', then 'shp' must be specified")
    if(!inherits(shp,
                 what = c("sf","data.frame"), which = FALSE)) {
      stop("The object passed to 'shp' must be an object of class 'sf'")
    }

    if(is.null(col_names)) {
      shp$region <- paste("reg",1:n_reg, sep="")
      col_names <- "region"
      names_reg <- shp$region
    } else {
      names_reg <- shp[[col_names]]
      if(n_reg != length(names_reg)) {
        stop("The names in the column identified by 'col_names' do not
         provide a unique set of names, but there are duplicates")
      }
    }
    shp <- st_transform(shp, st_crs(obj_sim$lp_grid_sim))
    inter <- st_intersects(shp, obj_sim$lp_grid_sim)
  }

  for(i in 1:n_models) {
    if(messages) message("Model: ", paste(model_names[i]),"\n")

    if_i <- interpret.formula(models[[i]])
    rhs_terms <- attr(terms(if_i$pf), "term.labels")
    # Check if there are any covariates
    if (length(rhs_terms) == 0) {
      predictors_i <- NULL
    } else {
      predictors_i <- obj_sim$lp_grid_sim
    }
    for(j in 1:n_sim) {
      if(messages) message("Processing simulation no.", j)
      f_i <- update(models[[i]], y ~ .)
      if(messages) message("Estimation")
      fits[[paste(model_names[i])]][[j]] <- glgpm(formula = f_i,
                                                  den = units_m,
                                                  family = obj_sim$family,
                                                  data = obj_sim$data_sim[[j]],
                                                  control_mcmc = control_mcmc,
                                                  messages = FALSE)

      if(messages) message("Prediction over the grid")
      preds[[paste(model_names[i])]][[j]] <-
        pred_over_grid(fits[[paste(model_names[i])]][[j]],
                       grid_pred = st_as_sfc(obj_sim$lp_grid_sim),
                       predictors = predictors_i,
                       type = type, messages = FALSE)
    }
  }

  n_samples <- (control_mcmc$n_sim-control_mcmc$burnin)/control_mcmc$thin
  n_pred <- nrow(obj_sim$lp_grid_sim)


  out <- list(pred_objective = list())

  if(any(pred_objective=="mse")) {
    out$pred_objective$mse <- array(NA, c(n_models, n_sim))
    rownames(out$pred_objective$mse) <- model_names
    colnames(out$pred_objective$mse) <- paste0("sim_",1:n_sim)
  }

  if(any(pred_objective=="classify")) {

    # Ensure categories are unique and strictly increasing
    categories <- unique(sort(categories))


    # Assign classification to the output object
    out$pred_objective$classify <- setNames(vector("list", length(model_names)), model_names)

    # Correctly generate breaks and labels
    breaks <- categories  # Use categories directly as breaks
    categories_class <- factor(paste0("(", head(categories, -1), ",",
                                      categories[-1], "]"))  # Labels to match intervals



    for(i in 1:n_models) {
      out$pred_objective$classify[[model_names[i]]] <- list(by_cat = list(),
                                                            across_cat = list())
      out$pred_objective$classify[[model_names[i]]]$by_cat <- vector("list", n_sim)
      for(j in 1:n_sim) {
        out$pred_objective$classify[[paste(model_names[i])]]$by_cat[[j]] <-
          data.frame(
            Class = categories_class,
            Sensitivity = NA,
            Specificity = NA,
            PPV = NA,
            NPV = NA,
            CC = NA
          )
      }
      out$pred_objective$classify[[model_names[i]]]$CC <- rep(NA,n_sim)
    }
  }
  lp_true_sim <- st_drop_geometry(obj_sim$lp_grid_sim[, grepl("lp_sim",
                                                              names(obj_sim$lp_grid_sim))])

  if(spatial_scale == "grid") {
    true_target_sim <- f_grid_target(lp_true_sim)
  } else if(spatial_scale == "area") {
    true_target_sim <- matrix(NA, nrow = n_reg, ncol = n_sim)
    true_target_grid_sim <- f_grid_target(lp_true_sim)
    for(i in 1:n_reg) {
      for(j in 1:n_sim) {
        if(length(inter[[i]])==0) {
          warning(paste("No points on the grid fall within", shp[[col_names]][h],
                        "and no predictions are carried out for this area"))
          no_comp <- c(no_comp, h)
        } else {
          true_target_sim[i,j] <- f_area_target(true_target_grid_sim[inter[[i]],j])
        }
      }
    }
  }



  for(i in 1:n_models) {
    for(j in 1:n_sim) {
      obj_pred_ij <- preds[[paste(model_names[i])]][[j]]
      if(length(obj_pred_ij$mu_pred)==1 && obj_pred_ij$mu_pred==0 &&
         include_covariates) {
        stop("Covariates have not been provided; re-run pred_over_grid
         and provide the covariates through the argument 'predictors'")
      }


      if(!include_covariates) {
        mu_target <- 0
      } else {

        if(is.null(obj_pred_ij$mu_pred)) stop("the output obtained from 'pred_S' does not
                                     contain any covariates; if including covariates
                                     in the predictive target these shuold be included
                                     when running 'pred_S'")
        mu_target <- obj_pred_ij$mu_pred
      }

      if(!include_cov_offset) {
        cov_offset <- 0
      } else {
        if(length(obj_pred_ij$cov_offset)==1) {
          stop("No covariate offset was included in the model;
           set include_cov_offset = FALSE, or refit the model and include
           the covariate offset")
        }
        cov_offset <- obj_pred_ij$cov_offset
      }

      if(include_nugget) {
        if(is.null(obj_pred_ij$par_hat$tau2)) stop("'include_nugget' cannot be
                                                   set to TRUE if this has not been included
                                                   in the fit of the model")
        Z_sim <- matrix(rnorm(n_samples*n_pred,
                              sd = sqrt(obj_pred_ij$par_hat$tau2)),
                        ncol = n_samples)
        obj_pred_ij$S_samples <- obj_pred_ij$S_samples+Z_sim
      }

      if(is.matrix(mu_target)) {
        lp_samples_ij <- sapply(1:n_samples,
                                function(h)
                                  mu_target[,h] + cov_offset +
                                  obj_pred_ij$S_samples[,h])
      } else {
        lp_samples_ij <- sapply(1:n_samples,
                                function(h)
                                  mu_target + cov_offset +
                                  obj_pred_ij$S_samples[,h])
      }

      target_samples_ij <- f_grid_target(lp_samples_ij)

      if(spatial_scale == "grid") {
        mean_target_ij <- apply(target_samples_ij, 1, mean)
      } else if(spatial_scale == "area") {
        target_area_samples_ij <- matrix(NA, nrow = n_reg, ncol = n_samples)
        mean_target_ij <- rep(NA,n_reg)
        for(h in 1:n_reg) {
          if(length(inter[[h]])==0) {
            warning(paste("No points on the grid fall within", shp[[col_names]][h],
                          "and no predictions are carried out for this area"))
            no_comp <- c(no_comp, h)
          } else {
            ind_grid_h <- inter[[h]]
            target_area_samples_ij[h,] <-  apply(target_samples_ij[ind_grid_h,], 2,
                                                 f_area_target)
            mean_target_ij[h] <- mean(target_area_samples_ij[h,])
          }
        }
      }

      if(any(pred_objective=="mse")) {
        out$pred_objective$mse[i,j] <- mean((mean_target_ij-true_target_sim[,j])^2)
      }

      if(any(pred_objective=="classify")) {
        true_class_ij <- cut(true_target_sim[,j], breaks = categories)
        n_categories <- length(categories)-1
        if(spatial_scale == "grid") {
          prob_cat_ij <- matrix(0, nrow=n_pred, ncol = n_categories)
        } else if(spatial_scale == "area") {
          prob_cat_ij <- matrix(0, nrow=n_reg, ncol = n_categories)
        }
        for(h in 1:(n_categories)) {
          if(spatial_scale == "grid") {
            prob_cat_ij[,h] <- apply(categories[h] < target_samples_ij &
                                       categories[h+1] > target_samples_ij, 1, mean)
          } else if(spatial_scale == "area") {
            prob_cat_ij[,h] <- apply(categories[h] < target_area_samples_ij &
                                       categories[h+1] > target_area_samples_ij, 1, mean)
          }
        }

        pred_class_ij <- apply(prob_cat_ij, 1, function(x) categories_class[which.max(x)])

        # Define the confusion matrix
        conf_matrix <- table(true_class_ij, pred_class_ij)

        # Calculate metrics for each class
        for (h in 1:nrow(conf_matrix)) {
          TP <- conf_matrix[h, h]  # True Positives: diagonal entry
          FP <- sum(conf_matrix[, h]) - conf_matrix[h, h]  # False Positives: column sum minus diagonal
          FN <- sum(conf_matrix[h, ]) - conf_matrix[h, h]  # False Negatives: row sum minus diagonal
          TN <- sum(conf_matrix) - sum(conf_matrix[h, ]) - sum(conf_matrix[, h]) + conf_matrix[h, h]


          # Handle cases where division by zero could occur
          out$pred_objective$classify[[paste(model_names[[i]])]]$by_cat[[j]][h,]$Sensitivity <-
            ifelse((TP + FN) == 0, NA, TP / (TP + FN))
          out$pred_objective$classify[[paste(model_names[[i]])]]$by_cat[[j]][h,]$Specificity <-
            ifelse((TN + FP) == 0, NA, TN / (TN + FP))
          out$pred_objective$classify[[paste(model_names[[i]])]]$by_cat[[j]][h,]$PPV <-
            ifelse((TP + FP) == 0, NA, TP / (TP + FP))
          out$pred_objective$classify[[paste(model_names[[i]])]]$by_cat[[j]][h,]$NPV <-
            ifelse((TN + FN) == 0, NA, TN / (TN + FN))
          out$pred_objective$classify[[paste(model_names[[i]])]]$by_cat[[j]][h,]$CC <-
            ifelse((TP + FN) == 0, NA, (TP ) / (TP + FN))
        }
        out$pred_objective$classify[[paste(model_names[[i]])]]$CC[j] <-
          mean(true_class_ij==pred_class_ij)
      }
    }
  }
  if(any(pred_objective=="classify")) out$pred_objective$classify$Class <- categories_class
  class(out) <- "RiskMap.sim.res"
  return(out)
}

##' @title Summarize Simulation Results
##'
##' @description Summarizes the results of model evaluations from a `RiskMap.sim.res` object. Provides average metrics for classification by category and overall correct classification (CC) summary.
##'
##' @param object An object of class `RiskMap.sim.res`, as returned by `assess_sim`.
##' @param ... Additional arguments (not used).
##'
##' @return A list containing summary data for each model:
##' - `by_cat_summary`: A data frame with average sensitivity, specificity, PPV, NPV, and CC by category.
##' - `CC_summary`: A numeric vector with mean, 2.5th percentile, and 97.5th percentile for CC across simulations.
##'
##' @method summary RiskMap.sim.res
##' @export
summary.RiskMap.sim.res <- function(object, ...) {
  stopifnot(inherits(object, "RiskMap.sim.res"))

  # Initialize results
  results <- list()

  # Check for "mse" in pred_objective
  if ("mse" %in% names(object$pred_objective)) {
    mse_data <- object$pred_objective$mse

    # Check if mse_data is a matrix
    if (is.matrix(mse_data)) {
      # Compute mean and SD for each model (row)
      mse_summary <- data.frame(
        Model = rownames(mse_data),
        MSE_mean = rowMeans(mse_data, na.rm = TRUE),
        MSE_sd = apply(mse_data, 1, sd, na.rm = TRUE)
      )

      results$mse <- mse_summary
    } else {
      stop("mse_data must be a matrix.")
    }
  }

  # Check for "classify" in pred_objective
  if ("classify" %in% names(object$pred_objective)) {
    classify_data <- object$pred_objective$classify

    # Loop over each model (e.g., M1, M2)
    n_models <- length(classify_data)-1
    name_models <- names(classify_data)[1:n_models]
    results$classify <- list()
    for(i in 1:n_models) {

      model_data <- classify_data[[i]]

      n_sim <- length(model_data$by_cat)
      res_class <- model_data$by_cat[[1]][,-1]
      den <- 0
      for(j in 2:n_sim) {
        if(!any(is.na(model_data$by_cat[[j]][,-1]))) {
          den <- den + 1
          res_class <- res_class+model_data$by_cat[[j]][,-1]
        }
      }
      res_class <- data.frame(res_class/den)
      res_class$Class <- model_data$by_cat[[1]][,1]

      cc_summary <- list(mean = mean(model_data$CC, na.rm = TRUE),
                         lower = quantile(model_data$CC, 0.025, na.rm = TRUE),
                         upper = quantile(model_data$CC, 0.975, na.rm = TRUE))

      results$classify[[paste(name_models[i])]] <- list(classify_res = res_class,
                                            cc_summary = list(mean = mean(model_data$CC, na.rm = TRUE),
                                                                     lower = quantile(model_data$CC, 0.025, na.rm = TRUE),
                                                                     upper = quantile(model_data$CC, 0.975, na.rm = TRUE)))
    }
  }

  # Assign class for S3 print method
  class(results) <- "summary.RiskMap.sim.res"
  return(results)
}



##' @title Print Simulation Results
##'
##' @description Prints a concise summary of simulation results from a `RiskMap.sim.res` object, including average metrics by category and a summary of overall correct classification (CC).
##'
##' @param x An object of class `summary.RiskMap.sim.res`, as returned by `summary.RiskMap.sim.res`.
##' @param ... Additional arguments (not used).
##'
##' @return Invisibly returns `x`.
##'
##'
##' Print Simulation Results
##'
##' Prints a concise summary of simulation results from a `summary.RiskMap.sim.res` object,
##' including average metrics by category and a summary of overall correct classification (CC).
##'
##' @param x An object of class `summary.RiskMap.sim.res`, as returned by `summary.RiskMap.sim.res`.
##' @param ... Additional arguments (not used).
##'
##' @return Invisibly returns `x`.
##'
##' @method print summary.RiskMap.sim.res
##' @export
print.summary.RiskMap.sim.res <- function(x, ...) {
  cat("Summary of Simulation Results\n\n")

  if (!is.null(x$mse)) {
    cat("Mean Squared Error (MSE):\n")
    print(x$mse)
    cat("\n")
  }

  if (!is.null(x$classify)) {
    cat("Classification Results:\n")

    # Iterate over each model in classify results
    for (model_name in names(x$classify)) {
      model_data <- x$classify[[model_name]]

      cat(sprintf("\nModel: %s\n", model_name))

      cat("\nAverages across simulations by Category:\n")
      print(model_data$classify_res)

      cat("\nProportion of Correct Classification (CC) across categories:\n")
      cc_summary <- model_data$cc_summary
      cat(sprintf("Mean: %.3f, 95%% CI: [%.3f, %.3f]\n",
                  cc_summary$mean, cc_summary$lower, cc_summary$upper))
    }
    cat("\n")
  }

  invisible(x)
}

