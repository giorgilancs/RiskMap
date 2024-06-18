##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @importFrom Matrix Matrix forceSymmetric
##' @export
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
  grid_pred <- st_transform(grid_pred,crs = object$crs)

  grp <- st_coordinates(grid_pred)
  n_pred <- nrow(grp)
  par_hat <- coef(object)
  p <- ncol(object$D)

  if(length(object$cov_offset) > 1) {
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
    intercept_only <- prod(object$D[,1]==1)
  } else {
    intercept_only <- FALSE
  }

  n_re <- length(object$re)
  if(n_re > 0 & type=="marginal") {
    warning("Predictions for the unstructured random effects will not be perfomed if type = 'marginal';
            if you wish to also include the predictions for the unstructured random
            effects then set type = 'joint'")
    n_re <- 0
  }

  if(!is.null(predictors)) {

    if(!is.data.frame(predictors)) stop("'predictors' must be an object of class 'data.frame'")
    if(nrow(predictors)!=n_pred) stop("the values provided for 'predictors' do not match the prediction grid passed to 'grid_pred'")

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
      if(!is.data.frame(re_predictors)) stop("'re_predictors' must be an object of class 'data.frame'")
      if(nrow(re_predictors)!=n_pred) stop("the values provided for 'predictors' do not match the prediction grid passed to 'grid_pred'")
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
  diag(R) <- diag(R)+nu2

  if(object$family!="gaussian") {
    Sigma <- par_hat$sigma2*R
    Sigma_inv <- solve(Sigma)
    A <- C%*%Sigma_inv

    simulation <-
      Laplace_sampling_MCMC(y = object$y, units_m = object$units_m, mu = mu, Sigma = Sigma,
                            sigma2_re = par_hat$sigma2_re,
                            ID_coords = object$ID_coords, ID_re = object$ID_re,
                            family = object$family, control_mcmc = control_sim,
                            messages = messages)
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
      diag(Sigma_pred) <- diag(Sigma_pred)+10e-10
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

  if(n_re>0) {
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
  }

  out$inter_f <- inter_f
  out$family <- object$family
  out$par_hat <- par_hat
  out$cov_offset <- pred_cov_offset
  out$type <- type
  class(out) <- "RiskMap.pred.re"
  return(out)
}

##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @importFrom Matrix Matrix forceSymmetric
##' @export
pred_target_grid <- function(object,
                        include_covariates = TRUE,
                        include_nugget = FALSE,
                        include_cov_offset = FALSE,
                        include_re = FALSE,
                        f_target = NULL,
                        pd_summary = NULL) {
  if(!inherits(object,
               what = "RiskMap.pred.re", which = FALSE)) {
    stop("The object passed to 'object' must be an output of
         the function 'pred_S'")
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


  out$lp_samples <- sapply(1:n_samples,
                           function(i)
                             mu_target + cov_offset +
                             object$S_samples[,i])

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

  for(i in 1:n_f) {
    target_samples_i <-
      f_target[[i]](out$lp_samples)
    out$target[[paste(names_f[i])]] <- list()
    for(j in 1:n_summaries) {
      out$target[[paste(names_f[i])]][[paste(names_s[j])]] <-
        apply(target_samples_i, 1, pd_summary[[j]])
    }
  }
  out$grid_pred <- object$grid_pred
  out$f_target <- names(f_target)
  out$pd_summary <- names(pd_summary)
  class(out) <- "RiskMap.pred.target.grid"
  return(out)
}

##' @importFrom terra rast plot as.data.frame
##' @method plot RiskMap.pred.target.grid
##' @export
##'
plot.RiskMap.pred.target.grid <- function(object,
                                     which_target,
                                     which_summary, ...) {
  t_data.frame <-
    terra::as.data.frame(cbind(st_coordinates(object$grid_pred),
                               object$target[[which_target]][[which_summary]]),
                         xy= TRUE)
  raster_out <- terra::rast(t_data.frame, crs = st_crs(object$grid_pred)$input)

  terra::plot(raster_out)
}

##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @export
pred_target_shp <- function(object, shp, shp_target,
                            weights = NULL,standardize_weights = FALSE,
                            col_names=NULL,
                            include_covariates = TRUE,
                            include_nugget = FALSE,
                            include_cov_offset = FALSE,
                            include_re = FALSE,
                            f_target = NULL,
                            pd_summary = NULL) {
  if(!inherits(object,
               what = "RiskMap.pred.re", which = FALSE)) {
    stop("The object passed to 'object' must be an output of
         the function 'pred_S'")
  }
  if(object$type!="joint") {
    stop("To run predictions with a shape file, joint predictions must be used;
         rerun 'pred_over_grid' and set type='joint'")
  }

  if(!inherits(shp,
               what = c("sf","data.frame"), which = FALSE)) {
    stop("The object passed to 'shp' must be an output of
         the function 'pred_S'")
  }

  if(!inherits(shp_target,
               what = c("function"), which = FALSE)) {
    stop("The object passed to 'shp_function' must be a function
         that takes as an input a vector and returns a single numeric value")
  }
  n_pred <- nrow(object$S_samples)

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
  out$lp_samples <- sapply(1:n_samples,
                           function(i)
                             mu_target + cov_offset +
                             object$S_samples[,i])
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

  for(h in 1:n_reg) {
    cat("Computing predictive target for:",shp[[col_names]][h])
    if(length(inter[[h]])==0) {
      stop(paste("No points on the grid fall within", shp[[col_names]][h]))
    }
    ind_grid_h <- inter[[h]]
    if(standardize_weights & !no_weights) {
      weights_h <- weights[ind_grid_h]/sum(weights[ind_grid_h])
    } else {
      weights_h <- weights[ind_grid_h]
    }
    for(i in 1:n_f) {

      target_samples_i <-
        apply(f_target[[i]](out$lp_samples[ind_grid_h,]),
              2, function(x) shp_target(weights_h*x))
      out$target[[paste(names_reg[h])]][[paste(names_f[i])]] <- list()
      for(j in 1:n_summaries) {
        out$target[[paste(names_reg[h])]][[paste(names_f[i])]][[paste(names_s[j])]] <-
          pd_summary[[j]](target_samples_i)
      }
    }
    cat(" \n")
  }

  for(i in 1:n_f) {
    for(j in 1:n_summaries) {
      name_ij <- paste(names_f[i],"_",paste(names_s[j]),sep="")
      shp[[name_ij]] <- rep(NA, n_reg)
      for(h in 1:n_reg) {
        which_reg <- which(shp[[col_names]]==names_reg[h])
        shp[which_reg,][[name_ij]] <-
          out$target[[paste(names_reg[h])]][[paste(names_f[i])]][[paste(names_s[j])]]
      }
    }
  }
  out$shp <- shp
  out$f_target <- names(f_target)
  out$pd_summary <- names(pd_summary)
  out$grid_pred <- object$grid_pred
  class(out) <- "RiskMap.pred.target.shp"
  return(out)
}


##' @importFrom ggplot2 ggplot geom_sf aes scale_fill_distiller
##' @method plot RiskMap.pred.target.shp
##' @export
##'
plot.RiskMap.pred.target.shp <- function(object,
                                          which_target,
                                          which_summary, ...) {
  col_shp_name <- paste(which_target,"_",which_summary,sep="")

  out <- ggplot(object$shp) +
    geom_sf(aes(fill = .data[[col_shp_name]])) +
    scale_fill_distiller(...)
  return(out)
}
