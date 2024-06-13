##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @importFrom Matrix Matrix forceSymmetric
##' @export
pred_S <- function(object, predictors = NULL,
                   grid_pred,
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

  inter_f <- interpret.formula(object$formula)
  inter_lt_f <- inter_f
  inter_lt_f$pf <- update(inter_lt_f$pf, NULL ~.)

  if(p==1) {
    intercept_only <- prod(object$D[,1]==1)
  } else {
    intercept_only <- FALSE
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

  out <- list()

  out$mu_pred <- mu_pred
  out$grid_pred <- grid_pred
  out$par_hat <- object$par_hat

  if(object$scale_to_km) grp <- grp/1000
  U_pred <- t(sapply(1:n_pred,
                     function(i) sqrt((object$coords[,1]-grp[i,1])^2+
                                        (object$coords[,2]-grp[i,2])^2)))
  U <- dist(object$coords)
  C <- par_hat$sigma2*matern_cor(U_pred, phi = par_hat$phi,
                                 kappa = object$kappa)
  Sigma <- par_hat$sigma2*matern_cor(U, phi = par_hat$phi,
                                     kappa = object$kappa,
                                     return_sym_matrix = TRUE)
  Sigma_inv <- solve(Sigma)
  A <- C%*%Sigma_inv
  mu <- as.numeric(object$D%*%par_hat$beta)

  if(control_sim$linear_model) {
    n_samples <- control_sim$n_sim
  } else {
    n_samples <- (control_sim$n_sim-control_sim$burnin)/control_sim$thin
  }

  if(object$family!="gaussian") {
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

      Sigma_cond <- Sigma_pred - A%*%t(C)
      Sigma_cond_sroot <- t(chol(Sigma_cond))
      out$S_samples <- sapply(1:n_samples,
                              function(i)
                                mu_cond_S[,i]+
                                Sigma_cond_sroot%*%rnorm(n_pred))
    }
  } else {
    mu_cond_S <- as.numeric(A%*%(object$y-mu))
    if(type=="marginal") {
      sd_cond_S <- sqrt(par_hat$sigma2-diag(A%*%t(C)))
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

  out$inter_f <- inter_f
  out$family <- object$family
  out$par_hat <- par_hat
  out$cov_offset <- pred_cov_offset
  class(out) <- "RiskMap.pred.re"
  return(out)
}

##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @importFrom Matrix Matrix forceSymmetric
##' @export
pred_target <- function(object,
                        include_covariates = TRUE,
                        include_nugget = FALSE,
                        include_cov_offset = FALSE,
                        f_target = NULL,
                        pd_summary = NULL,
                        shp = NULL) {
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

  out <- list()
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
  class(out) <- "RiskMap.pred.target"
  return(out)
}

##' @importFrom terra rast plot as.data.frame
##' @method coef RiskMap.pred.target
##' @export
##'
plot.RiskMap.pred.target <- function(object,
                                     which_target,
                                     which_summary, ...) {
  t_data.frame <-
    terra::as.data.frame(cbind(st_coordinates(object$grid_pred),
                               object$target[[which_target]][[which_summary]]),
                         xy= TRUE)
  raster_out <- terra::rast(t_data.frame, crs = st_crs(object$grid_pred)$input)

  terra::plot(raster_out)
}
