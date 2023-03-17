
##' @title Empirical variogram
##' @description Computes the empirical variogram using ``bins'' of distance provided by the user.
##' @param data an object of class \code{sf} containing the variable for which the variogram
##' is to be computed and the coordinates
##' @param variable a character indicating the name of variable for which the variogram is to be computed.
##' @param bins a vector indicating the `bins` to be used to define the classes of distance used in the computation of the variogram.
##' By default \code{bins=NULL} and bins are then computed as \code{seq(0, d_max/2, length=15)} where \code{d_max} is the maximum distance observed in the data.
##' @param n_permutation a non-negative integer indicating the number of permutation used to compute the 95% confidence
##' level envelope under the assumption of spatial independence. By default \code{n_permutation=0}, and no envelope is generated.
##' @param convert_to_utm a logical value, indicating if the conversion to UTM shuold be performed (\code{convert_to_utm = TRUE}) or
##' the coordinate reference system of the data must be used without any conversion (\code{convert_to_utm = FALSE}).
##' By default \code{convert_to_utm = TRUE}. Note: if \code{convert_to_utm = TRUE} the conversion to UTM is performed using
##' the epsg provided by \code{\link{propose_utm}}.
##' @param scale_to_km a logical value, indicating if the distances used in the variogram must be scaled
##' to kilometers (\code{scale_to_km = TRUE}) or left in meters (\code{scale_to_km = FALSE}).
##' By default \code{scale_to_km = FALSE}
##'
##'
##' @return an object of class 'variogram' which is a list containing the following components
##' @return \code{variogram} a data-frame containing the following columns: \code{mid_points},
##' the middle points of the classes of distance provided by \code{bins};
##' \code{obs_vari} the values of the observed variogram; \code{obs_vari} the number of pairs.
##' If \code{n_permutation > 0}, the data-frame also contains \code{lower_bound} and \code{upper_bound}
##' corresponding to the lower and upper bounds of the 95% confidence intervals
##' used to assess the departure of the observed variogram from the assumption of spatial independence.
##'
##' @return \code{scale_to_km} the value passed to \code{scale_to_km}
##' @return \code{n_permutation} the number of permutations
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##'
##' @importFrom sf st_transform st_coordinates
##' @importFrom stats dist quantile
##' @export
##'

s_variogram <- function(data, variable, bins = NULL,
                      n_permutation = 0,
                      convert_to_utm = TRUE,
                      scale_to_km = FALSE) {

  if(class(data)[1]!="sf") stop("'data' must be an object of class 'sf'")
  if(n_permutation < 0 |
     n_permutation != round(n_permutation)) stop("n_permutation must be a positive integer number")

  if(!convert_to_utm) cat("The distances of the variogram are computed assuming
                          that the CRS of the data gives distances in meters or kilometers")
  data <- st_transform(data, crs = 4326)
  data <- st_transform(data, crs = propose_utm(data))
  coords <- st_coordinates(data)
  d <- as.numeric(dist(coords))
  if(scale_to_km) d <- d/1000
  v <- (as.numeric(dist(data[[variable]])) ^ 2) / 2
  vario_df <- data.frame(d=d, v=v)
  if(is.null(bins)) {
    max.v <- max(d)/2
    bins <- seq(0, max.v, length=15)
  } else {
    max.v <- max(bins)
  }
  if(max.v > max(d)) stop("the provided distances in 'bins' go beyond
                          the maximum observed distance")
  mid_points <- bins[-length(bins)] + (bins[2] - bins[1]) / 2

  vario_df <- vario_df[vario_df$d < max.v,]
  if(nrow(vario_df)==0) stop("the values provided in 'bins' do not match the
  scale of the obsererved distances; consider setting scale_to_km = TRUE")
  vario_df$dist_class <- cut(vario_df$d, breaks = bins,
                             include.lowest = TRUE, right = TRUE)
  variogram <- data.frame(mid_points = mid_points)
  variogram$obs_vari <- tapply(vario_df$v, vario_df$dist_class, mean)
  variogram$n_obs <- tapply(vario_df$v, vario_df$dist_class, length)
  variogram$n_obs[is.na(variogram$n_obs)] <- 0
  if(n_permutation > 0) {
    v_perm <- matrix(NA, nrow=length(variogram$obs_vari),
                     ncol = n_permutation)
    n <- nrow(data)
    ind_perm <- sample(1:n)
    obs_are <- variogram$n_obs>0
    for(i in 1:n_permutation) {
      ind_perm <- sample(1:n)
      v_i <- (as.numeric(dist(data[[variable]][ind_perm])) ^ 2) / 2

      vario_df_i <- data.frame(d=d, v=v_i)
      vario_df_i <- vario_df_i[vario_df_i$d < max.v,]
      vario_df_i$dist_class <- vario_df$dist_class

      v_perm[,i] <- tapply(vario_df_i$v, vario_df_i$dist_class, mean)
    }
    variogram$lower_bound <- NA
    variogram$lower_bound[obs_are] <- apply(v_perm[obs_are,], 1,
                                            function(x) quantile(x, 0.025))
    variogram$upper_bound <- NA
    variogram$upper_bound[obs_are] <- apply(v_perm[obs_are,], 1,
                                            function(x) quantile(x, 0.025))
  }
  result <- list(variogram = variogram)
  result$scale_to_km <- scale_to_km
  result$n_permutation <- n_permutation

  class(result) <- "RiskMap_variogram"
  return(result)
}


##' @title Plotting the empirical variogram
##' @description ...
##' @param variogram ...
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##' @method plot variogram.RiskMap
##' @importFrom methods is
##' @export
##'

plot.variogram.RiskMap <- function(variogram) {
  if(!is(variogram, "RiskMap_variogram")) stop("variogram must be an object of class 'RiskMap_variogram'")

}
