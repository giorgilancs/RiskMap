
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
