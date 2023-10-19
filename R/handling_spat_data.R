##' @title Create grid of points within shape file
##' @export
create_grid <- function(shp, spat_res,
                        grid_crs = NULL) {

  if(class(shp)[1]!="sf") stop("'shp' must be an object of class 'sf'")

  if(is.na(st_crs(shp))) stop("The CRS for 'shp' is missing")

  if(is.null(grid_crs)) {
    grid_crs <- st_crs(shp)
  } else {
    shp <- st_transform(shp, crs = grid_crs)
  }
  grid_box <- st_make_grid(shp,
                           cellsize = spat_res*1000,
                           what="centers")

  liberia.inout <- st_intersects(grid_box,
                                 shp,
                                 sparse = FALSE)
  grid_out <- grid_box[liberia.inout]
  return(grid_out)
}
