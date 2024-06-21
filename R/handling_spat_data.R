##' @title Create Grid of Points Within Shapefile
##'
##' @description
##' Generates a grid of points within a given shapefile. The grid points are created based on a specified spatial resolution.
##'
##' @param shp An object of class 'sf' representing the shapefile within which the grid of points will be created.
##' @param spat_res Numeric value specifying the spatial resolution in kilometers for the grid.
##' @param grid_crs Coordinate reference system for the grid. If NULL, the CRS of 'shp' is used. The shapefile 'shp' will be transformed to this CRS if specified.
##'
##' @details
##' This function creates a grid of points within the boundaries of the provided shapefile ('shp'). The grid points are generated using the specified spatial resolution ('spat_res'). If a coordinate reference system ('grid_crs') is provided, the shapefile is transformed to this CRS before creating the grid.
##'
##' @return
##' An 'sf' object containing the generated grid points within the shapefile.
##'
##' @importFrom sf st_make_grid st_intersects st_transform st_crs
##' @export
##'
##' @examples
##' library(sf)
##'
##' # Example shapefile data
##' nc <- st_read(system.file("shape/nc.shp", package="sf"))
##'
##' # Create grid with 10 km spatial resolution
##' grid <- create_grid(nc, spat_res = 10)
##'
##' # Plot the grid
##' plot(st_geometry(nc))
##' plot(grid, add = TRUE, col = 'red')
##'
##' @seealso
##' \code{\link[sf]{st_make_grid}}, \code{\link[sf]{st_intersects}}, \code{\link[sf]{st_transform}}, \code{\link[sf]{st_crs}}
##'
##' @author Emanuele Giorgi \email{e.giorgi@@lancaster.ac.uk}
##' @author Claudio Fronterre \email{c.fronterr@@lancaster.ac.uk}
##'
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
