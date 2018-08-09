
#'@param data a data object of class 'ppj_data'
#'@param mesh.args. list of extra mesh arguments. guess this is a silly implementation actually. 

build_mesh <- function(data, mesh.args = NULL){

  stopifnot(inherits(data, 'ppj_data'))


  pars <- list(convex = -0.01,
               concave = -0.5,
               resolution = 300,
               max.edge = c(0.4, 8), 
               cut = 0.4, 
               offset = c(1, 15))

  pars[names(mesh.args)] <- mesh.args

  outline <- unionSpatialPolygons(data$shapefiles, IDs = rep(1, length(data$shapefiles)))
  plot(outline)

  # Join coords from all polygons
  coords <- list()
  for(i in seq_len(length(outline@polygons[[1]]@Polygons))){
    coords[[i]] <- outline@polygons[[1]]@Polygons[[i]]@coords
  }
  coords <- do.call(rbind, coords)
  coords <- rbind(coords, as.matrix(data$pr[, c('longitude', 'latitude')]))
    
  # Find a nonconvex hull around points.
  outline.hull <- inla.nonconvex.hull(coords, 
                                      convex = pars$convex, 
                                      concave = pars$concave,
                                      resolution = pars$resolution)



  mesh <- inla.mesh.2d( 
    boundary = outline.hull,
    max.edge = pars$max.edge, 
    cut = pars$cut, 
    offset = pars$offset)


    message('Number of nodes: ', mesh$n)


  p <- 
    ggplot_projection_shapefile(raster = NULL, 
                                spatialpolygons = outline, 
                                mesh = mesh, 
                                shapecol = 'red') +
      scale_colour_manual(values = c('#00000000', 'red', 'white', 'white', 'white'))

  print(p)


  return(mesh)


}
