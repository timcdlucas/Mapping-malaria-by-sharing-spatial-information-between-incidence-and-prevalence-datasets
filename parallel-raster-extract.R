
#' Parallel extraction of means of rasters by shape file.
#' 
#' Parallelisation is performed across rasters, not shapes. 
#' So this function is only useful if you are extracting 
#' data from many raster layers.
#' As the overhead for parallel computation in windows is high
#' it only makes sense to parallelise in this way.
#'
#' 
#' @param raster A raster brick or stack
#' @param shape A shape object 
#' @param fun The function used to aggregate the pixel data. If NULL, raw pixel data is returned.
#' @param id Name of column in shape object to be used to bind an ID column to output.
#' @param ... Other arguments to raster::extract
#' @importFrom foreach %dopar%
#' @importFrom snow stopCluster
#' @importFrom snow makeCluster
#' @importFrom doSNOW registerDoSNOW
#'
#' @export
#' @examples 
#' 
#' # Create raster stack
#' r <- raster::raster(ncol=36, nrow=18)
#' r[] <- 1:raster::ncell(r)
#' r <- stack(r, r)
#' 
#' # Create polygon
#' cds1 <- rbind(c(-180,-20), c(-160,5), c(-60, 0), c(-160,-60), c(-180,-20))
#' cds2 <- rbind(c(80,0), c(100,60), c(120,0), c(120,-55), c(80,0))
#' polys <- raster::spPolygons(cds1, cds2)
#' 
#' #plot(r)
#' #plot(polys, add=TRUE)
#' 
#' # Standard raster::extract
#' v <- raster::extract(r, polys)
#' v <- raster::extract(r, polys, fun = mean, df = TRUE)
#' 
#' # Parallel extract
#' #  Register the cluster first (use doMC on linux)
#' cl <- snow::makeCluster(2) #change the 2 to your number of CPU cores  
#' doSNOW::registerDoSNOW(cl)  

#' v2 <- parallelExtract(r, polys, fun = mean)
#' snow::stopCluster(cl)  

parallelExtract <- function(raster, shape, fun = mean, id = 'OBJECTID',  ...){
  
  if (!requireNamespace("foreach", quietly = TRUE)) {
    stop("foreach needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  shape@data[, id] <- as.character(shape@data[, id])
  
  i <- NULL
  # Run extract in parallel.
  values <- foreach::foreach(i = seq_along(shape)) %dopar% {  
    raster::extract(raster, shape[i, ], fun = fun, na.rm = TRUE, cellnumbers = TRUE, ...)
  } 
  if(!is.null(fun)){
    # If a summary function was given, just bind everything together and add ID column
    df <- data.frame(do.call(rbind, values))
    if(class(shape) == 'SpatialPolygonsDataFrame'){
      df <- cbind(ID = as.data.frame(shape)[, id], df)
    } else{
      df <- cbind(ID = names(shape), df)
      id <- 'id'
    }
    
    names(df) <- c(id, names(raster))
    
    return(df)
  } else {
    # If no summary was given we get a list of length n.shapes
    #   each entry in the list is a dataframe with n.covariates columns
    #   Want to make covariates columns, rbind shapes, and add shape and cell id columns.
    
    # list of vectors, one for each covariate
    values_id <- lapply(seq_along(values), function(x) cbind(as.numeric(shape$area_id[x]), values[[x]][[1]]))
    
    
    df <- data.frame(do.call(rbind, values_id))
    names(df) <- c(id, 'cellid', names(raster))
    
    return(df)
  }
  
}



#extract(raster, shape, df = TRUE, fun = mean, na.rm = TRUE)

