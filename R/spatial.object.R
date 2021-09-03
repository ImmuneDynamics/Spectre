#' spatial object
#' 
#' @export

### Spatial object

spatial <- setClass(Class = 'spatial', 
                    slots = c(RASTERS = 'RasterStack', 
                              MASKS = 'list',
                              DATA = 'list'
                    ))
