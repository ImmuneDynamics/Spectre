#' spatial object
#' 
#' @export

### Spatial object

library('Spectre')
require('data.table')
require('raster')

spatial <- setClass(Class = 'spatial', 
                    slots = c(RASTERS = 'RasterStack', 
                              MASKS = 'list',
                              DATA = 'list'
                    ))
