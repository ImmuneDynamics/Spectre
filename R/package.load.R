#' package.load - a function to load (library) all required packages.
#'
#' This function allows you to load all of the common use packages dependencies for Spectre.
#'
#' @return loads all the common use package libraries.
#'
#' @param type DEFAULT = "general". If "general", then loads packages required for general Spectre usage. If "spatial", then loads  additional packages required for spatial analysis. If "ML", then loads additional packages required for machine-learing functionality.
#'
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}
#'
#' @usage See \url{https://sydneycytometry.org.au/spectre} for usage instructions and vignettes.
#'
#' @examples
#' package.load()
#'
#' @export

package.load <- function(type = "general")
{
  if(type == "general"){
    require('devtools')
    require('data.table')
    require('plyr')
    require('dplyr')
    require('tidyr')
    require('rstudioapi')
    require('Rtsne')
    require('umap')
    require('reticulate')
    require('ggplot2')
    require('ggthemes')
    require('scales')
    require('colorRamps')
    require('RColorBrewer')
    require('gridExtra')
    require('ggpointdensity')

    require('flowCore')
    require('Biobase')
    require('flowViz')
    require('FlowSOM')
  }

  if(type == "spatial"){
    require('raster')
    require('tiff')
    require('rgeos')
  }

  if(type == "ML"){
    require('caret')
  }

}
