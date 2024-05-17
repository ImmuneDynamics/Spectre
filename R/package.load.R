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

package.load <- function(type = "general"){
  
    require('colorRamps')
    require('data.table')
    require('dendsort')
    require('factoextra')
    require('flowCore')
    require('FlowSOM')
    require('ggplot2')
    require('ggpointdensity')
    require('ggpubr')
    require('ggthemes')
    require('gridExtra')
    require('gtools')
    require('irlba')
    require('parallel')
    require('patchwork')
    require('pheatmap')
    require('RColorBrewer')
    require('rstudioapi')
    require('rsvd')
    require('Rtsne')
    require('scales')
    require('scattermore')
    require('umap')
    require('uwot')
    require('viridis')
  
  if(type == "spatial"){
    require('raster')
    require('tiff')
    require('exactextractr')
    require('sp')
    require('sf')
    require('stars')
    require('qs')
    require('s2')
    require('rhdf5')
    require('HDF5Array')
  }
}
