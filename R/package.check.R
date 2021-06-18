#' package.check - a function to check the installation of all required packages.
#'
#' This function allows you to check to see if all the common use packages dependencies for Spectre are installed.
#'
#' @return returns an error message if one of the common use packages are not installed. Proceeds in order of package importance, and only the first error message encountered will be returned.
#'
#' @param type DEFAULT = "general". If "general", then checks for the packages required for general Spectre usage. If "spatial", then checks for additional packages required for spatial analysis. If "ML", then checks for additional packages required for machine-learing functionality.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}
#'
#' @usage See \url{https://sydneycytometry.org.au/spectre} for usage instructions and vignettes.
#'
#' @examples
#' package.check()
#'
#' @export

package.check <- function(type = "general"){
  
    res.list <- list()
  
    if(!is.element('data.table', installed.packages()[,1])){
      message('data.table is required but not installed')
      res.list[['data.table']] <- 'data.table'
    } 
      
    if(!is.element('plyr', installed.packages()[,1])) {
      message('plyr is required but not installed')
      res.list[['plyr']] <- 'plyr'
    } 
    
    if(!is.element('dplyr', installed.packages()[,1])) {
      message('dplyr is required but not installed')
      res.list[['data.table']] <- 'data.table'
    } 
    
    if(!is.element('tidyr', installed.packages()[,1])) {
      message('tidyr is required but not installed')
      res.list[['tidyr']] <- 'tidyr'
    } 
    
    if(!is.element('rstudioapi', installed.packages()[,1])) {
      message('rstudioapi is required but not installed')
      res.list[['rstudioapi']] <- 'rstudioapi'
    } 
    
    if(!is.element('Rtsne', installed.packages()[,1])) {
      message('Rtsne is required but not installed')
      res.list[['Rtsne']] <- 'Rtsne'
    } 
    
    if(!is.element('umap', installed.packages()[,1])) {
      message('umap is required but not installed')
      res.list[['umap']] <- 'umap'
    } 
    
    if(!is.element('reticulate', installed.packages()[,1])) {
      message('reticulate is required but not installed')
      res.list[['reticulate']] <- 'reticulate'
    } 
    
    if(!is.element('ggplot2', installed.packages()[,1])) {
      message('ggplot2 is required but not installed')
      res.list[['ggplot2']] <- 'ggplot2'
    } 
    
    if(!is.element('ggthemes', installed.packages()[,1])) {
      message('ggthemes is required but not installed')
      res.list[['ggthemes']] <- 'ggthemes'
    } 
    
    if(!is.element('scales', installed.packages()[,1])) {
      message('scales is required but not installed')
      res.list[['scales']] <- 'scales'
    } 
    
    if(!is.element('colorRamps', installed.packages()[,1])) {
      message('colorRamps is required but not installed')
      res.list[['colorRamps']] <- 'colorRamps'
    } 
    
    if(!is.element('RColorBrewer', installed.packages()[,1])) {
      message('RColorBrewer is required but not installed')
      res.list[['RColorBrewer']] <- 'RColorBrewer'
    } 
    
    if(!is.element('gridExtra', installed.packages()[,1])) {
      message('gridExtra is required but not installed')
      res.list[['gridExtra']] <- 'gridExtra'
    } 
    
    if(!is.element('ggpointdensity', installed.packages()[,1])) {
      message('ggpointdensity is required but not installed')
      res.list[['ggpointdensity']] <- 'ggpointdensity'
    } 
    
    if(!is.element('pheatmap', installed.packages()[,1])) {
      message('pheatmap is required but not installed')
      res.list[['pheatmap']] <- 'pheatmap'
    } 
    
    if(!is.element('ggpubr', installed.packages()[,1])) {
      message('ggpubr is required but not installed')
      res.list[['ggpubr']] <- 'ggpubr'
    } 
    
    if(!is.element('caret', installed.packages()[,1])) {
      message('caret is required for machine learning analysis but is not installed')
      res.list[['caret']] <- 'caret'
    } 
    
    if(!is.element('class', installed.packages()[,1])) {
      message('class is required for machine learning analysis but is not installed')
      res.list[['class']] <- 'class'
    } 
    
    if(!is.element('flowCore', installed.packages()[,1])) {
      message('flowCore is required but not installed. Please install from BioConductor.')
      res.list[['flowCore']] <- 'flowCore'
    } 
    
    if(!is.element('Biobase', installed.packages()[,1])) {
      message('Biobase is required but not installed. Please install from BioConductor.')
      res.list[['Biobase']] <- 'Biobase'
    } 
    
    if(!is.element('flowViz', installed.packages()[,1])) {
      message('flowViz is required but not installed. Please install from BioConductor.')
      res.list[['flowViz']] <- 'flowViz'
    } 
    
    if(!is.element('FlowSOM', installed.packages()[,1])) {
      message('FlowSOM is required but not installed. Please install from BioConductor.')
      res.list[['FlowSOM']] <- 'FlowSOM'
    } 
    

  if(type == "spatial"){
    if(!is.element('raster', installed.packages()[,1])) {
      message('raster is required for SPATIAL analysis but is not installed')
      res.list[['raster']] <- 'raster'
    } 
    
    if(!is.element('tiff', installed.packages()[,1])) {
      message('tiff is required for SPATIAL analysis but is not installed')
      res.list[['tiff']] <- 'tiff'
    } 
    
    if(!is.element('rgeos', installed.packages()[,1])) {
      message('rgeos is required for SPATIAL analysis but is not installed')
      res.list[['rgeos']] <- 'rgeos'
    } 
    
    if(!is.element('exactextractr', installed.packages()[,1])) {
      message('exactextractr is required for SPATIAL analysis but is not installed')
      res.list[['exactextractr']] <- 'exactextractr'
    } 
    
    if(!is.element('sp', installed.packages()[,1])) {
      message('sp is required for SPATIAL analysis but is not installed')
      res.list[['sp']] <- 'sp'
    } 
    
    if(!is.element('sf', installed.packages()[,1])) {
      message('sf is required for SPATIAL analysis but is not installed')
      res.list[['sf']] <- 'sf'
    } 
    
    if(!is.element('stars', installed.packages()[,1])) {
      message('stars is required for SPATIAL analysis but is not installed')
      res.list[['stars']] <- 'stars'
    } 
    
    if(!is.element('qs', installed.packages()[,1])) {
      message('qs is required for SPATIAL analysis but is not installed')
      res.list[['qs']] <- 'qs'
    } 
    
    if(!is.element('s2', installed.packages()[,1])) {
      message('s2 is required for SPATIAL analysis but is not installed')
      res.list[['s2']] <- 's2'
    } 
  }
  
    
  ### Final message
    
    if(length(res.list) == 0){
      message("All packages successfully installed")
    }

}
