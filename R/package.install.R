#' package.install - a function to install packages required for Spectre.
#'
#' This function allows you to check to see if all the common use packages dependencies for Spectre are installed. Will only install if the package has not been installed, will not update packages.
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
#' package.install()
#' package.install(type = 'spatial')
#' package.install(type = 'ML')
#'
#' @export

package.install <- function(type = "general")
{
    if(!require('data.table')) {install.packages('data.table')}
    if(!require('plyr')) {install.packages('plyr')}
    if(!require('dplyr')) {install.packages('dplyr')}
    if(!require('tidyr')) {install.packages('tidyr')}
    if(!require('rstudioapi')) {install.packages('rstudioapi')}
    if(!require('Rtsne')) {install.packages('Rtsne')}
    if(!require('umap')) {install.packages('umap')}
    if(!require('ggplot2')) {install.packages('ggplot2')}
    if(!require('ggthemes')) {install.packages('ggthemes')}
    if(!require('scales')) {install.packages('scales')}
    if(!require('colorRamps')) {install.packages('colorRamps')}
    if(!require('RColorBrewer')) {install.packages('')}
    if(!require('gridExtra')) {install.packages('gridExtra')}
    if(!require('ggpointdensity')) {install.packages('ggpointdensity')}
    if(!require('pheatmap')) {install.packages('pheatmap')}
    if(!require('ggpubr')) {install.packages('ggpubr')}
    if(!require('factoextra')) {install.packages('factoextra')}
    if(!require('reticulate')) {install.packages('reticulate')}

    ## Install BiocManager to download packages from Bioconductor
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

    ## Download additional BioConductor packages
    if(!require('flowCore')) {BiocManager::install('flowCore')}
    if(!require('Biobase')) {BiocManager::install('Biobase')}
    if(!require('flowViz')) {BiocManager::install('flowViz')}
    if(!require('FlowSOM')) {BiocManager::install('FlowSOM')}

  if(type == "spatial"){
    if(!require('raster')) {install.packages('raster')}
    if(!require('tiff')) {install.packages('tiff')}
    if(!require('rgeos')) {install.packages('rgeos')}
    if(!require('velox')) {
      require('devtools')
      install_github("hunzikp/velox")
    }
    if(!require('sp')) {install.packages('sp')}
    if(!require('sf')) {install.packages('sf')}
    if(!require('stars')) {install.packages('stars')}
    if(!require('qs')) {install.packages('qs')}
  }

  if(type == "ML"){
    if(!require('caret')) {install.packages('caret')}
    if(!require('class')) {install.packages('class')}
  }

}
