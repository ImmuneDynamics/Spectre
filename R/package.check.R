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

package.check <- function(type = "general")
{
    if(!is.element('data.table', installed.packages()[,1])) message('data.table is required but not installed')
    if(!is.element('plyr', installed.packages()[,1])) message('plyr is required but not installed')
    if(!is.element('dplyr', installed.packages()[,1])) message('dplyr is required but not installed')
    if(!is.element('tidyr', installed.packages()[,1])) message('tidyr is required but not installed')
    if(!is.element('rstudioapi', installed.packages()[,1])) message('rstudioapi is required but not installed')
    if(!is.element('Rtsne', installed.packages()[,1])) message('Rtsne is required but not installed')
    if(!is.element('umap', installed.packages()[,1])) message('umap is required but not installed')
    if(!is.element('reticulate', installed.packages()[,1])) message('reticulate is required but not installed')
    if(!is.element('ggplot2', installed.packages()[,1])) message('ggplot2 is required but not installed')
    if(!is.element('ggthemes', installed.packages()[,1])) message('ggthemes is required but not installed')
    if(!is.element('scales', installed.packages()[,1])) message('scales is required but not installed')
    if(!is.element('colorRamps', installed.packages()[,1])) message('colorRamps is required but not installed')
    if(!is.element('RColorBrewer', installed.packages()[,1])) message('RColorBrewer is required but not installed')
    if(!is.element('gridExtra', installed.packages()[,1])) message('gridExtra is required but not installed')
    if(!is.element('ggpointdensity', installed.packages()[,1])) message('ggpointdensity is required but not installed')
    if(!is.element('pheatmap', installed.packages()[,1])) message('ggpointdensity is required but not installed')
    if(!is.element('ggpubr', installed.packages()[,1])) message('ggpointdensity is required but not installed')

    if(!is.element('flowCore', installed.packages()[,1])) message('flowCore is required but not installed. Please install from BioConductor.')
    if(!is.element('Biobase', installed.packages()[,1])) message('Biobase is required but not installed. Please install from BioConductor.')
    if(!is.element('flowViz', installed.packages()[,1])) message('flowViz is required but not installed. Please install from BioConductor.')
    if(!is.element('FlowSOM', installed.packages()[,1])) message('FlowSOM is required but not installed. Please install from BioConductor.')

  if(type == "spatial"){
    if(!is.element('raster', installed.packages()[,1])) message('raster is required for SPATIAL analysis but is not installed')
    if(!is.element('tiff', installed.packages()[,1])) message('tiff is required for SPATIAL analysis but is not installed')
    if(!is.element('rgeos', installed.packages()[,1])) message('rgeos is required for SPATIAL analysis but is not installed')
    if(!is.element('velox', installed.packages()[,1])) message('velox is required for SPATIAL analysis but is not installed')
    if(!is.element('sp', installed.packages()[,1])) message('sp is required for SPATIAL analysis but is not installed')
    if(!is.element('sf', installed.packages()[,1])) message('sf is required for SPATIAL analysis but is not installed')
    if(!is.element('stars', installed.packages()[,1])) message('stars is required for SPATIAL analysis but is not installed')
  }

  if(type == "ML"){
    if(!is.element('caret', installed.packages()[,1])) message('caret is required for machine learning analysis but is not installed')
    if(!is.element('class', installed.packages()[,1])) message('class is required for machine learning analysis but is not installed')
  }

}
