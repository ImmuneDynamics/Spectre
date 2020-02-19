#' check.packages
#'
#' This function allows you to check to see if all the common use packages dependencies for Spectre are installed.
#'
#' @return returns an error message if one of the common use packages are not installed. Proceeds in order of package importance, and only the first error message encountered will be returned.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}
#'
#' @usage See \url{https://sydneycytometry.org.au/spectre} for usage instructions and vignettes.
#'
#' @examples
#' check.packages()
#'
#' @export

check.packages <- function()
{
  if(!is.element('devtools', installed.packages()[,1])) stop('devtools is required by not installed')
  if(!is.element('data.table', installed.packages()[,1])) stop('data.table is required by not installed')
  if(!is.element('plyr', installed.packages()[,1])) stop('plyr is required by not installed')
  if(!is.element('dplyr', installed.packages()[,1])) stop('dplyr is required by not installed')
  if(!is.element('tidyr', installed.packages()[,1])) stop('tidyr is required by not installed')
  if(!is.element('rstudioapi', installed.packages()[,1])) stop('rstudioapi is required by not installed')
  if(!is.element('Rtsne', installed.packages()[,1])) stop('Rtsne is required by not installed')
  if(!is.element('umap', installed.packages()[,1])) stop('umap is required by not installed')
  if(!is.element('ggplot2', installed.packages()[,1])) stop('ggplot2 is required by not installed')
  if(!is.element('ggthemes', installed.packages()[,1])) stop('ggthemes is required by not installed')
  if(!is.element('scales', installed.packages()[,1])) stop('scales is required by not installed')
  if(!is.element('colorRamps', installed.packages()[,1])) stop('colorRamps is required by not installed')
  if(!is.element('RColorBrewer', installed.packages()[,1])) stop('RColorBrewer is required by not installed')
  if(!is.element('gridExtra', installed.packages()[,1])) stop('gridExtra is required by not installed')

  if(!is.element('flowCore', installed.packages()[,1])) stop('gridExtra is required by not installed')
  if(!is.element('Biobase', installed.packages()[,1])) stop('gridExtra is required by not installed')
  if(!is.element('flowViz', installed.packages()[,1])) stop('gridExtra is required by not installed')
  if(!is.element('FlowSOM', installed.packages()[,1])) stop('gridExtra is required by not installed')
  }
