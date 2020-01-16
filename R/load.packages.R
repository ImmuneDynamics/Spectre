#' load.packages
#'
#' This function allows you to load all of the common use packages dependencies for Spectre.
#' 
#' @return loads all the common use package libraries.
#' 
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}
#'
#' @usage See \url{https://sydneycytometry.org.au/spectre} for usage instructions and vignettes.
#'
#' @examples
#' load.packages()
#'
#' @export

load.packages <- function()
{
  require('devtools')
  require('data.table')
  require('plyr')
  require('dplyr')
  require('tidyr')
  require('rstudioapi')
  require('Rtsne')
  require('umap')
  require('ggplot2')
  require('ggthemes')
  require('scales')
  require('colorRamps')
  require('RColorBrewer')
  require('gridExtra')
  require('caret')
  }