#' ilastik.prep
#'
#' @usage ilastik.prep(from, to, files)
#'
#' @param from NO DEFAULT. Directory containing FOLDERS (ROIs) or OME.TIFF files
#' @param to NO DEFAULT. Directory where you would like to copy selected OME.TIFF files for use in Ilastik
#' @param files NO DEFAULT. A list of image names (e.g. c("148Nd_148Nd_CD11b.ome.tiff", "152Sm_152Sm_CD3e.ome.tiff") indicating which images you wish to use in Ilastik.
#'
#' @return ...
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}.
#'
#' @examples ilastik.prep()
#'
#' @export

ilastik.prep <- function(from, to, files){
  
  ## Test data
      # from <- PrimaryDirectory
      # to <- OutputDirectory
      # files <- for.ilastik
  
  ## Copy folder names from the 'from' to the the 'to' directory
      folders <- list.dirs(from, recursive = FALSE, full.names = FALSE)
      setwd(to)
      for(i in folders){
        dir.create(i, showWarnings = FALSE)
      }
  
  ## Copy the selected files across
      for(a in folders){
        for(i in files){
          file.copy(from = paste0(PrimaryDirectory, "/", a, "/", i), to = paste0(OutputDirectory, "/", a, "/", i))
        }
      }
  
}




########################################################################
# copy.folders <- function(from, to){
#   folders <- list.dirs(from, recursive = FALSE, full.names = FALSE)
#   setwd(to)
#   for(i in folders){
#     dir.create(i, showWarnings = FALSE)
#   }
# }
# 
# copy.files <- function(from, to) {
#   todir <- dirname(to)
#   if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
#   #file.rename(from = from,  to = to)
#   file.copy(from = from, to = to)
# }


