#' file.merge - ...
#'
#' @usage file.merge(x, ...)
#'
#' @param x List of data frames
#' @param remove.duplicates Do you want to remove duplicates?
#'
#' @return A combined dataframe.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}. Helpful examples at \url{https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/}
#'
#' @examples
#' file.merge()
#'
#' @export

file.merge <- function(x,remove.duplicates = TRUE)

{
  cell.dat <- plyr::rbind.fill(x)

  if(remove.duplicates == TRUE){
    cell.dat <- cell.dat[!duplicated(cell.dat), ]  # remove rows containing duplicate values within rounding
  }

  assign("cell.dat", cell.dat, envir = globalenv())
  ifelse(exists("cell.dat"), "All files merged into 'cell.dat'", "ERROR in merging files into 'cell.dat'")
}
