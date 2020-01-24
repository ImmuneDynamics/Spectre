#' file.merge - ...
#'
#' @usage file.merge(x, ...)
#'
#' @param x List of data frames or data tables
#' @param remove.duplicates Do you want to remove duplicates?
#'
#' @return A combined data.table
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}.
#'
#' @examples
#' file.merge()
#'
#' @export

file.merge <- function(x,remove.duplicates = TRUE)

{

  ## Data table (fastest) row binding solution
    cell.dat <- data.table::rbindlist(x, fill = TRUE)

  ## Plyr (slower) row binding solution (then need to convert from DF to DT)
    # cell.dat <- plyr::rbind.fill(x) # slower plyr solution
    # cell.dat <- as.data.table(cell.dat)

  ## Duplicates
  if(remove.duplicates == TRUE){
    cell.dat <- cell.dat[!duplicated(cell.dat), ]  # remove rows containing duplicate values within rounding
  }

  ## Assign to global environment
  #assign("cell.dat", cell.dat, envir = globalenv())
  return(cell.dat)

  #ifelse(exists("cell.dat"), "All files merged into 'cell.dat'", "ERROR in merging files into 'cell.dat'")
  ifelse(exists("cell.dat"), message("All files merged into 'cell.dat'"), stop("ERROR in merging files into 'cell.dat'"))

  }
