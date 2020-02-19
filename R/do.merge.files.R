#' do.merge.files - Function to merge a list of data.tables (one data.table per 'sample') into a single large data.table.
#'
#' @usage file.merge(x, remove.duplicates)
#'
#' @param x NO DEFAULT. List of data.tables (or data.frames)
#' @param remove.duplicates DEFAULT = TRUE. Do you want to remove duplicates?
#'
#' @return Returns a combined data.table.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}.
#'
#' @examples
#' cell.dat <- file.merge(x = data.list, remove.duplicates = TRUE)
#'
#' @export

file.merge <- function(x, remove.duplicates = TRUE){

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
  ifelse(exists("cell.dat"), message("All files merged into a single data.table"), stop("ERROR in merging files into a single data.table"))

}

