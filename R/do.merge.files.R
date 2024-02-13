#' do.merge.files - Function to merge a list of data.tables (one data.table per 'sample') into a single large data.table.
#'
#' @usage do.merge.files(dat, remove.duplicates)
#'
#' @param dat NO DEFAULT. List of data.tables (or data.frames)
#' @param remove.duplicates DEFAULT = TRUE. Do you want to remove duplicates?
#'
#' @return Returns a combined data.table.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @examples
#' cell.dat <- do.merge.files(dat = data.list, remove.duplicates = TRUE)
#'
#' @export

do.merge.files <- function(dat, remove.duplicates = TRUE) {
    # require: data.table
  
    ## Data table (fastest) row binding solution
    cell.dat <- data.table::rbindlist(dat, fill = TRUE)
  
    ## Plyr (slower) row binding solution (then need to convert from DF to DT)
    # cell.dat <- plyr::rbind.fill(dat) # slower plyr solution
    # cell.dat <- as.data.table(cell.dat)
  
    ## Duplicates
    if (remove.duplicates == TRUE) {
      # cell.dat <- cell.dat[!duplicated(cell.dat), ]  # remove rows containing duplicate values within rounding
      cell.dat <- unique(cell.dat) # removes duplicate rows
    }
  
    ## Assign to global environment
    # assign("cell.dat", cell.dat, envir = globalenv())
    return(cell.dat)
  
    # ifelse(exists("cell.dat"), "All files merged into 'cell.dat'", "ERROR in merging files into 'cell.dat'")
    # TODO: remove me. Will never get here because of return.
    # ifelse(exists("cell.dat"), message("All files merged into a single data.table"), stop("ERROR in merging files into a single data.table"))
}
