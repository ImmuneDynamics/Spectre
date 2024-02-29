#' do.merge.files
#' 
#' Function to merge a list of data.tables (one data.table per 'sample') 
#' into a single large data.table.
#'
#' @usage do.merge.files(dat)
#'
#' @param dat NO DEFAULT. List of data.tables (or data.frames)
#' @param remove.duplicates DEFAULT = TRUE. Do you want to remove duplicated cells?
#'
#' @return A combined data.table.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}.
#'
#' @examples
#' cell.dat <- do.merge.files(dat = data.list, remove.duplicates = TRUE)
#'
#' @export

do.merge.files <- function(dat, remove.duplicates = TRUE) {

    ## Data table (fastest) row binding solution
    cell.dat <- data.table::rbindlist(dat, fill = TRUE)
    
    ## Duplicates
    if (remove.duplicates == TRUE) {
      # cell.dat <- cell.dat[!duplicated(cell.dat), ]  # remove rows containing duplicate values within rounding
      cell.dat <- unique(cell.dat) # removes duplicate rows
    }
    
    
    return(cell.dat)
}
