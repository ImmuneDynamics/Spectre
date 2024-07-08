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
#' @references \url{https://sydneycytometry.org.au/spectre}.
#'
#' @examples
#' data.list <- list(
#' "01_Mock_01" = Spectre::demo.clustered[Sample == '01_Mock_01'],
#' "02_Mock_02" = Spectre::demo.clustered[Sample == '02_Mock_02']
#' )
#' cell.dat <- do.merge.files(dat = data.list, remove.duplicates = TRUE)
#' head(cell.dat)
#' 
#' @usage do.merge.files(dat, remove.duplicates = TRUE)
#'
#' @export do.merge.files

do.merge.files <- function(dat, remove.duplicates = TRUE){
  
  ## Check that necessary packages are installed
    if(!is.element('Spectre', installed.packages()[,1])) stop('Spectre is required but not installed')
    if(!is.element('data.table', installed.packages()[,1])) stop('data.table is required but not installed')

  ## Require packages
    
    require(data.table)

  ## Data table (fastest) row binding solution
    cell.dat <- data.table::rbindlist(dat, fill = TRUE)

  ## Plyr (slower) row binding solution (then need to convert from DF to DT)
    # cell.dat <- plyr::rbind.fill(dat) # slower plyr solution
    # cell.dat <- as.data.table(cell.dat)

  ## Duplicates
  if(remove.duplicates == TRUE){
    # cell.dat <- cell.dat[!duplicated(cell.dat), ]  # remove rows containing duplicate values within rounding
    cell.dat <- unique(cell.dat) # removes duplicate rows
  }

  ## Assign to global environment
  #assign("cell.dat", cell.dat, envir = globalenv())
  return(cell.dat)

  #ifelse(exists("cell.dat"), "All files merged into 'cell.dat'", "ERROR in merging files into 'cell.dat'")
  ifelse(exists("cell.dat"), message("All files merged into a single data.table"), stop("ERROR in merging files into a single data.table"))

}

