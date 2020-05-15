#' do.split.data - Function to split a data.tables into multiple files based on time point.
#'
#' @usage do.split.data(cell.dat. time.point.col)
#'
#' @param cell.dat NO DEFAULT. List of data.tables (or data.frames)
#' @param time.point.col NO DEFAULT. Column denoting the time point.
#' @param clustering.cols NO DEFAULT. The columns to be used for clustering.
#'
#'
#' @author Givanna Putri, \email{ghar1821@@uni.sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}.
#'
#' @examples
#' do.split.data(cell.dat = cell.dat, time.point.col = 'TimePoint', clustering.cols = c("SCA1", "FITC Ly6C"))
#'
#' @export

do.split.data <- function(cell.dat, time.point.col, clustering.cols) {
  if(!is.element('Spectre', installed.packages()[,1])) stop('Spectre is required but not installed')
  require(Spectre)
  
  time.points <- unique(cell.dat[[time.point.col]])
  sapply(time.points, function(time.point) {
    cell.dat.subset <- cell.dat[cell.dat[[time.point.col]] == time.point, clustering.cols, with=FALSE]
    Spectre::write.files(cell.dat.subset, time.point)
  })
} 
