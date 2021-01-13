#' do.zscore - Calculate z.score for a dataset
#'
#' @usage do.zscore(dat, use.cols)
#'
#' @param dat NO DEFAULT. A data.table (or data.frame) containing the data to be converted to z-scores. Z-score transformed values will be added as new columns.
#' @param use.cols NO DEFAULT. The columns to be used for z-score calculations.
#' @param append.name DEFAULT = '_zscore'. Text to be appended to the end of the new z-score transformed columns.
#' 
#' @return Returns a new data.table with z-score calculations for each selected column
#'
#' @examples
# res <- do.zscore(dat = Spectre::demo.sum,
#                  use.cols = names(Spectre::demo.sum)[c(4:15)])
#' 
#' @import data.table
#'
#' @export

do.zscore <- function(dat,
                      use.cols,
                      append.name = '_zscore') 

{
  
  ### Packages
  
      require('data.table')
  
  ### Testing
  
      # dat <- demo.sum
      # use.cols <- names(demo.sum)[c(4:15)]
  
  ### Using R's scale()
  
      res <- dat
  
      sc <- scale(dat[,use.cols, with = FALSE])
      sc <- as.data.table(sc)
      names(sc) <- paste0(names(sc), append.name)
      
      res <- cbind(dat, sc)
      
  ### Manual calculation
      
      # zscore <- function(x) {
      #   z <- (x - mean(x)) / sd(x)
      #   return(z)
      # }
      
  ### Return
      
      return(res)
}
