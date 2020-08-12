#' do.zscore - Calculate z.score for a dataset
#'
#' @usage do.zscore(dat, use.cols)
#'
#' @param dat NO DEFAULT. A data.table (or data.frame) containing the data to be converted to z-scores.
#' @param use.cols NO DEFAULT. The columns to be used for z-score calculations.
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
                      use.cols) 

{
  
  ### Testing
      # dat <- demo.sum
      # use.cols <- names(demo.sum)[c(4:15)]
  
  ### Using R's scale()
      res <- dat
  
      sc <- scale(dat[,use.cols, with = FALSE])
      sc <- as.data.table(sc)
      res[,use.cols] <- sc
      
  ### Manual calculation
      # zscore <- function(x) {
      #   z <- (x - mean(x)) / sd(x)
      #   return(z)
      # }
      
  ### Return
      
      return(res)
  
}
