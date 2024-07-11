#' do.zscore - Calculate z.score for a dataset
#'
#' @usage do.zscore(dat, use.cols)
#'
#' @param dat NO DEFAULT. A data.table (or data.frame) containing the data to be converted to z-scores. Z-score transformed values will be added as new columns.
#' @param use.cols NO DEFAULT. The columns to be used for z-score calculations.
#' @param append.name DEFAULT = '_zscore'. Text to be appended to the end of the new z-score transformed columns.
#' @param replace DEFAULT = FALSE. If FALSE, appends new columns to the data.table. If TRUE, replaces the values in the existing columns with the z-score tranformed values, and does not change the column names.
#'
#' @return Returns a new data.table with z-score calculations for each selected column
#'
#' @examples
#' do.zscore(dat = Spectre::demo.clustered, use.cols = c("NK11", "CD4"))
#' 
#' @import data.table
#'
#' @export do.zscore

do.zscore <- function(dat,
                      use.cols,
                      append.name = "_zscore",
                      replace = FALSE) {
  
  ### Testing

  # dat <- demo.sum
  # use.cols <- names(demo.sum)[c(4:15)]

  ### Using R's scale()

  sc <- scale(dat[, use.cols, with = FALSE])
  sc <- as.data.table(sc)

  ### Replace or not

  if (isTRUE(replace)) {
    res <- dat
    res[, use.cols] <- sc
  } else {
    names(sc) <- paste0(names(sc), append.name)
    res <- cbind(dat, sc)
  }

  ### Manual calculation

  # zscore <- function(x) {
  #   z <- (x - mean(x)) / sd(x)
  #   return(z)
  # }

  ### Return

  return(res)
}
