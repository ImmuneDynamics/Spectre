#' do.rescale - Re-scale data in selected columns between two values, usually 0 and 1.
#'
#' This function allows you to re-scale the data in selected columns between two values, usually 0 and 1.
#'
#' @seealso \url{https://sydneycytometry.org.au/spectre} for usage instructions and vignettes.
#' @references \url{https://sydneycytometry.org.au/spectre}
#'
#' @param dat NO DEFAULT. data.table Input sample.
#' @param use.cols NO DEFAULT. Vector of character column names -- these columns will be normalised and added to the data.table as new columns.
#' @param new.min DEFAULT = 0. The new minimum value.
#' @param new.max DEFAULT = 1. The new maximum value.
#' @param append.name DEFAULT = '_rescaled'. Text to be appended to the column names of re-scaled data.
#' 
#' @return A data.table with new columns added, that contain the re-scaled data.
#'
#' @import data.table
#'
#' @export

do.rescale <- function(dat,
                       use.cols,
                       new.min = 0,
                       new.max = 1,
                       append.name = '_rescaled'
                       ){ 

  ### Packages
      require('data.table')
  
  ### Test data
      # dat <- Spectre::demo.asinh
      # use.cols <- names(dat)[c(11:19)]
      # new.min = -1
      # new.max = 1
      # append.name = '_rescaled'
      
  ### Create normalisation function
      #norm.fun <- function(x) {(x - min(x, na.rm=TRUE))/(max(x,na.rm=TRUE) -min(x, na.rm=TRUE))}
      norm.fun <- function(x) {(x - min(x))/(max(x)-min(x)) * (new.max - new.min) + new.min}
      #norm.fun <- function(x){scales::rescale(value, to=c(new.min,new.max))}

  ### Establish dataset
      value <- dat[,use.cols,with = FALSE]

  ### Normalise between new values
      res <- as.data.table(lapply(value, norm.fun)) # by default, removes the names of each row
      names(res) <- paste0(names(res), append.name)
      dat <- cbind(dat, res)
      
  ### Return
      return(dat)
}
