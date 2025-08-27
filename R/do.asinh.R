#' Transform data using ArcSinH transformation
#'
#' Function to transform data in selected columns using ArcSinH transformation
#' with a specified co-factor.
#' 
#'
#' @param dat NO DEFAULT. data.table object. Input data.
#' @param use.cols NO DEFAULT. Vector of character column names.
#' The columns to transform.
#' These columns will be transformed and added to the data.table as new columns.
#' @param cofactor DEFAULT = 5. Co-factor to use for arcsinh transformation.
#' Can be string of co-factors that align with columns in 'use.cols', just make sure each column is assigned a co-factor (repeat values if necessary).
#' @param append.cf DEFAULT = FALSE. 
#' Appends the co-factor used to the end of the name of the transformed columns.
#' @param reduce.noise DEFAULT = FALSE. 
#' This is an experimental calculation which should reduce noise 
#' from negative values. Use with caution.
#'
#' @return A data.table with new columns added, that contain asinh 
#' transformed data.
#'
#' @import data.table
#' 
#' @usage do.asinh(dat, use.cols, cofactor=5, append.cf=FALSE, 
#' reduce.noise=FALSE)
#'
#' @export do.asinh
#' 
do.asinh <- function(dat,
                     use.cols,
                     cofactor = 5,
                     append.cf = FALSE,
                     reduce.noise = FALSE) {
  
  
  value <- dat[, use.cols, with = FALSE]
  
  .check_numeric_columns(value, use.cols)
  
  ### Optional noise reduction
  # https://github.com/JinmiaoChenLab/cytofkit/issues/71
  if (reduce.noise) {
    warning("This noise reduction function is experimental. Use with caution!")
    value <- value - 1
    loID <- which(value < 0)
    if (length(loID) > 0) {
      value[loID] <- rnorm(length(loID), mean = 0, sd = 0.01)
    }
  }
  
  ### Arcsinh calculation
  if(length(cofactor) == 1) {
    value <- as.matrix(value)
    value <- value / cofactor
    value <- asinh(value) # value <- log(value + sqrt(value^2 + 1))
    
    value <- data.table::as.data.table(value)
  } else {
    # Check number of co-factors is equal to number of columns
    if (length(cofactor) != length(use.cols)) {
      message("Number of cofactors is not equal to number of columns. Repeat cofactor values if needed or only input one cofactor for all columns.")
      stop("do.asinh stopped")
    }
    
    value <- as.matrix(value)
    
    # Divide each column by the relevant optimised cofactor
    #https://stackoverflow.com/questions/48151278/dividing-columns-of-a-matrix-by-elements-of-a-vector
    value <- sweep(value, 2, cofactor, FUN = '/')
    
    # Calculate asinh
    value <- asinh(value)
    
    value <- data.table::as.data.table(value)
    
  }
  
  ### Options to append the CF used
  
  if (append.cf == TRUE) {
    if (length(use.cols) > 1) {
      names(value) <- paste0(names(value), "_asinh_cf", cofactor)
    }
    if (length(use.cols) == 1) {
      names(value) <- paste0(use.cols, "_asinh_cf", cofactor)
    }
  }
  
  ### Wrap up
  dat <- cbind(dat, value)
  return(dat)
  
}
