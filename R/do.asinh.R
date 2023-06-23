#' Run ArcSinh transformation
#'
#' Transform data in selected columns using ArcSinh transformation with a specified co-factor.
#'
#' @seealso \url{https://sydneycytometry.org.au/spectre} for usage instructions and vignettes.
#' @references \url{https://sydneycytometry.org.au/spectre}
#'
#' @param dat NO DEFAULT. data.table input sample.
#' @param use.cols NO DEFAULT. Vector of character column names -- these columns will be transformed and added to the data.table as new columns.
#' @param cofactor DEFAULT = 5. Co-factor to use for arcsinh transformation.
#' @param append.cf DEFAULT = FALSE. Appends the co-factor used to the end of the name of the transformed columns.
#' @param reduce.noise DEFAULT = FALSE. This is an experimental calculation which should reduce noise from negative values. Use with caution.
#' @param digits DEFAULT = NULL. Number of decimal places as a limit, not used if NULL. Values beyond will be rounded. Equal to the number or less (i.e. if 9 is used, but only 5 digits are present, then 5 digits will be used). Important to control for small floating point error differences in different OS.
#'
#' @return A data.table with new columns added, that contain the asinh transformed data.
#'
#' @usage do.asinh(dat, use.cols)
#'
#' @import data.table
#'
#' @export do.asinh

do.asinh <- function(dat,
                     use.cols,
                     cofactor = 5,
                     append.cf = FALSE,
                     reduce.noise = FALSE,
                     digits = NULL) {

  ### Setup data
  value <- dat[, use.cols, with = FALSE]

  ### Numeric checks

  if (isFALSE(all(sapply(value, is.numeric)))) {
    message("It appears that one column in your dataset is non numeric")
    print(sapply(value, is.numeric))
    stop("do.asinh stopped")
  }

  ### Optional noise reduction

  # https://github.com/JinmiaoChenLab/cytofkit/issues/71
  if (reduce.noise == TRUE) {
    message("This noise reduction function is experimental, and should be used with caution")
    value <- value - 1
    loID <- which(value < 0)
    if (length(loID) > 0) {
      value[loID] <- rnorm(length(loID), mean = 0, sd = 0.01)
    }
  }

  ### Arcsinh calculation

  value <- as.matrix(value)
  value <- value / cofactor
  value <- asinh(value) # value <- log(value + sqrt(value^2 + 1))
  
  if(!is.null(digits)){
    value <- round(value, digits = digits)
  }
  
  value <- as.data.table(value)
  
  ### Options to append the CF used

  if (append.cf == TRUE) {
    if (length(use.cols) > 1) {
      names(value) <- paste0(names(value), "_asinh_cf", cofactor)
    }
    if (length(use.cols) == 1) {
      names(value) <- paste0(use.cols, "_asinh_cf", cofactor)
    }
  }

  if (append.cf == FALSE) {
    if (length(use.cols) > 1) {
      names(value) <- paste0(names(value), "_asinh")
    }
    if (length(use.cols) == 1) {
      names(value) <- paste0(use.cols, "_asinh")
    }
  }

  ### Wrap up
  dat <- cbind(dat, value)
  return(dat)
}
