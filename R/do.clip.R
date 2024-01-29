#' Clip data
#'
#' Clips data using a specified lower and upper value.
#' E.g. if the upper value is set to 1000, then any values above 1000 will be converted to 1000. No rows are lost, rather the values in those rows are converted.
#'
#' @seealso \url{https://sydneycytometry.org.au/spectre} for usage instructions and vignettes.
#' @references \url{https://sydneycytometry.org.au/spectre}
#'
#' @param dat NO DEFAULT. data.table Input sample.
#' @param use.cols NO DEFAULT. Columns to clip
#' @param min.value DEFAULT = NULL. Specify the value for which values below will be converted to. If set to NULL, no clipping will occur.
#' @param max.value DEFAULT = NULL. Specify the value for which values above will be converted to. If set to NULL, no clipping will occur.
#' @param append.name DEFAULT = '_clipped'. Text appended to new column names after clipping.
#'
#' @return A data.table with new columns added, that contain the clipped data.
#'
#' @usage do.clip(dat, use.cols)
#'
#' @import data.table
#'
#' @export do.clip

do.clip <- function(dat,
                    use.cols,
                    min.value = NULL,
                    max.value = NULL,
                    append.name = "_clipped") {

  ### Packages

  # Require: data.table

  ### Demo data

  # dat <- demo.asinh
  # use.cols <- names(demo.asinh)[c(11:19)]
  # min.value <- 1
  # max.value <- 3
  # append.name = '_clipped'

  ### Setup data

  value <- dat[, use.cols, with = FALSE]

  ### Numeric checks

  if (isFALSE(all(sapply(value, is.numeric)))) {
    message("It appears that one column in your dataset is non numeric")
    print(sapply(value, is.numeric))
    stop("do.clip stopped")
  }

  ### Lower clipping

  if (!is.null(min.value)) {
    for (a in use.cols) {
      # a <- use.cols[[1]]

      temp <- value[, a, with = FALSE]
      temp[temp[[a]] < min.value, ] <- min.value
      value[, a] <- temp
    }
  }

  ### Upper clipping

  if (!is.null(max.value)) {
    for (a in use.cols) {
      # a <- use.cols[[1]]

      temp <- value[, a, with = FALSE]
      temp[temp[[a]] > max.value, ] <- max.value
      value[, a] <- temp
    }
  }

  ### Wrap up

  names(value) <- paste0(names(value), append.name)

  dat <- cbind(dat, value)
  return(dat)
}
