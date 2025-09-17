#' Transform data using ArcSinh transformation
#'
#' Function to transform data in selected columns using ArcSinh transformation
#' with a specified co-factor.
#'
#' @param dat NO DEFAULT. data.table object. Input data.
#' @param use.cols NO DEFAULT. Vector of character column names.
#' The columns to transform.
#' These columns will be transformed and added to the data.table as new columns.
#' @param cofactor DEFAULT = 5. Co-factor to use for arcsinh transformation.
#' Can be list of co-factors that matches the columns in `use.cols`.
#' If specifying a list of co-factors make sure each column in `use.cols` is assigned 
#' a co-factor (repeat values if necessary).
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
    value <- as.matrix(value)
    if (length(cofactor) == 1) {
        value <- value / cofactor
    } else {
        # Check number of co-factors is equal to number of columns
        if (length(cofactor) != length(use.cols)) {
            message("Number of cofactors is not equal to number of columns. Repeat cofactor values if needed or only input one cofactor for all columns.")
            stop("do.asinh stopped")
        }
        # Divide each column by the relevant optimised cofactor
        # https://stackoverflow.com/questions/48151278/dividing-columns-of-a-matrix-by-elements-of-a-vector
        value <- sweep(value, 2, cofactor, FUN = "/")
    }
    # Calculate asinh
    value <- asinh(value)
    value <- data.table::as.data.table(value)

    ### Options to append the CF used
    dat_out <- data.table::copy(dat)

    for (i in seq_along(use.cols)) {
        marker <- use.cols[i]
        new_marker_name <- paste0(marker, "_asinh")
        if (append.cf) {
            if (length(cofactor) == 1) {
                new_marker_name <- paste0(marker, "_asinh_cf", cofactor)
            } else {
                new_marker_name <- paste0(marker, "_asinh_cf", cofactor[i])
            }
        }
        dat_out <- .add_or_replace_column(
            dt = dat_out,
            col_name = new_marker_name,
            values = value[[marker]]
        )
    }

    return(dat_out)
}
