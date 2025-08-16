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
    value <- asinh(value / cofactor)

    ### Options to append the CF used

    dat_out <- data.table::copy(dat)

    for (marker in use.cols) {
        new_marker_name <- paste0(marker, "_asinh")
        if (append.cf) {
            new_marker_name <- paste0(marker, "_asinh_cf", cofactor)
        }
        dat_out <- .add_or_replace_column(
            dt = dat_out,
            col_name = new_marker_name,
            values = value[[marker]]
        )
    }
    
    return(dat_out)
}
