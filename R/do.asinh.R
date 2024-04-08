#' Do ArcSinh transformation
#'
#' Transform data in selected columns using ArcSinh transformation with a specified co-factor.
#'
#' @seealso \url{https://immunedynamics.io/spectre/} for usage instructions and vignettes.
#' @references \url{https://immunedynamics.io/spectre/}
#'
#' @param dat NO DEFAULT. 
#' Either a data.table or a Spectre object to apply arc-sinh transformation to.
#' @param use.cols NO DEFAULT. 
#' A vector of character column names to apply arc-sinh transformation to.
#' @param cofactor DEFAULT = 5. Co-factor to use for arcsinh transformation.
#' Can be vector of co-factors that align with columns in 'use.cols'.
#' See details for more information.
#' @param append.cf DEFAULT = FALSE. Appends the co-factor used to the end of 
#' the name of the transformed columns.
#' @param reduce.noise DEFAULT = FALSE. This is an experimental calculation 
#' which should reduce noise from negative values. Use with caution.
#' @param digits DEFAULT = NULL. Number of decimal places as a limit, not used if NULL. 
#' Values beyond will be rounded. 
#' Equal to the number or less (i.e. if 9 is used, but only 5 digits are present, 
#' then 5 digits will be used). Important to control for small floating point 
#' error differences in different OS.
#' @param verbose DEFAULT = TRUE.
#' If TRUE, the function will print progress updates as it executes.
#' 
#' @details
#' ## Specifying different cofactor for different marker
#' 
#' If use.cols is `c("CD3", "CD4")`, and cofactor is `c(5,6)`, co-factor
#' 5 will be applied to CD3 and 6 will be applied to CD4.
#' 
#' In the case where there is more than one cofactor specified, but there are more or less
#' cofactors than the number of markers in use.cols, the function will give you an error. 
#' Hence, if you intend to specify different co-factor for different marker,
#' make sure each marker is assigned a co-factor (repeat values if necessary),
#' 
#' An exception to this is when you have multiple markers in use.cols but you only specify
#' one cofactor as either a numerical value or a vector of size one, 
#' e.g., `cofactor = 5` or `cofactor = c(5)`.
#' In this case, the function will apply a cofactor of 5 to all the markers in use.cols.
#' Unfortunately there is no way in R to differentiate a numerical variable and a vector.
#' 
#' 
#' @examples
#' library(data.table)
#' # Assuming dat is a data.table
#' dat <- data.table(NK11=rnorm(10, 2), CD3=rnorm(10,1))
#' # Default co-factor 5 will be used for both NK11 and CD3
#' dat_asinh <- do.asinh(dat, use.cols=c("NK11", "CD3"))
#' dat_asinh
#' 
#' # Apply asinh with default co-factor to only NK11
#' dat_asinh <- do.asinh(dat, use.cols=c("NK11"))
#' dat_asinh
#' 
#' # Apply different co-factor to the markers
#' dat_asinh <- do.asinh(dat, use.cols=c("NK11", "CD3"), cofactor=c(5,10))
#' dat_asinh
#' 
#' 
#' @export
setGeneric("do.asinh", function(dat,
                                use.cols,
                                cofactor = 5,
                                append.cf = FALSE,
                                reduce.noise = FALSE,
                                digits = NULL,
                                verbose = TRUE,
                                ...) {
    standardGeneric("do.asinh")
    
})

#' @param data_source Character. The name of the data in Spectre object to 
#' apply arc-sinh transformation to.
#' Only used if dat is a Spectre object.
#' @param output_name Character. What name should the arc-sinh transformed 
#' data be stored under in the Spectre object.
#' Only used if dat is a Spectre object.
#'
#' @exportMethod do.asinh
#' @rdname do.asinh
setMethod("do.asinh", "Spectre", function(
        dat,
        use.cols,
        data_source,
        output_name,
        cofactor = 5,
        append.cf = FALSE,
        reduce.noise = FALSE,
        digits = NULL,
        verbose = TRUE) {
    
    dat_to_apply_asinh <- dat[[data_source]]
    
    if (verbose) {
        message(paste(
            "Performing arcsinh transformation to", data_source, "containing",
            nrow(dat_to_apply_asinh), "cells and",
            ncol(dat_to_apply_asinh), "features"))
        
        message(paste(
            "Applying co-factor of",
            paste(cofactor, collapse = ", "),
            "to the following markers (in order):",
            paste(use.cols, collapse = ", ")
        ))
    }
        
    asinh_dat <- do_actual_transformation(
        dat = dat_to_apply_asinh,
        use.cols = use.cols,
        cofactor = cofactor,
        append.cf = append.cf,
        reduce.noise = reduce.noise,
        digits = digits,
        verbose = verbose
    )
    
    # just so the cell id column is first!
    cell_id_col <- dat@cell_id_col
    asinh_dat <- data.table(cell_id = dat_to_apply_asinh[[cell_id_col]], asinh_dat)
    setnames(asinh_dat, "cell_id", cell_id_col)
    
    dat <- add.new.data(dat, asinh_dat, output_name)
    
    return(dat)
    
})

#' @rdname do.asinh
#' @exportMethod do.asinh
setMethod("do.asinh", "data.table", function(
        dat,
        use.cols,
        cofactor = 5,
        append.cf = FALSE,
        reduce.noise = FALSE,
        digits = NULL,
        verbose = TRUE) {
    
    if (verbose) {
        message(paste("Performing arcsinh transformation for data.table containing",
                      nrow(dat), "cells and",
                      ncol(dat), "features"))
        
        message(paste(
            "Applying co-factor of",
            paste(cofactor, collapse = ", "),
            "to the following markers (in order):",
            paste(use.cols, collapse = ", ")
        ))
    }
        
    asinh_dat <- do_actual_transformation(
        dat = dat,
        use.cols = use.cols,
        cofactor = cofactor,
        append.cf = append.cf,
        reduce.noise = reduce.noise,
        digits = digits,
        verbose = verbose
    )
    dat <- cbind(dat, asinh_dat)
    
    return(dat)
    
})


# Internal function that actually do the asinh transformation
do_actual_transformation <- function(dat,
                                     use.cols,
                                     cofactor = 5,
                                     append.cf = FALSE,
                                     reduce.noise = FALSE,
                                     digits = NULL,
                                     verbose = TRUE) {
    
    if (verbose) {
        message("Doing some checks on columns in use.cols.")
    }
    
    # Check if the columns exist in the data.table
    columns_exist <- all(use.cols %in% colnames(dat))
    if (!all(columns_exist)) {
        error(paste(
            "Some values in use.cols do not exist as columns in dat!",
            paste0("Columns: ", paste(use.cols, collapse = ", ")),
            paste0("Exist? ", paste(columns_exist, collapse = ", ")),
            "do.asinh STOP!",
            sep = "\n"
        ))
    }
    
    
    ### Numeric columns checks
    
    are_columns_numeric <- sapply(use.cols, function(col) is.numeric(dat[[col]]))
    
    if (!all(are_columns_numeric)) {
        error(paste(
            "It appears that some columns in the dataset are not numeric!",
            paste0("Columns: ", paste(use.cols, collapse = ", ")),
            paste0("Is numeric? ", paste(are_columns_numeric, collapse = ", ")),
            "do.asinh STOP!",
            sep = "\n"
        ))
    }
    
    if (verbose) {
        message("Doing some checks on cofactors.")
    }
    
    # Check that co-factors have been specified correctly.
    if (length(cofactor) > 1 & length(cofactor) != length(use.cols)) {
        error(paste(
            "You have specified more than one co-factor, but",
            "the number of markers to apply arc-sinh to (", length(use.cols), "markers )",
            "does not match the number of specified co-factors (", length(cofactor), "cofactors )."
        ))
    }
    
    ### Setup data
    value <- dat[, use.cols, with = FALSE]
    
    ### Optional noise reduction
    
    # https://github.com/JinmiaoChenLab/cytofkit/issues/71
    if (reduce.noise) {
        warning("This noise reduction function is experimental, and should be used with caution!!")
        value <- value - 1
        loID <- which(value < 0)
        if (length(loID) > 0) {
            value[loID] <- rnorm(length(loID), mean = 0, sd = 0.01)
        }
    }
    
    ### Arcsinh calculation
    
    if (length(cofactor) == 1) {
        value <- value[, lapply(.SD, function(x) asinh(x / cofactor)), .SDcols = use.cols]
        ### Options to append the CF used
        if (append.cf)
            names(value) <- paste0(names(value), "_asinh_cf", cofactor)
        else
            names(value) <- paste0(names(value), "_asinh")
    } else {
        
        # Apply different co-factor to different markers
        
        # Do the actual asinh transformation..
        # Divide each column by the relevant optimised cofactor
        #https://stackoverflow.com/questions/48151278/dividing-columns-of-a-matrix-by-elements-of-a-vector
        value <- sweep(value, 2, cofactor, FUN = '/')
        value <- asinh(value)
        
        ### Options to append the CF used
        if (append.cf)
            names(value) <- paste0(use.cols, "_asinh_cf", cofactor)
        else
            names(value) <- paste0(use.cols, "_asinh")
        
        value <- do.call(cbind, value)
        value <- data.table(value)
    }
    
    
    if(!is.null(digits)){
        value <- round(value, digits = digits)
    }
    
    return(value)
    
}







