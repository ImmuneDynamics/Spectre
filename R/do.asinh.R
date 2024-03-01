#' @export
#' @rdname SpectreObject-class
setMethod("do.asinh", "SpectreObject", function(
        dat,
        use.cols,
        slot_name = NULL,
        cofactor = 5,
        append.cf = FALSE,
        reduce.noise = FALSE,
        digits = NULL,
        verbose = TRUE) {
    
    if (is.null(slot_name))
        stop("slot is NULL. Don't know which slot to apply arcsinh transformation to.")
    
    dat_to_apply_asinh <- slot(dat, slot_name)
    
    if (verbose)
        message(paste(
            "Performing arcsinh transformation for slot", slot, "containing",
            nrow(dat_to_apply_asinh), "cells and",
            ncol(dat_to_apply_asinh), "features"))
    asinh_dat <- do_actual_transformation(
        dat = dat_to_apply_asinh,
        use.cols = use.cols,
        cofactor = cofactor,
        append.cf = append.cf,
        reduce.noise = reduce.noise,
        digits = digits
    )
    dat_to_apply_asinh <- cbind(dat_to_apply_asinh, asinh_dat)
    
    slot(dat, slot_name) <- dat_to_apply_asinh
    
    return(dat)
    
})

#' @export
setMethod("do.asinh", "data.table", function(
        dat,
        use.cols,
        slot_name = NULL,
        cofactor = 5,
        append.cf = FALSE,
        reduce.noise = FALSE,
        digits = NULL,
        verbose = TRUE) {
    
    if (verbose)
        message(paste("Performing arcsinh transformation for data.table containing",
                      nrow(dat), "cells and",
                      ncol(dat), "features"))
    asinh_dat <- do_actual_transformation(
        dat = dat,
        use.cols = use.cols,
        cofactor = cofactor,
        append.cf = append.cf,
        reduce.noise = reduce.noise,
        digits = digits
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
                                     digits = NULL) {
    
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
    
    ### Optional noise reduction
    
    # https://github.com/JinmiaoChenLab/cytofkit/issues/71
    if (reduce.noise == TRUE) {
        warning("This noise reduction function is experimental, and should be used with caution!!")
        value <- value - 1
        loID <- which(value < 0)
        if (length(loID) > 0) {
            value[loID] <- rnorm(length(loID), mean = 0, sd = 0.01)
        }
    }
    
    ### Arcsinh calculation
    
    value <- dat[, lapply(.SD, function(x) asinh(x / cofactor)), .SDcols = use.cols]
    
    if(!is.null(digits)){
        value <- round(value, digits = digits)
    }
    
    ### Options to append the CF used
    if (append.cf)
        names(value) <- paste0(names(value), "_asinh_cf", cofactor)
    else
        names(value) <- paste0(names(value), "_asinh")
    
    return(value)
    
}







