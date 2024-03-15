#' @param input_data_name Character. The name of the data in Spectre object to 
#' apply arc-sinh transformation to.
#' Only used if dat is a Spectre object.
#' @param output_data_name Character. What name should the arc-sinh transformed 
#' data be stored under in the Spectre object.
#' Only used if dat is a Spectre object.
#'
#' @exportMethod do.asinh
#' @rdname do.asinh
setMethod("do.asinh", "Spectre", function(
        dat,
        use.cols,
        input_data_name,
        output_data_name,
        cofactor = 5,
        append.cf = FALSE,
        reduce.noise = FALSE,
        digits = NULL,
        verbose = TRUE) {
    
    dat_to_apply_asinh <- dat[[input_data_name]]
    
    if (verbose) {
        message(paste(
            "Performing arcsinh transformation to", input_data_name, "containing",
            nrow(dat_to_apply_asinh), "cells and",
            ncol(dat_to_apply_asinh), "features"))
    }
        
    asinh_dat <- do_actual_transformation(
        dat = dat_to_apply_asinh,
        use.cols = use.cols,
        cofactor = cofactor,
        append.cf = append.cf,
        reduce.noise = reduce.noise,
        digits = digits
    )
    
    # just so the cell id column is first!
    cell_id_col <- dat@cell_id_col
    asinh_dat <- data.table(cell_id = dat_to_apply_asinh[[cell_id_col]], asinh_dat)
    setnames(asinh_dat, "cell_id", cell_id_col)
    
    dat <- add.new.data(dat, asinh_dat, output_data_name)
    
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
    }
        
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
    
    if (!is.data.table(cofactor)) {
        value <- value[, lapply(.SD, function(x) asinh(x / cofactor)), .SDcols = use.cols]
        ### Options to append the CF used
        if (append.cf)
            names(value) <- paste0(names(value), "_asinh_cf", cofactor)
        else
            names(value) <- paste0(names(value), "_asinh")
    } else {
        # Apply different co-factor to different markers
        
        # R you jerk. Why can't you differentiate a vector from numeric!?
        cofactor_array <- cofactor$cofactor
        names(cofactor_array) <- cofactor$marker
        
        # check the use.cols and cofactor has exactly the same element.
        in_cofactor_not_in_usecols <- setdiff(names(cofactor_array), use.cols)
        
        if (length(in_cofactor_not_in_usecols) > 0) {
            error(paste(
                "These columns have co-factors specified but are not present in use.cols:",
                paste0(in_cofactor_not_in_usecols, collapse = ", "),
                "do.asinh STOP! Please make sure the markers use.cols and the marker column in cofactor data.table are consistent!"
            ))
        }
        
        in_usecols_not_in_cofactor <- setdiff(use.cols, names(cofactor_array))
        
        if (length(in_usecols_not_in_cofactor) > 0) {
            error(paste(
                "These columns are present in use.cols but missing co-factor:",
                paste0(in_usecols_not_in_cofactor, collapse = ", "),
                "do.asinh STOP! Please make sure the markers use.cols and the marker column in cofactor data.table are consistent!"
            ))
        }
        
        
        
        # Do the actual asinh transformation..
        value <- lapply(names(cofactor_array), function(marker) {
            asinh(value[[marker]] / cofactor_array[marker])
        })
        
        ### Options to append the CF used
        if (append.cf)
            names(value) <- paste0(names(cofactor_array), "_asinh_cf", cofactor_array)
        else
            names(value) <- paste0(names(cofactor_array), "_asinh")
        
        value <- do.call(cbind, value)
        value <- data.table(value)
    }
    
    
    if(!is.null(digits)){
        value <- round(value, digits = digits)
    }
    
    return(value)
    
}







