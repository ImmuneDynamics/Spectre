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
#' Can be vector of co-factors that align with columns in 'use.cols'
#' See details for more information.
#' It can also be set to NULL to get the function to automatically infer the co-factors.
#' @param cofactor_inference_method DEFAULT flowVS.
#' Whether to automatically calculates co-factor for each marker.
#' Only used if cofactor is NULL, and will force 'append.cf' to TRUE.
#' Accepted options: 'flowVS' or 'top10'.
#' flowVS uses the flowVS::estParamFlowVS (be patient as it takes time). 
#' 'top10' uses the 90th percentile of values that are negative (based on absolute value) 
#' as the co-factor, and if values are all positive it will use the 10th percentile as co-factor. 
#' @param append.cf DEFAULT = FALSE. Appends the co-factor used to the end of 
#' the name of the transformed columns.
#' Regardless of whether this is TRUE or FALSE, the cofactors used for each marker
#' will be added as attributes to the returned data.table.
#' Use `attributes(<returned_data_table>$cofactors)` to check.
#' @param reduce.noise DEFAULT = FALSE. This is an experimental calculation 
#' which should reduce noise from negative values. Use with caution.
#' @param digits DEFAULT = NULL. Number of decimal places as a limit, not used if NULL. 
#' Values beyond will be rounded. 
#' Equal to the number or less (i.e. if 9 is used, but only 5 digits are present, 
#' then 5 digits will be used). Important to control for small floating point 
#' error differences in different OS.
#' @param verbose DEFAULT = TRUE.
#' If TRUE, the function will print progress updates as it executes.
#' @param add_to_table DEFAULT = TRUE 
#' Whether to add the arcsinh-ed columns into the input data (TRUE) or not (FALSE).
#' If FALSE, a new data.table containing the arcsinh-ed columns will be returned.
#' If TRUE and there happen to be some existing columns in your input data that share
#' the same name as the archsinh-ed markers, they will be overwritten.
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
#' dat <- data.table(NK11=rnorm(10, 2), CD3=rnorm(10,1), cell_id=paste0("cell_", seq(10)))
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
#' @return A data.table with new columns added, that contain the asinh transformed data.
#'
#' @import data.table
#' 
#' @usage do.asinh(dat, use.cols, cofactor=5, append.cf=FALSE, reduce.noise=FALSE)
#'
#' @export do.asinh

do.asinh <- function(dat,
                     use.cols,
                     cofactor = 5,
                     cofactor_inference_method = c("flowVS", "top10"),
                     append.cf = FALSE,
                     reduce.noise = FALSE,
                     digits = NULL,
                     add_to_table = TRUE,
                     verbose = TRUE) {
    
    if (verbose) {
        message("Doing some checks on columns in use.cols.")
    }
    
    # Check if the columns exist in the data.table
    columns_exist <- all(use.cols %in% colnames(dat))
    if (!all(columns_exist)) {
        stop(paste(
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
        stop(paste(
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
    
    # Infer co-factor or checking the specified co-factor is acceptable.
    if (is.null(cofactor)) {
        cofactor_inference_method <- match.arg(cofactor_inference_method)
        if (verbose) {
            message(paste(
                "Inferring co-factors using",
                cofactor_inference_method,
                "algorithm to the following markers:",
                paste(use.cols, collapse = ", ")
            ))
            
            message("Please check your plots to make sure the inferred co-factors are appropriate!")
        }
        cofactor <- infer_cofactor(
            value = value, 
            cofactor_inference_method = cofactor_inference_method,
            use.cols = use.cols
        )
        # automatically append the co-factor.
        # append.cf <- TRUE
        
    } else {
        # This is where the co-factor is manually specified
        if (verbose) {
            message("Doing some checks on cofactors.")
        }
        
        # Check that co-factors have been specified correctly.
        if (length(cofactor) > 1 & length(cofactor) != length(use.cols)) {
            stop(paste(
                "You have specified more than one co-factor, but",
                "the number of markers to apply arc-sinh to (", length(use.cols), "markers )",
                "does not match the number of specified co-factors (", length(cofactor), "cofactors )."
            ))
        }
        
        # will only get here if the checks are all green!
        if (verbose) {
            message(paste(
                "Applying co-factor of",
                paste(cofactor, collapse = ", "),
                "to the following markers (in order):",
                paste(use.cols, collapse = ", ")
            ))
        }
        
        # One co-factor rules it all.
        if (length(cofactor) == 1) {
            # value <- value[, lapply(.SD, function(x) asinh(x / cofactor)), .SDcols = use.cols]
            cofactor <- rep(cofactor[1], length(use.cols))
            
        } 
    }
    
    # Do the actual asinh transformation..
    # Divide each column by the relevant optimised cofactor
    #https://stackoverflow.com/questions/48151278/dividing-columns-of-a-matrix-by-elements-of-a-vector
    value <- sweep(value, 2, cofactor, FUN = '/')
    value <- asinh(value)
    
    ### Options to append the CF used
    if (append.cf) {
        names(value) <- paste0(use.cols, "_asinh_cf", cofactor)
    } else {
        names(value) <- paste0(use.cols, "_asinh")
    }
    
    if(!is.null(digits)){
        value <- round(value, digits = digits)
    }
    
    # needed so add_to_data_table works properly
    value <- data.table(value)
    
    names(cofactor) <- use.cols
    
    if (add_to_table) {
        res_to_return <- add_to_data_table(dat, value)
        attr(res_to_return, "cofactors") <- cofactor
        return(res_to_return)
    } else {
        attr(value, "cofactors") <- cofactor
        return(value)
    }
    
}

#' Internal function which infer the co-factor using flowVS or top10 method.
#'
#' @param value a data.table containing only the markers to be transformed.
#' @param cofactor_inference_method which inferrence method to use? flowVS or top10.
#' @param use.cols A vector of character column names to apply arc-sinh transformation to.
#'
#' @author Felix Marsh-Wakefield
#' 
infer_cofactor <- function(value, cofactor_inference_method, use.cols) {
    if (cofactor_inference_method == 'flowVS') {
        check_packages_installed(c("flowVS", "Biobase"))
        
        ## Create flowFrame metadata (column names with descriptions) plus flowFrame
        metadata <- data.frame(name = dimnames(value)[[2]], desc = paste("column", dimnames(value)[[2]], "from dataset"))
        dat.ff <- new("flowFrame",
                      exprs = as.matrix(value), # in order to create a flow frame, data needs to be read as matrix
                      parameters = Biobase::AnnotatedDataFrame(metadata)
        )
        
        # Needs to be a 'flowSet' object...
        dat.ff <- flowCore::flowSet(dat.ff)
        
        ## Calculate cofactors
        cofactors <- flowVS::estParamFlowVS(dat.ff, channels = use.cols)
    } else if (cofactor_inference_method == 'top10') {
        ## Calculate cofactors
        cofactors <- lapply(value, function(x) {
            if(min(x)<0) {
                abs.dat <- abs(x[x<0])
                
                quantile(abs.dat,
                         probs = 0.9,
                         na.rm = TRUE)
                
            } else {
                # median(x)
                quantile(x,
                         probs = 0.1,
                         na.rm = TRUE)
            }
        })
        
        cofactors <- unlist(cofactors)
    }
    return(cofactors)
    
}







