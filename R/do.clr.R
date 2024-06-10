#' Do CLR transformation
#'
#' Perform CLR transformation.
#' 
#' @param dat A spectre object.
#' @param data_source Character. The name of the data in Spectre object to 
#' apply CLR transformation to.
#' @param output_name Character. What name should the CLR transformed 
#' data be stored under in the Spectre object.
#' @param use_cols A vector of character column names to apply CLR transformation to.
#' @param how Whether to apply the CLR transformation across cells ("across_cells")
#' which essentially apply CLR for each feature, or across features ("across_features")
#' which essentially apply CLR for each cell.
#' @param verbose DEFAULT = TRUE.
#' If TRUE, the function will print progress updates as it executes.
#'
#' @return An updated spectre object
#' @export
#'
do.clr <- function(dat, data_source, output_name, 
                   use_cols, how = c("across_cells", "across_features"),
                   verbose = FALSE) {
    
    # TODO need to add check whether dat is a spectre object
    
    dat_to_apply_clr <- dat[[data_source]]
    
    are_columns_numeric <- sapply(use_cols, function(col) is.numeric(dat_to_apply_clr[[col]]))
    
    if (!all(are_columns_numeric)) {
        stop(paste(
            "It appears that some columns in the dataset are not numeric!",
            paste0("Columns: ", paste(use.cols, collapse = ", ")),
            paste0("Is numeric? ", paste(are_columns_numeric, collapse = ", ")),
            "do.clr STOP!",
            sep = "\n"
        ))
    }
    
    how <- match.arg(how)
    
    # From Seurat
    clr_trans <- function(x) {
        log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x))))
    }

    if (how == 'across_cells') {
        if (verbose) {
            message("Applying CLR across cells (repeat CLR per feature).")
        }
        
        dat_clr_applied <- dat_to_apply_clr[, use_cols, with = FALSE]
        
        # https://stackoverflow.com/questions/16846380/apply-a-function-to-every-specified-column-in-a-data-table-and-update-by-referen
        for (j in use_cols) {
            set(dat_clr_applied, j = j, value = clr_trans(dat_clr_applied[[j]]))
        }
        
        
    } else {
        if (verbose) {
            message("Applying CLR across features (repeat CLR per cell).")
        }
        
        dat_clr_applied <- dat_to_apply_clr[, use_cols, with = FALSE]
        dat_clr_applied <- apply(dat_clr_applied, 1, clr_trans)
        
        # weird how i have to transpose it, but whatever
        dat_clr_applied <- data.table(t(dat_clr_applied))
        
    }
    
    # add the cell_id back
    dat_clr_applied[, cell_id := dat_to_apply_clr[[dat@cell_id_col]]]
    setnames(dat_clr_applied, "cell_id", dat@cell_id_col)
    
    
    dat <- add.new.data(
        spectre_obj = dat,
        dat = dat_clr_applied,
        dat_name = output_name,
        metadata = list("clr_trans" = how)
    )
    
    return(dat)
}
