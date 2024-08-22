#' Run cyCombine
#' 
#' Run cyCombine algorithm to do batch correction.
#' 
#' @references Pedersen, C.B., Dam, S.H., Barnkob, M.B. et al. 
#' cyCombine allows for robust integration of single-cell cytometry datasets 
#' within and across technologies. 
#' Nat Commun 13, 1698 (2022). 
#' https://doi.org/10.1038/s41467-022-29383-5
#' 
#' @references Ashhurst, T.M., Marsh‐Wakefield, F., Putri, G.H. et al.
#' Integration, exploration, and analysis of high‐dimensional single‐cell 
#' cytometry data using Spectre. 
#' Cytometry Part A, 101(3), pp.237-253 (2022).
#' https://doi.org/10.1002/cyto.a.24350
#' 
#' @param dat NO DEFAULT. 
#' A Spectre object containing the data to apply batch correction to.
#' @param use_cols NO DEFAULT. 
#' A vector of character column names to apply batch correction to.
#' @param batch_col NO default. 
#' The column that denotes the batch or dataset that each cell belongs to
#' @param xdim DEFAULT 8. X-dimension size of the SOM grid.
#' @param ydim DEFAULT 8. Y-dimension size of the SOM grid.
#' @param rlen DEFAULT 10. Number of times the data is presented to the SOM network.
#' Higher values are recommended if 10 does not appear to perform well.
#' @param parametric DEFAULT TRUE.
#' If TRUE, the parametric version of ComBat is used. 
#' If FALSE, the non-parametric version is used.
#' @param seed DEFAULT 42. Seed to use when creating the SOM.
#' @param covar DEFAULT NULL. The covariate ComBat should use. 
#' Can be a vector or a column name in the Spectre object (element data_source). 
#' If NULL, no covar will be used
#' @param anchor DEFAULT NULL. 
#' Experimental: A column or vector specifying which samples are replicates and 
#' which are not. 
#' If specified, this column will be used as a covariate in ComBat. 
#' Be aware that it may be confounded with the condition.
#' @param norm_method DEFAULT scale. 
#' Normalization method. Should be either 'rank', 'scale' or 'qnorm'.
#' @param ties_method DEFAULT average.
#' The method to handle ties, when using rank. 
#' See ?rank for other options.
#' @param verbose DEFAULT TRUE 
#' Whether to print progress messages. TRUE to print, FALSE to suppress.
#' @param cell_id_col Character. The column in `dat` denoting
#' the unique identifier of the cells.
#'
#' @return batch corrected data.table
#' 
#' @export
#'
run.cycombine <- function(
        dat,
        use_cols,
        batch_col,
        cell_id_col,
        xdim = 8,
        ydim = 8,
        rlen = 10,
        parametric = TRUE,
        seed = 42,
        covar = NULL,
        anchor = NULL,
        norm_method = "scale",
        ties_method = "average",
        verbose = TRUE) {
    
    # TODO cycombine has label parameter which according to the code
    # can be used as a substitute to clustering the data using SOM.
    # Need to look into this if we want to offer it as a parameter?
    #
    # TODO this is giving out a weird warnings Unknown or uninitialised column: `sample`.
    # No idea where it comes from..
    
    check_packages_installed(c("cyCombine"))
    
    if (verbose) {
        message("Running cyCombine")
        message("(1/3) Preparing data.")
    }
    
    
    # have to copy because otherwise it will change the column name of the
    # spectre object to "batch".
    # TODO maybe I don't need to do this?
    df <- copy(dat)
    
    # cyCombine "batch" column
    setnames(df, batch_col, "batch")
    
    # TODO reassess if we need this
    # setnames(df, "sample_id", "sample")
    # df <- df[, c(use_cols, "batch"), with = FALSE]
    
    if (verbose) {
        message("(2/3) Running cyCombine")
    }
    
    cycombine_res <- cyCombine::batch_correct(
        df = df,
        label = NULL,
        xdim = xdim,
        ydim = ydim,
        rlen = rlen,
        parametric = parametric,
        seed = seed,
        covar = covar,
        anchor = anchor,
        markers = use_cols,
        norm_method = norm_method,
        ties.method = ties_method
    )
    
    if (verbose) {
        message("(3/3) Constructing final data")
    }
    
    cycombine_res <- data.table(cycombine_res)
    
    cols_in_df <- names(df)
    
    # a bit of wrangling so we return the same set of columns as input
    res_to_return <- cycombine_res[, cols_in_df, with = FALSE]
    setnames(res_to_return, "batch", batch_col)
    
    # save the metadata
    
    param_metadata <- data.table(
        label = NULL,
        xdim = xdim,
        ydim = ydim,
        rlen = rlen,
        parametric = parametric,
        seed = seed,
        covar = paste(covar, collapse = ", "),
        anchor = paste(covar, collapse = ", "),
        markers = paste(use_cols, collapse = ", "),
        norm_method = norm_method,
        ties.method = ties_method
    )
    param_metadata <- as.data.table(t(param_metadata), keep.rownames = TRUE)
    setnames(param_metadata, c("rn", "V1"), c("parameter", "value"))
    
    cycombine_extra_info <- c(cell_id_col, setdiff(names(cycombine_res), names(df)))
    
    metadata <- list(
        "parameter" = param_metadata,
        "cycombine_extra_info" = cycombine_res[, cycombine_extra_info, with = FALSE]
        
    )
    
    for (meta_name in names(metadata)) {
        # meta_name = names(metadata)[1]
        attr(res_to_return, meta_name) <- metadata[[meta_name]]
        
    }
    
    return(res_to_return)
    
}










