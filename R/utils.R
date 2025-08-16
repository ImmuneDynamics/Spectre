#' Internal function to check if columns in use.cols are numeric
#' 
#' @param dat data.frame. Input data.
#' @param use.cols character. Columns to check.
#' 
#' @return TRUE if all columns in use.cols are numeric.
#' 
#' @noRd 
#' @keywords internal
#' 
.check_numeric_columns <- function(dat, use.cols) {
    # Check markers in use.cols exist in dat
    if (!all(use.cols %in% colnames(dat))) {
        missing_cols <- setdiff(use.cols, colnames(dat))
        stop(sprintf(
            "The following columns are not found in the data: %s", 
            paste(missing_cols, collapse = ", ")
        ))
    }
    
    # Check if all columns in use.cols are numeric
    non_numeric_cols <- use.cols[
        !vapply(dat[, ..use.cols], is.numeric, logical(1))
    ]
    
    if (length(non_numeric_cols) > 0) {
        stop(sprintf(
            "Non-numeric columns are in use.cols: %s.", 
            paste(non_numeric_cols, collapse = ", ")
        ))
    }
    
    return(TRUE)
}

#' Add or replace a column in a data.table
#' 
#' @param dt data.table. Input data.
#' @param col_name character. Name of the column to add or replace.
#' @param values character. Values to add or replace in the column.
#' 
#' @return data.table with the new column added or replaced.
#' 
#' @noRd 
#' @keywords internal
#'
.add_or_replace_column <- function(dt, col_name, values) {
    if (col_name %in% names(dt)) {
        warning(sprintf(
            "Column '%s' already exists and will be replaced.", 
            col_name
        ))
    }
    return(dt[, (col_name) := values])
}
