#' Check packages are installed
#' 
#' To check required packages are installed
#' 
check_packages_installed <- function(packages_name) {
    for (p in packages_name) {
        if (!is.element(p, installed.packages()[, 1])) 
            stop(paste(p, "is required but not installed"))
    }
}

#' Add columns in one data table (dat_to_add) to another (dat).
#' Replace columns with the same name
add_to_data_table <- function(dat, dat_to_add) {
    cols <- names(dat_to_add)
    
    same_cols <- intersect(cols, names(dat))
    if (length(same_cols) > 0) {
        warning(paste(
            "These columns in input data will be replaced!",
            paste(same_cols, collapse = ", ")
        ))
    }
    
    dat[, (cols) := dat_to_add[, cols, with=FALSE]]
    
    return(dat)
}
