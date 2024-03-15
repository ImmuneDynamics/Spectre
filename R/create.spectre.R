#' Create Spectre object
#' 
#' @param cell_id_col Character. The column in the list of data.tables denoting
#' the unique identifier of the cells.
#'
#' @export
create.spectre.object <- function(cell_id_col) {
    return(new("Spectre", cell_id_col=cell_id_col))
}

#' Add new data to Spectre object
#' 
#' Add a new data to the list in Spectre object.
#'
#' @param spectre_obj A Spectre object.
#' @param dat A data.table to add to the object.
#' @param dat_name The name of the data.
#'
#' @return An updated Spectre object.
#' @export
#'
#' @examples
#' dat <- create.spectre.object("cell_id")
#' x <- Spectre::demo.asinh
#' x$cell_id <- paste0("cell_", seq(1, nrow(x)))
#' dat <- add.new.data(dat, x, "cyto_asinh")
#' dat
add.new.data <- function(spectre_obj, dat, dat_name) {
    
    if (!is(spectre_obj, "Spectre"))
        stop("spectre_obj must be of class Spectre")
    
    if (!spectre_obj@cell_id_col %in% names(dat)) {
        stop(paste("dat must contain", spectre_obj@cell_id_col, "column which uniquely identify the cells!"))
    }
    
    spectre_obj[[dat_name]] <- dat
    return(spectre_obj)
}