#' Create Spectre object
#' 
#' @param cell_id_col Character. The column in the list of data.tables denoting
#' the unique identifier of the cells.
#' @param name Character vector. DEFAULT = NULL.
#' Describes the experiment/study associated with the object.
#'
#' @export
create.spectre.object <- function(cell_id_col, name = NULL) {
    
    if (is.null(name)) {
        return(new("Spectre", cell_id_col=cell_id_col))
    }
    
    return(new("Spectre", cell_id_col=cell_id_col, name=name))
}

#' Add new data to Spectre object
#' 
#' Add a new data to the list in Spectre object.
#'
#' @param spectre_obj A Spectre object.
#' @param dat A data.table to add to the object.
#' @param dat_name The name of the data.
#' @param metadata A named List. DEFAULT = NULL.
#' Additional information (metadata) about the data, including its description.
#' Each element in the list must be associated with a name.
#' If this is NULL, a new element will be created in `spectre_obj`'s metadata
#' slot, with `dat_name` as the name.
#'
#' @return An updated Spectre object.
#' @export
#'
#' @examples
#' dat <- create.spectre.object("cell_id")
#' x <- Spectre::demo.asinh
#' x$cell_id <- paste0("cell_", seq(1, nrow(x)))
#' dat <- add.new.data(dat, x, "cyto_asinh", metadata = list("description"= "test"))
#' dat
add.new.data <- function(spectre_obj, dat, dat_name, metadata = NULL) {
    
    if (!is(spectre_obj, "Spectre")) {
        stop("spectre_obj must be of class Spectre")
    }
        
    if (!spectre_obj@cell_id_col %in% names(dat)) {
        stop(paste("dat must contain", spectre_obj@cell_id_col, "column which uniquely identify the cells!"))
    }
    
    spectre_obj[[dat_name]] <- dat
    
    if (is.null(metadata)) {
        spectre_obj@metadata[[dat_name]] <- list()
    } else {
        spectre_obj@metadata[[dat_name]] <- metadata
    }
    
    
    return(spectre_obj)
}

#' Add new metadata to Spectre object
#' 
#' Add new information for a given dataset to the "metadata" slot in Spectre object.
#' Please see details for more information about how edge cases are handled.
#'
#' @param spectre_obj A Spectre object.
#' @param metadata An object containing the metadata to add.
#' @param metadata_name A name associated with the metadata.
#' For example, if the metadata describe cofactors used to do asinh transformation,
#' an appropriate value for this parameter will be something like "asinh_confactor".
#' @param dataset_name The name of the dataset to associate the metadata with.
#'
#' @return An updated Spectre object.
#' 
#' @details
#' If a metadata_name for the dataset_name already exists in the metadata
#' slot, it will be overwritten.
#' 
#' If dataset_name does not exist, i.e. you have either not added the dataset
#' associated with dataset_name or you have manually added the dataset but forgot
#' to update the metadata slot (i.e. not using add.new.data function),
#' it will produce an error.
#' 
#' @export
#'
#' @examples
#' dat <- create.spectre.object("cell_id")
#' metadata_dt <- data.table(marker=paste0("marker_", seq(1,10)), cofactors = seq(1,10))
#' dat <- add.new.metadata(dat, metadata_dt, "asinh_cofactors")
#' dat
add.new.metadata <- function(spectre_obj, metadata, metadata_name, dataset_name) {
    if (!is(spectre_obj, "Spectre")) {
        stop("spectre_obj must be of class Spectre")
    }
        
    # check if the dataset associated with the metadata exist.
    if (! dataset_name %in% names(spectre_obj)) {
        stop(paste(
            dataset_name, "is not present in the Spectre object!\n",
            "Only the following datasets are present:\n",
            paste(names(spectre_obj), collapse = ", ")
        ))
    }
    
    # check if the metadata already exist?
    if (metadata_name %in% names(spectre_obj@metadata[[dataset_name]])) {
        warning(paste(
            "A metadata with name", metadata_name,
            "is already present for dataset", dataset_name,
            ".\n Overwriting", metadata_name, "with new metadata."
        ))
    }
    
    spectre_obj@metadata[[dataset_name]][[metadata_name]] <- metadata
    return(spectre_obj)
}






