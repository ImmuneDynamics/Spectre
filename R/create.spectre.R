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
    
    # check if the data already exist?
    if (dat_name %in% names(spectre_obj)) {
        warning(paste(
            "A dataset with name", dat_name,
            "is already present.\n",
            "Overwriting it with the new data."
        ))
    }
    
    spectre_obj[[dat_name]] <- dat
    
    if (!is.null(metadata)) {
        spectre_obj <- add.new.metadata(
            spectre_obj = spectre_obj,
            metadata = metadata,
            dataset_name = dat_name
        )
    }
    
    
    return(spectre_obj)
}

#' Add new metadata to Spectre object
#' 
#' Add new information for a given dataset to the "metadata" slot in Spectre object.
#' Please see details for more information about how edge cases are handled.
#'
#' @param spectre_obj A Spectre object.
#' @param metadata A list of metadata to add to the Spectre object.
#' Each element must be associated with a name.
#' For example, if one of the elements stores the cofactors used to do asinh transformation,
#' this list must be something like `list("cofactors" = 5)`.
#' @param dataset_name The name of the dataset to associate the metadata with.
#'
#' @return An updated Spectre object.
#' 
#' @details
#' If a metadata_name for the dataset_name already exists, it will be overwritten.
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
    
    for (meta_name in names(metadata)) {
        # meta_name = names(metadata)[1]
        if (meta_name %in% names(attributes(spectre_obj[[dataset_name]]))) {
            warning(paste(
                dataset_name,
                "already has attribute with name",
                meta_name,
                ". Overwriting."
            ))
        }
        
        attr(spectre_obj[[dataset_name]], meta_name) <- metadata[[meta_name]]
        
    }
    
    return(spectre_obj)
}






