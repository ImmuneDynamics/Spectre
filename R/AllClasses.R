#' Spectre Class
#' 
#' Spectre class is an extension the R native list class designed to store and 
#' manage cytometry data and the results of analysis associated with the data.
#' It is specifically developed to store them as a list of data.tables.
#' The types of data and analysis results that can be stored in this list includes
#' (but not limited to) raw data, asinh-transformed data, and results from downstream
#' analyses such as clustering or UMAP.
#' See details for more information.
#' 
#' @details
#' Each element within the list must be a data.table.
#' All of these elements *must* include a column which uniquely identify each cell.
#' The name of this column must be the same for all data.tables in the list and
#' must be speified in the cell_id_col slot. 
#' Elements in the list should be named descriptively (e.g., "cyto" for raw data, "cyto_asinh" for asinh-transformed data)
#' but custom naming is supported to accommodate different preferences (e.g., "cyto_raw" for raw data).
#'
#' @slot cell_id_col Character. Specifies the name of the column that uniquely identifies each cell within the data tables. 
#' This column must be present in all the data.tables in the list.
#' @slot name Character vector. Describes the dataset stored in the object.
#' This slot can be used to store the metadata for the data stored within the object.
#' @slot spatial Reserved for future use to store spatial data.
#' @slot other List. For storing any additional data that does not fit within the list 
#' or other defined slots.
#' 
#' 
#' @export
#' 
#' @name Spectre-class
#' @docType class
#' @rdname Spectre-class
#'
setClass("Spectre", 
         representation('list',
                        cell_id_col='character',
                        name='vector',
                        spatial="list",
                        other="list")
)

#' Print SpectreObject
#'
#' @param SpectreObject A SpectreObject 
#'
#' @export
#'
setMethod("show", "Spectre", function(object) {
    if(length(object@name < 1)){
        message(paste0("Spectre object: '", object@name, "'"))
    } else {
        message("Spectre object")
    }
    
    message("  ---")
    
    
    ### Datasets details
    
    for(i in names(object)){
        message(paste0('  $', i, ' (', nrow(object[[i]]), ' rows, ', ncol(object[[i]]), ' cols)'))
    }
    
    message("  ---")
    message(paste("  @cell_id_col:", object@cell_id_col))
    
    ### Spatial slot
    # TODO uncomment me when we have sorted spatial data.
    # if(length(object@spatial)){
    #     message('  ---')
    #     
    #     l <- length(names(object@spatial))
    #     
    #     if(l <= 3){
    #         for(i in names(object@spatial)){
    #             message(paste0('  @spatial$', i))
    #         }
    #     } else {
    #         for(i in c(1:3)){
    #             a <- names(object@spatial)[i]
    #             message(paste0('  @spatial$', a))
    #         }
    #         message(paste0('   +', l-3))
    #     }
    # }
    
    ### Other slot
    
    if(length(object@other)){
        message('  ---')
        
        for(i in names(object@other)){
            message(paste0('  @other$', i))
        }
    }
    
})



