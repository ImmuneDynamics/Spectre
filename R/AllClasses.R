#' SpectreObject Class
#'
#' @slot cytometry_data A data.table slot storing cytometry data. 
#' Rows are cells, columns are markers/features.
#' @slot citeseq_data A data.table slot storing CITEseq data.
#' Rows are cells, columns are antibodies.
#' 
#' 
#' @examples
#' # Create an instance of SpectreObject with empty slots
#' 
#' @details
#' You should use this object if you want to jointly analyse cytometry and 
#' CITEseq data.
#' 
#' If you only want to analyse cytometry data, read the file in as data.table
#' and run Spectre functions as per normal.
#' There is no need to use this object.
#' 
#' Both `cytometry_data` and `citeseq_data` slots must be filled.
#' By default, if you don't pass any data.table into both slots, they will
#' both be empty. 
#' If you specify something for `cytometry_data` you MUST specify a data.table
#' for `citeseq_data`.
#' 
#' @export
#' 
#' @name SpectreObject-class
#' @docType class
#' @rdname SpectreObject-class
#'
setClass("SpectreObject",
         representation(
             cytometry_data = "data.table", 
             citeseq_data = "data.table"
             ),
         prototype(citeseq_data = data.table())
         
)


#' @export
#' @rdname SpectreObject-class
SpectreObject <- function(cytometry_data, citeseq_data = data.table()) {
    if (!is.data.table(cytometry_data)) {
        stop("cytometry_data slot must be a data.table object")
    }
    if (!is.data.table(citeseq_data)) {
        stop("citeseq_data slot must be a data.table object")
    }
    return(new("SpectreObject", 
               cytometry_data = cytometry_data, 
               citeseq_data = citeseq_data))
}


#' Print SpectreObject
#'
#' @param SpectreObject A SpectreObject 
#'
#' @export
#'
setMethod("show", "SpectreObject", function(object) {
    cat("class: SpectreObject\n")
    
    cat("cytometry_data slot:", class(object@cytometry_data), "with",
        ncol(object@cytometry_data), "features", 
        nrow(object@cytometry_data), "cells")
    
    cat("\n")
    
    cat("citeseq_data slot:", class(object@citeseq_data), "with",
        ncol(object@citeseq_data), "features", 
        nrow(object@citeseq_data), "cells")
    
    # if (is.null(object@cytometry_data)) {
    #     cat("cytometry_data slot: NULL with 0 features, 0 cells")
    # } else {
    #     cat("cytometry_data slot:", class(object@cytometry_data), "with",
    #         ncol(object@cytometry_data), "features", 
    #         nrow(object@cytometry_data), "cells")
    # }
    # 
    # cat("\n")
    # 
    # if (is.null(object@citeseq_data)) {
    #     cat("citeseq_data slot:", class(object@citeseq_data)," with 0 features, 0 cells")
    # } else {
    #     cat("citeseq_data slot:", class(object@citeseq_data), "with",
    #         ncol(object@citeseq_data), "features", 
    #         nrow(object@citeseq_data), "cells")
    # }
})

# Need to ask Tom what this is doing..
# Spectre <- setClass(Class = "Spectre",
#                     slots = c(
#                       data = "list",
#                       analysis = "list",
#                       other = "list",
#                       key ='vector'
#                     )
# 
# Spectre <- setClass(Class = "Spectre",representation('list',
#                                                      name='vector',
#                                                      spatial="list",
#                                                      other="list"
# ))
# setMethod("show",
#           "Spectre",
#           function(object){
#             
#             ### Initial
#             
#             require('data.table')
#             
#             if(length(dat@name < 1)){
#               message(paste0("Spectre object: '", dat@name, "'"))
#             } else {
#               message("Spectre object")
#             }
#             
#             ### Datasets
#             
#             for(i in names(dat)){
#               r <- nrow(dat[[i]])
#               if(is.null(r)){
#                 r <- 0
#               }
#               c <- ncol(dat[[i]])
#               if(is.null(c)){
#                 c <- 0
#               }
#               message(paste0('  $', i, ' (', r, ' rows, ', c, ' cols)'))
#             }
#             
#             ### Spatial
#             
#             if(length(dat@spatial)){
#               message('  --')
#               
#               l <- length(names(dat@spatial))
# 
#               if(l <= 3){
#                 for(i in names(dat@spatial)){
#                   message(paste0('  @spatial$', i))
#                 }
#               } else {
#                 for(i in c(1:3)){
#                   a <- names(dat@spatial)[i]
#                   message(paste0('  @spatial$', a))
#                 }
#                 message(paste0('   +', l-3))
#               }
#             }
#             
#             ### Other
#             
#             if(length(dat@other)){
#               message('  --')
#               
#               for(i in names(dat@other)){
#                 message(paste0('  @other$', i))
#               }
#             }
#             
#           })

