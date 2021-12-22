#' create.dt - convert a Seurat or SingleCellExperiment to a data.table
#'
#' This function converts a Seurat or SingleCellExperiment object into a list
#' containing a data.table, with vectors of gene and dimensionality reduction
#'
#' @param dat NO DEFAULT. A Seurat or SingleCellExperiment or flowFrame object.
#' @param from DEFAULT = NULL. By default, the class of object will be detected
#'   automatically, but this can be overwritten using from. Can be from =
#'   'Seurat' or 'SingleCellExperiment' or 'flowFrame'.
#' 
#' @return A list containing several elements.
#' The first element is a data.table which concatenates all components of dat into a single data.table.
#' The remaining elements contain the metadata of the dat such as the gene names,
#' cell barcodes, top x highly expressed genes, assays, etc.
#' 
#'
#' @usage create.dt(dat, from)
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://github.com/ImmuneDynamics/Spectre}.
#'
#' @import data.table
#'
#' @export create.dt

create.dt <- function(dat, from = NULL)
{
    ######################################################################################################
    ### Setup
    ######################################################################################################
    
    ### Packages
    
    require('data.table')
    
    ### Demo data
    
    # dat <- pbmc
    # dat <- pbmc.sce
    # from <- NULL
    
    # dat <- Spectre::demo.asinh
    # from <- NULL
    #
    # metadata <- data.frame(name=dimnames(dat)[[2]], desc=paste('column',dimnames(dat)[[2]],'from dataset'))
    # dat.ff <- new("flowFrame",
    #               exprs=as.matrix(dat), # in order to create a flow frame, data needs to be read as matrix
    #               parameters=Biobase::AnnotatedDataFrame(metadata))
    #
    # head(flowCore::exprs(dat.ff))
    # dat <- dat.ff
    #
    # rm(dat.ff)
    # rm(metadata)
    
    ### Determine class of object
    
    supported_type <- c("Seurat", "SingleCellExperiment", "flowFrame")
    error_msg <- paste("Could not determine the object type.", 
                       "Currently only Seurat objects are supported.", 
                       "You can also try manually specifying the object type using the 'from' argument.",
                       "(e.g. from = 'Seurat' or from = 'SingleCellExperiment')")
    
    if (!is.null(from) && !(from %in% supported_type)) {
        # The object type is not supported then
        stop(error_msg)
    }
    
    # Will only get here if the given from is supported
    if (!is.null(from)) {
        object.type <- from
    }
    else {
        # Automatically detect
        message("Inferring dat type")
        if(class(dat)[1] == 'Seurat'){
            object.type <- "Seurat"
        }
        
        else if(class(dat)[1] == 'SingleCellExperiment'){
            object.type <- "SingleCellExperiment"
        }
        
        else if(class(dat)[1] == 'flowFrame'){
            object.type <- "flowFrame"
        }
        
        else {
            stop(error_msg)
        }
    }

    ######################################################################################################
    ### OPTION: Seurat objects
    ######################################################################################################
    
    if (object.type == "Seurat") {
        message(object.type, ' detected')
        
        ### Packages
        
        require('dplyr')
        require('Seurat')
        require('patchwork')
        
        ### Extract
        
        a <- GetAssayData(object = dat)
        assays <- names(dat@assays)
        dim.reds <- names(dat@reductions)
        
        var.features <- VariableFeatures(dat)
        
        if (length(var.features) > 0) {
            var.features.top10 <- head(VariableFeatures(dat), 10)
        }
        
        ### Genes and cells
        
        geneNames <- a@Dimnames[[1]]
        cellNames <- a@Dimnames[[2]]
        
        ### Start data.table
        
        res <- as.data.table(cellNames)
        names(res) <- 'cellNames'
        
        ### Add metadata
        
        if (!is.null(dat@meta.data)) {
            col.meta <- dat@meta.data
            col.meta <- as.data.table(col.meta)
            meta.cols <- names(col.meta)
            
            res <- cbind(res, col.meta)
        }
        
        ### Add data from slots
        
        for (i in assays) {
            # i <- assays[[1]]
            
            types <- vector()
            
            if (ncol(dat@assays[[i]]@counts) > 0) {
                types <- c(types, 'counts')
                
                x1 <-
                    GetAssayData(object = dat,
                                 assay = i,
                                 slot = 'counts')
                x1 <- as.data.table(x1)
                x1 <- data.table::transpose(x1)
                names(x1) <- geneNames
                names(x1) <-
                    paste0(names(x1), "_", i, "_", 'counts')
                res <- cbind(res, x1)
            }
            
            if (ncol(dat@assays[[i]]@data) > 0) {
                types <- c(types, 'data')
                
                x2 <-
                    GetAssayData(object = dat,
                                 assay = i,
                                 slot = 'data')
                x2 <- as.data.table(x2)
                x2 <- data.table::transpose(x2)
                names(x2) <- geneNames
                names(x2) <- paste0(names(x2), "_", i, "_", 'data')
                res <- cbind(res, x2)
            }
            
            if (ncol(dat@assays[[i]]@scale.data) > 0) {
                types <- c(types, 'scale.data')
                
                x3 <-
                    GetAssayData(object = dat,
                                 assay = i,
                                 slot = 'scale.data')
                x3 <- as.data.table(x3)
                x3 <- data.table::transpose(x3)
                names(x3) <- geneNames
                names(x3) <-
                    paste0(names(x3), "_", i, "_", 'scale.data')
                res <- cbind(res, x3)
            }
            
            rm(i)
            rm(x1)
            rm(x2)
            rm(x3)
        }
        
        ### Add dim reductions
        
        for (i in dim.reds) {
            # i <- dim.reds[[2]]
            
            tmp <- dat@reductions[[i]]@cell.embeddings
            tmp <- as.data.table(tmp)
            
            names(tmp) <- paste0(i, '_', names(tmp))
            res <- cbind(res, tmp)
        }
        
        ### Wrap up and return
        
        final.res <- list()
        
        final.res$data.table <- res
        final.res$geneNames <- geneNames
        final.res$cellNames <- cellNames
        
        final.res$meta.data <- meta.cols
        
        if (length(var.features) > 0) {
            final.res$var.features <- var.features
            final.res$var.features.top10 <- var.features.top10
        }
        
        final.res$assays <- paste0('_', assays)
        final.res$slots <- paste0('_', types)
        final.res$dim.reds <- paste0(dim.reds, '_')
        
        message(paste0(
            "Converted a ",
            object.type,
            " object into a data.table stored in a list"
        ))
        return(final.res)
    }
    
    ######################################################################################################
    ### OPTION: SingleCellExperiment
    ######################################################################################################
    
    if (object.type == "SingleCellExperiment") {
        message(object.type, ' detected')
        
        ### Packages
        
        require('SingleCellExperiment')
        
        ###
        
        geneNames <- rownames(dat) # genes
        geneNames
        
        cellNames <- colnames(dat) # cells
        cellNames
        
        ###
        assays <- names(dat@assays@data)
        assays
        
        dim.reds <- names(reducedDims(dat))
        dim.reds
        
        ###
        
        res <- as.data.table(cellNames)
        names(res) <- 'cellNames'
        
        ### Add metadata
        
        message('-- Adding metadata')
        
        if (!is.null(dat@colData)) {
            col.meta <- dat@colData
            col.meta <- as.data.table(col.meta)
            meta.cols <- names(col.meta)
            
            res <- cbind(res, col.meta)
            res
        }
        
        ### Assay data
        
        message('-- Adding assay data')
        
        for (i in assays) {
            # i <- assays[[1]]
            
            tmp <- dat@assays@data[[i]]
            tmp <- as.matrix(tmp)
            # dim(tmp)
            # colnames(tmp) # columns = cells
            tmp <- as.data.table(tmp)
            tmp <- data.table::transpose(tmp)
            names(tmp) <- geneNames
            names(tmp) <- paste0(names(tmp), "_", i)
            
            res <- cbind(res, tmp)
            
            rm(i)
            rm(tmp)
        }
        
        ### DimRed data
        
        message('-- Adding DimRed data')
        
        for (i in dim.reds) {
            # i <- dim.reds[[2]]
            
            tmp <- reducedDims(dat)[[i]]
            # dim(tmp)
            # colnames(tmp) # columns = DRs
            tmp <- as.data.table(tmp) # Dont transpose
            
            new.names <- paste0(i, "_", c(1:length(names(tmp))))
            names(tmp) <- new.names
            
            res <- cbind(res, tmp)
            
            rm(i)
            rm(tmp)
        }
        
        ### Wrap up and return
        
        message('-- Finalising')
        
        final.res <- list()
        
        final.res$data.table <- res
        final.res$geneNames <- geneNames
        final.res$cellNames <- cellNames
        final.res$meta.data <- meta.cols
        final.res$assays <- paste0('_', assays)
        final.res$dim.reds <- paste0(dim.reds, '_')
        
        message(paste0(
            "Converted a ",
            object.type,
            " object into a data.table stored in a list"
        ))
        return(final.res)
    }
    
    ######################################################################################################
    ### OPTION: flowFrames
    ######################################################################################################
    
    if (object.type == "flowFrame") {
        message(object.type, ' detected')
        
        ### Packages
        
        require('flowCore')
        require('data.table')
        
        ### Extract 'exprs'
        
        res <- exprs(dat)
        res <- res[1:nrow(res), 1:ncol(res)]
        res <- as.data.table(res)
        
        for (i in names(res)) {
            # i <- names(res)[1]
            
            if (!any(is.na(as.numeric(as.character(res[[i]]))))) {
                res[[i]] <- as.numeric(res[[i]])
            }
            
        }
        
        ### Setup list
        
        final.res <- list()
        
        final.res$data.table <- res
        final.res$parameters <- dat@parameters
        final.res$description <- dat@description
        
        ### Return
        
        message(paste0(
            "Converted a ",
            object.type,
            " object into a data.table stored in a list"
        ))
        return(final.res)
        
    }
    
}
