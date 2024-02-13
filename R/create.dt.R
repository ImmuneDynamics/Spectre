#' Create a data.table
#'
#' Function to convert a Seurat or SingleCellExperiment or flowFrame object into a list
#' containing a data.table object storing the data, as well as a bunch of vectors grouping
#' the various features (columns) in the data.table, e.g. dim.reds vector will list
#' the kind of dimensionality reduction (by PCA or UMAP for instance) that have been applied
#' to the data and henceforth the columns in the data.table associated with them.
#'
#' @param dat NO DEFAULT. A Seurat or SingleCellExperiment or flowFrame object.
#' @param from DEFAULT = NULL. The class of the dat parameter. By default, this will be detected
#' automatically, but this can be overwritten using from. Can be from =
#' 'Seurat' or 'SingleCellExperiment' or 'flowFrame'.
#'
#' @return A list containing several elements.
#' The first element is a data.table which concatenates all components of dat into a single data.table.
#' The remaining elements contain the metadata of the dat such as the gene names,
#' cell barcodes, top x highly expressed genes, assays, etc.
#'
#'
#' @examples
#' devtools::install_github("satijalab/seurat-data")
#' library(SeuratData)
#' InstallData("pbmc3k")
#' pbmc <- LoadData("pbmc3k", type = "pbmc3k.final")
#' cell.dat <- create.dt(pbmc)
#'
#' library(scRNAseq)
#' sce <- ReprocessedAllenData("tophat_counts")
#  counts(sce) <- assay(sce, "tophat_counts")
#  cell.dat <- create.dt(sce)
#'
#' @author
#' Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#' Givanna Putri
#'
#'
#' @export 

create.dt <- function(dat, from = NULL) {
    # Umcomment for development.
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

    supported_type <-
        c("Seurat", "SeuratObject", "SingleCellExperiment", "flowFrame")
    error_msg <- paste(
        "Could not determine the object type.",
        "Currently only Seurat objects are supported.",
        "You can also try manually specifying the object type using the 'from' argument.",
        "(e.g. from = 'Seurat' or from = 'SingleCellExperiment')"
    )

    if (!is.null(from) && !(from %in% supported_type)) {
        # The object type is not supported then
        stop(error_msg)
    }

    # Will only get here if the given from is supported
    if (!is.null(from)) {
        object.type <- from
    } else {
        # Automatically detect
        message("Inferring dat type")
        if (class(dat)[1] == "Seurat") {
            object.type <- "Seurat"
        } else if (class(dat)[1] == "SingleCellExperiment") {
            object.type <- "SingleCellExperiment"
        } else if (class(dat)[1] == "flowFrame") {
            object.type <- "flowFrame"
        } else {
            stop(error_msg)
        }
    }


    message(object.type, " detected")
    if (object.type == "Seurat") {
        res <- convert.seurat.object(dat = dat)
    } else if (object.type == "SingleCellExperiment") {
        res <- convert.sce(dat = dat)
    } else if (object.type == "flowFrame") {
        res <- convert.flowframe(dat = dat)
    }

    message(paste0(
        "Converted a ",
        object.type,
        " object into a data.table stored in a list"
    ))

    return(res)
}


#' Convert Seurat object
#'
#' Convert Seurat object to data.table
#'
#' @param dat Seurat object to convert
#'
#' @return A list of data.table objects.
#'
#' @noRd
#' 
convert.seurat.object <- function(dat) {
    
    # require: Seurat
    check_packages_installed(c("Seurat"))

    # a <- Seurat::GetAssayData(object = dat)
    assays <- names(dat@assays)
    dim.reds <- names(dat@reductions)

    # get HVGs.
    var.features <- Seurat::VariableFeatures(dat)
    if (length(var.features) > 0) {
        var.features.top10 <- head(Seurat::VariableFeatures(dat), 10)
    }

    # get gene names and cell names
    # geneNames <- a@Dimnames[[1]]
    # cellNames <- a@Dimnames[[2]]

    ### Start data.table

    res <- data.table::data.table(cellNames = colnames(dat))

    ### Add metadata
    if (!is.null(dat@meta.data)) {
        col.meta <- dat@meta.data
        col.meta <- data.table::as.data.table(col.meta)
        meta.cols <- names(col.meta)

        res <- cbind(res, col.meta)
    }

    ### Add counts from different assays and layers

    # each assay will have its own feature names (genes or ADTs)
    # store it so we can return it in a list later
    all_assays_feature_names <- list()
    
    # we store the layers name for each assay as well
    all_assays_layers <- list()
    
    for (assay in assays) {
        
        # slot is seurat is deprecated now.
        # it is called layers.
        # This is Seurat's Layers...
        layers <- SeuratObject::Layers(dat[[assay]])
        all_assays_layers[[assay]] <- layers
        
        all_assays_feature_names[[assay]] <- list()
        
        for (layer in layers) {
            
            cnt_mat <- SeuratObject::LayerData(dat, assay=assay, layer=layer)
            
            # if there is scale.data layer, the feature names are unique to the layer
            # because not all genes will be in the layer
            feature_names <- rownames(cnt_mat)
            all_assays_feature_names[[assay]][[layer]] <- feature_names
            
            cnt_mat <- data.table::as.data.table(cnt_mat)
            cnt_mat <- data.table::transpose(cnt_mat)
            names(cnt_mat) <- paste(feature_names, assay, layer, sep = "_")
            
            res <- cbind(res, cnt_mat)
            
            rm(cnt_mat)
            rm(feature_names)
            
        }
    }

    ### Add dim reductions

    for (i in dim.reds) {
        # i <- dim.reds[[2]]

        tmp <- dat@reductions[[i]]@cell.embeddings
        tmp <- data.table::as.data.table(tmp)

        names(tmp) <- paste0(i, "_", names(tmp))
        res <- cbind(res, tmp)
    }

    ### Wrap up and return

    final.res <- list()

    final.res$data.table <- res
    final.res$featureNames <- all_assays_feature_names
    final.res$cellNames <- colnames(dat)

    final.res$meta.data <- meta.cols

    if (length(var.features) > 0) {
        final.res$var.features <- var.features
        final.res$var.features.top10 <- var.features.top10
    }
    
    # No need to check for NULL, it will still work anyway
    final.res$assays <- paste0("_", assays)
    final.res$layers <- all_assays_layers
    final.res$dim.reds <- paste0(dim.reds, "_")

    return(final.res)
}


#' Convert SCE object.
#'
#' Convert SCE object into data.table.
#'
#' @param dat SingleCellExperiment object to convert.
#'
#' @return A list of data.table objects.
#'
#' @noRd
#' 
convert.sce <- function(dat) {
    
    # require: SingleCellExperiment
    check_packages_installed(c("SingleCellExperiment"))
    
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

    res <- data.table::as.data.table(cellNames)
    names(res) <- "cellNames"

    ### Add metadata

    message("-- Adding metadata")

    if (!is.null(dat@colData)) {
        col.meta <- dat@colData
        col.meta <- data.table::as.data.table(col.meta)
        meta.cols <- names(col.meta)

        res <- cbind(res, col.meta)
        res
    }

    ### Assay data

    message("-- Adding assay data")

    for (i in assays) {
        # i <- assays[[1]]

        tmp <- dat@assays@data[[i]]
        tmp <- Matrix::as.matrix(tmp)
        # dim(tmp)
        # colnames(tmp) # columns = cells
        tmp <- data.table::as.data.table(tmp)
        tmp <- data.table::transpose(tmp)
        names(tmp) <- geneNames
        names(tmp) <- paste0(names(tmp), "_", i)

        res <- cbind(res, tmp)

        rm(i)
        rm(tmp)
    }

    ### DimRed data

    message("-- Adding DimRed data")

    for (i in dim.reds) {
        # i <- dim.reds[[2]]

        tmp <- reducedDims(dat)[[i]]
        # dim(tmp)
        # colnames(tmp) # columns = DRs
        tmp <- data.table::as.data.table(tmp) # Dont transpose

        new.names <- paste0(i, "_", c(1:length(names(tmp))))
        names(tmp) <- new.names

        res <- cbind(res, tmp)

        rm(i)
        rm(tmp)
    }

    ### Wrap up and return

    message("-- Finalising")

    final.res <- list()

    final.res$data.table <- res
    final.res$geneNames <- geneNames
    final.res$cellNames <- cellNames
    final.res$meta.data <- meta.cols
    final.res$assays <- paste0("_", assays)
    final.res$dim.reds <- paste0(dim.reds, "_")

    return(final.res)
}



#' Convert FlowFrame object.
#'
#' Convert FlowFrame object into data.table.
#' @param dat FlowFrame object to convert.
#'
#' @return A list of data.table objects.
#' 
#' @noRd
#' 
convert.flowframe <- function(dat) {
    
    # require: flowCore
    
    ### Extract 'exprs'

    res <- exprs(dat)
    res <- res[1:nrow(res), 1:ncol(res)]
    res <- data.table::as.data.table(res)

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
