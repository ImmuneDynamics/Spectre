#' Run CCA
#' 
#' Use Seurat CCA to do batch correction.
#' For more details on CCA, see https://satijalab.org/seurat/reference/runcca.
#'
#' @param dat NO DEFAULT. A data.table
#' @param use_cols NO DEFAULT. 
#' A vector of character column names to apply batch correction to.
#' @param batch_col Character. The column in the data in Spectre object that identifies
#' which batch each cell belongs to.
#' @param reference_batch DEFAULT NULL. Whether to align to batches to a given batch.
#' If yes, then supply this parameter with the name of the batch you want to align the other batches to.
#' @param k_anchor DEFAULT 5. Passed to Seurat's FindIntegrationAnchors function.
#' Essentially, it determines the number of neighbors (k) to use when 
#' Seurat's FindIntegrationAnchors is picking anchors.
#' @param seed DEFAULT 42. Seed used when running PCA. 
#' @param verbose DEFAULT = TRUE.
#' If TRUE, the function will print progress updates as it executes.
#' @param cell_id_col Character. The column in `dat` denoting
#' the unique identifier of the cells.
#'
#' @return batch corrected data.table
#' 
#' @usage run.cca(dat, use_cols, batch_col, cell_id_col,
#' reference_batch = NULL, k_anchor = 5, seed = 42, verbose = TRUE)
#'
#' @examples
#' library(data.table)
#' 
#' # download sample data
#' download.file(url = "https://github.com/ImmuneDynamics/data/blob/main/demo.batches.1.RData?raw=TRUE", destfile = 'demo.batches.1.RData', mode = 'w')
#' 
#' # this will store data in a variable call demo.batches.1
#' load("demo.batches.1.RData")
#' 
#' # clean up
#' unlink(c('demo.batches.1.RData'))
#' 
#' # cca can be slow. subsample so it doesn't take too long
#' demo.batches.1 <- Spectre::do.subsample(
#'     demo.batches.1, targets = rep(100, 2), 
#'     divide.by = "Batch"
#' )
#' 
#' # assign unique id to each cell
#' demo.batches.1[, cell_id := paste0("Cell_", seq(nrow(demo.batches.1)))]
#' 
#' markers <- c("CD45_chn", "CD48_chn", "CD117_chn",
#'              "CD11b_chn", "SiglecF_chn", "NK11_chn", "B220_chn",
#'              "CD8a_chn", "CD4_chn", "Ly6C_chn", "Ly6G_chn", "CD115_chn",
#'              "CD3e_chn", "CD16.32_chn", "MHCII_chn")
#' 
#' dat_batch_corrected = run.cca(
#'     dat = demo.batches.1,
#'     cell_id_col = "cell_id",
#'     use_cols = markers,
#'     batch_col = "Batch",
#'     verbose = FALSE,
#'     reference_batch = NULL
#' )
#' 
#' dat_batch_corrected
#' 
run.cca <- function(dat,
                    use_cols,
                    batch_col,
                    cell_id_col,
                    reference_batch = NULL,
                    k_anchor = 5,
                    seed = 42,
                    verbose = TRUE) {
    # for testing only
    # dat = Spectre::demo.batches
    # 
    # output_name = "cyto_batch_corrected"
    # use_cols = names(dat$cyto_batch)[1:15]
    # batch_col = "Batch"
    # reference_batch = NULL
    # k_anchor = 5
    # seed = 42
    # verbose = TRUE
    
    check_packages_installed("Seurat")
    
    if (verbose) {
        message('Running CCA.')
        message('(1/6) setting up data.')
    }
    
    # just incase thre are some non-standard naming and seurat obj complained
    new_col_name <- paste0("Col", seq(length(use_cols)))
    names(new_col_name) <- use_cols
    
    batches <- unique(dat[[batch_col]])
    
    seurat_objs <- lapply(batches, function(batch) {
        
        # batch <- batches[1]
        
        cnt_mtx <- dat[dat[[batch_col]] == batch, c(use_cols, cell_id_col), with = FALSE]
        setnames(cnt_mtx, names(new_col_name), new_col_name)
        
        sparse_cnt_mtx <- t(as.matrix(cnt_mtx[, new_col_name, with = FALSE]))
        colnames(sparse_cnt_mtx) <- cnt_mtx[[cell_id_col]]
        sparse_cnt_mtx <- Matrix::Matrix(sparse_cnt_mtx, sparse = TRUE)
        
        if (verbose) {
            message(paste('(2/6) creating Seurat object for batch', batch))
        }
        
        seurat_obj <- Seurat::CreateSeuratObject(counts = sparse_cnt_mtx, assay = 'cyto')
        # kind of useless as it will always return all the features? but needed for later
        seurat_obj <- Seurat::FindVariableFeatures(
            seurat_obj, selection.method = "vst", nfeatures = length(use_cols), 
            verbose=verbose
        )
        
        # have to do this as otherwise scale data will complain later
        seurat_obj[['cyto']]$data <- seurat_obj[['cyto']]$counts
        
        return(seurat_obj)
        
    })
    names(seurat_objs) <- batches
    
    ### Select integration features, scale data, and run PCA
    
    if (verbose) {
        message('(3/6) Performing scaling and PCA.')
    }
    
    # this will rank the features
    integ_features <- Seurat::SelectIntegrationFeatures(object.list = seurat_objs, verbose=verbose)
    
    seurat_objs <- lapply(seurat_objs, function(seurat_obj) {
        seurat_obj <- Seurat::ScaleData(seurat_obj, features = integ_features, 
                                        verbose = verbose, assay = 'cyto')
        # no need to set PCs as it will just default either 50 or less if we have less markers than 50
        # TODO not sure about the approx parameter to run standard svd instead. Note, the number of PCs will be the number of markers
        # if setting approx=FALSE. If approx is true, npcs will be number of markers - 1. Not sure why.
        # TODO manually setting the npcs rather than getting the function to infer it. Is this the best way?
        seurat_obj <- Seurat::RunPCA(seurat_obj, features = integ_features, verbose = verbose, 
                                     assay = 'cyto', seed.use = seed, npcs = length(use_cols))
        
        return(seurat_obj)
    })
    
    if (verbose) {
        message('(4/6) Finding integration anchors.')
    }
    
    if(is.null(reference_batch)){
        # TODO have to make sure the dims is 1 less than number of features we have. Otherwise it gives stupid error.
        immune_anchors <- Seurat::FindIntegrationAnchors(
            object.list = seurat_objs, 
            anchor.features = integ_features, 
            dims = seq(length(integ_features)-1), 
            k.anchor = k_anchor,
            reduction = 'cca',
            verbose = verbose)
        
        
    } else {
        # TODO have to make sure the dims is 1 less than number of features we have. Otherwise it gives stupid error.
        immune_anchors <- Seurat::FindIntegrationAnchors(
            object.list = seurat_objs, 
            anchor.features = integ_features, 
            dims = seq(length(integ_features)-1), 
            k.anchor = k_anchor,
            reduction = 'cca',
            reference = which(names(seurat_objs) == reference_batch),
            verbose = verbose)
    }
    
    ### Integrate the data
    
    if (verbose) {
        message('(5/6) Integrating data')
    }
    
    # no way to use standard svd. just have to put up with it for now.
    batch_corrected_seurat_obj <- Seurat::IntegrateData(
        anchorset = immune_anchors, 
        dims = seq(length(integ_features)-1),
        verbose = verbose
    )
    # just to be sure!
    Seurat::DefaultAssay(batch_corrected_seurat_obj) <- "integrated"
    
    ### Re-construct Spectre object
    
    if (verbose) {
        message('(6/6) Constructing final data')
    }
    
    batch_corrected_dat <- data.table::transpose(as.data.table(batch_corrected_seurat_obj[['integrated']]$data))
    names(batch_corrected_dat) <- rownames(batch_corrected_seurat_obj[['integrated']]$data)
    
    batch_corrected_dat[[cell_id_col]] <- Seurat::Cells(batch_corrected_seurat_obj)
    
    # rename the markers
    setnames(batch_corrected_dat, new_col_name, names(new_col_name))
    
    # order the cell id 
    cell_id_batch_info <- dat[, c(cell_id_col, batch_col), with = FALSE]
    # because sort is set to false, the order of cell_id_batch_info is preserved.
    # Purrrfect!
    batch_corrected_dat <- merge.data.table(
        cell_id_batch_info,
        batch_corrected_dat,
        by = cell_id_col,
        sort = FALSE
    )
    
    # check umap. very simple check just to see cca is correcting the batches.
    # umap_dat <- run.umap(dat$cyto_batch_corrected, use.cols = use_cols)
    # make.colour.plot(umap_dat, "UMAP_X", "UMAP_Y", "Batch", randomise.order = FALSE)
    # 
    # umap_pre_cor <- run.umap(dat$cyto_batch, use.cols = use_cols)
    # make.colour.plot(umap_pre_cor, "UMAP_X", "UMAP_Y", "Batch", randomise.order = FALSE)
    
    return(batch_corrected_dat)
    
    
}















