#' Run rPCA for batch correction
#' 
#' Use Seurat rPCA to do batch correction.
#' For more details on rPCA, see https://satijalab.org/seurat/articles/integration_rpca.html.
#'
#' @param dat NO DEFAULT. 
#' A Spectre object containing the data to apply batch correction to.
#' @param data_source NO DEFAULT. Character. The name of the data in Spectre object to 
#' apply batch correction to.
#' @param output_name NO DEFAULT. Character. The name of the data in Spectre object to 
#' store the batch corrected data to.
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
#'
#' @return Spectre object with new element.
#'
#' @examples
#' dat_raw = Spectre::demo.batches
#' dat_raw[, cell_id := paste0("Cell_", seq(nrow(dat_raw)))]
#' dat = create.spectre.object(cell_id_col = "cell_id")
#' dat = add.new.data(spectre_obj = dat, dat = dat_raw, "cyto_batch")
#' 
#' markers = c("CD45_chn", "CD48_chn", "CD117_chn", 
#' "CD11b_chn", "SiglecF_chn", "NK11_chn", "B220_chn", 
#' "CD8a_chn", "CD4_chn", "Ly6C_chn", "Ly6G_chn", "CD115_chn", 
#'          "CD3e_chn", "CD16.32_chn", "MHCII_chn")
#' 
#' dat = run.rpca.batch.correction(
#'  dat = dat,
#'  data_source = "cyto_batch",
#'  output_name = "cyto_batch_corrected",
#'  use_cols = markers,
#'  batch_col = "Batch",
#'  verbose = FALSE,
#'  reference_batch = NULL
#' )
#' 
#' @export
#' 
run.rpca.batch.correction <- function(
        dat,
        data_source,
        output_name,
        use_cols,
        batch_col,
        reference_batch = NULL,
        k_anchor = 5,
        seed = 42,
        verbose = TRUE
) {
    
    
    # for testing only
    # dat_raw = Spectre::demo.batches
    # dat_raw[, cell_id := paste0("Cell_", seq(nrow(dat_raw)))]
    # dat = create.spectre.object(cell_id_col = "cell_id")
    # dat = add.new.data(spectre_obj = dat, dat = dat_raw, "cyto_batch")
    # 
    # data_source = "cyto_batch"
    # output_name = "cyto_batch_corrected"
    # use_cols = names(dat$cyto_batch)[1:15]
    # batch_col = "Batch"
    # reference_batch = NULL
    # k_anchor = 5
    # seed = 42
    # verbose = TRUE
    
    if (!is(dat, "Spectre")) {
        stop("dat must be of class Spectre")
    }
    
    check_packages_installed("Seurat")
    
    if (verbose) {
        message('Running rPCA batch correction.')
        message('(1/6) setting up data.')
    }
    
    # TODO need to first check the markers exists in use_cols
    
    # just incase thre are some non-standard naming and seurat obj complained
    new_col_name <- paste0("Col", seq(length(use_cols)))
    names(new_col_name) <- use_cols
    
    batches <- unique(dat[[data_source]][[batch_col]])
    
    cell_id_col <- dat@cell_id_col
    
    seurat_objs <- lapply(batches, function(batch) {
        
        # batch <- batches[1]
        
        cnt_mtx <- dat[[data_source]]
        cnt_mtx <- cnt_mtx[cnt_mtx[[batch_col]] == batch, c(use_cols, cell_id_col), with = FALSE]
        setnames(cnt_mtx, names(new_col_name), new_col_name)
        
        sparse_cnt_mtx <- t(as.matrix(cnt_mtx[, new_col_name, with = FALSE]))
        colnames(sparse_cnt_mtx) <- cnt_mtx[[cell_id_col]]
        sparse_cnt_mtx <- Matrix::Matrix(sparse_cnt_mtx)
        
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
    integ_features <- Seurat::SelectIntegrationFeatures(
        object.list = seurat_objs, 
        verbose=verbose
    )
    
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
            reduction = 'rpca',
            verbose = verbose)
        
        
    } else {
        # TODO have to make sure the dims is 1 less than number of features we have. Otherwise it gives stupid error.
        immune_anchors <- Seurat::FindIntegrationAnchors(
            object.list = seurat_objs, 
            anchor.features = integ_features, 
            dims = seq(length(integ_features)-1), 
            k.anchor = k_anchor,
            reduction = 'rpca',
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
    cell_id_batch_info <- dat[[data_source]][, c(cell_id_col, batch_col), with = FALSE]
    # because sort is set to false, the order of cell_id_batch_info is preserved.
    # Purrrfect!
    batch_corrected_dat <- merge.data.table(
        cell_id_batch_info,
        batch_corrected_dat,
        by = cell_id_col,
        sort = FALSE
    )
    
    dat <- add.new.data(dat, batch_corrected_dat, output_name)
    
    # check umap. very simple check just to see rpca is correcting the batches.
    # umap_dat <- run.umap(dat$cyto_batch_corrected, use.cols = use_cols)
    # make.colour.plot(umap_dat, "UMAP_X", "UMAP_Y", "Batch", randomise.order = FALSE)
    # 
    # umap_pre_cor <- run.umap(dat$cyto_batch, use.cols = use_cols)
    # make.colour.plot(umap_pre_cor, "UMAP_X", "UMAP_Y", "Batch", randomise.order = FALSE)
    
    return(dat)
    
    
}













