#' Run rPCA for label transfer
#' 
#' Use Seurat rPCA to perform label transfer from a CITEseq data to a cytometry data.
#'
#' @param dat NO DEFAULT. 
#' A Spectre object containing the cytometry and CITESeq data to 
#' perform label transfer.
#' @param  cytometry_data_source NO DEFAULT. Character. The name of the cytometry 
#' data in dat.
#' @param citeseq_data_source NO DEFAULT. Character. The name of the CITEseq 
#' data in dat.
#' @param output_name NO DEFAULT. Character. What name do you want to store
#' the resulting predicted label under in dat?
#' @param celltype_label_col NO DEFAULT. The column in the CITEseq data which
#' denotes the cell type/label each cell in the CITEseq data.
#' @param use_cols DEFAULT = NULL. 
#' A vector of marker/adt names to use to perform label transfer.
#' If left as the default value NULL, the function will automatically infer
#' markers/adts that are common between the CITEseq and cytometry data using
#' base R `intersect` function.
#' We highly recommend you manually specify this because the same proteins can be
#' spelled in my different ways (e.g., PD1 or PD-1). 
#' Additionally, same protein can have different names (e.g., CCR7 is equivalent
#' to CD197).
#' See details for more information.
#' @param k_anchor DEFAULT 20. Passed to Seurat's FindIntegrationAnchors function.
#' Essentially, it determines the number of neighbors (k) to use when 
#' Seurat's FindIntegrationAnchors is picking anchors.
#' @param seed DEFAULT 42. Seed used when running PCA. 
#' @param verbose DEFAULT = TRUE.
#' If TRUE, the function will print progress updates as it executes.
#' 
#' @details
#' ## Missing markers in use_cols
#' If you supply markers name in use_cols, and it turns out that some of them are
#' not available in either the CITEseq data, the cytometry data, or both, these
#' markers will be abandoned.
#' 
#' ## No common markers/adts were found
#' Regardless of whether you supply some markers in use_cols or you get the function
#' to infer them, if as it turns out there are no markers that are common between the 
#' CITEseq data and the cytometry data (and use_cols if you supply it with values),
#' the function will stop running.
#' 
#'
#' @return Spectre object with new element stored in output_name
#'
#' @export
#' 
run.rpca.label.transfer <- function(
        dat,
        cytometry_data_source,
        citeseq_data_source,
        output_name,
        celltype_label_col,
        use_cols = NULL,
        k_anchor = 20,
        seed = 42,
        verbose = TRUE
) {
    
    
    # for testing only
    # dat = qs::qread("~/Documents/spectre/v2/label_transfer/20240410_spectre_obj.qs")
    # cytometry_data_source = "cytof_bm_asinh"
    # citeseq_data_source = "citeseq_bm_clr"
    # output_name = "cytof_bm_asinh_rpca"
    # celltype_label_col = "FlowSOM_metacluster"
    # seed = 42
    # verbose = TRUE
    # k_anchor = 20
    # use_cols = NULL
    
    if (!is(dat, "Spectre")) {
        stop("dat must be of class Spectre")
    }
    
    check_packages_installed(c("Seurat"))
    
    if (verbose) {
        message('Running rPCA label transfer.')
    }
    
    cell_id_col <- dat@cell_id_col
    
    # Find common markers if use_cols is NULL
    # and if it isn't, do some checks on the use_cols
    if (is.null(use_cols)) {
        if (verbose) {
            message(paste(
                '(1/5) Inferring common markers between',
                citeseq_data_source,
                'and',
                cytometry_data_source
            ))
        }
        citeseq_adts <- names(dat[[citeseq_data_source]])
        cytometry_markers <- names(dat[[cytometry_data_source]])
        
        common_markers <- intersect(citeseq_adts, cytometry_markers)
        # remove cell_id column
        common_markers <- common_markers [common_markers != cell_id_col]
        
        if (length(common_markers) == 0) {
            stop(paste(
                "No common markers were found between",
                citeseq_data_source,
                'and',
                cytometry_data_source,
                "!!\nCannot proceed with label transfer!"
            ))
        } else {
            if (verbose) {
                message(paste(
                    "The following common markers will be used:",
                    paste(common_markers, collapse = ", "),
                    sep = "\n"
                ))
            }
        }
    } else {
        if (verbose) {
            message(paste(
                '(1/5) Checking markers in use_cols exist in both',
                citeseq_data_source,
                'and',
                cytometry_data_source
            ))
        }
        citeseq_adts <- intersect(use_cols, names(dat[[citeseq_data_source]]))
        cytometry_markers <- intersect(use_cols, names(dat[[cytometry_data_source]]))
        
        # Some checks and warnings to see whether all markers in use_cols 
        # are present in the citeseq and cytometry data
        if (length(citeseq_adts) != length(use_cols)) {
            warning(paste(
                "The following markers are missing in the CITEseq data:",
                paste(setdiff(use_cols,names(dat[[citeseq_data_source]])), collapse = ","),
                sep = "\n"
            ))
        }
        
        if (length(cytometry_markers) != length(use_cols)) {
            warning(paste(
                "The following markers are missing in the cytometry data:",
                paste(setdiff(use_cols,names(dat[[cytometry_data_source]])), collapse = ","),
                sep = "\n"
            ))
        }
        
        common_markers <- intersect(citeseq_adts, cytometry_markers)
        
        if (length(common_markers) == 0) {
            stop(paste(
                "No common markers were found between",
                citeseq_data_source,
                'and',
                cytometry_data_source,
                "!!\nCannot proceed with label transfer!"
            ))
        } else {
            if (verbose) {
                message(paste(
                    "The following common markers will be used:",
                    paste(common_markers, collapse = ", "),
                    sep = "\n"
                ))
            }
        }
    }
    
    # For renaming the markers because Seurat is a pain in the arse if there
    # are non alphanumeric character in the marker names
    
    placeholder_common_markers <- paste0("Col", seq(length(common_markers)))
    names(placeholder_common_markers) <- common_markers
    
    # subset citeseq to common markers only
    if (verbose) {
        message('(2/5) setting up data.')
    }
    
    citeseq_dat <- dat[[citeseq_data_source]][, common_markers, with=FALSE]
    setnames(citeseq_dat, common_markers, placeholder_common_markers)
    citeseq_dat <- Matrix::Matrix(t(as.matrix(citeseq_dat)), sparse=TRUE)
    colnames(citeseq_dat) <- dat[[citeseq_data_source]][[cell_id_col]]
    
    citeseq_dat <- Seurat::CreateSeuratObject(counts = citeseq_dat, data = citeseq_dat, assay = "ADT")
    
    # TODO I'm not sure if we should do this, but I can't think of any harm in doing it.
    # but again it is then putting the data in the scale.data slot which is not needed
    # for finding integration anchors.
    # citeseq_dat <- Seurat::ScaleData(citeseq_dat, features = placeholder_common_markers)
    
    citeseq_dat <- Seurat::AddMetaData(
        object = citeseq_dat,
        metadata = rep("citeseq", ncol(citeseq_dat)),
        col.name = "technology"
    )
    
    # subset citeseq to common markers only
    cyto_dat <- dat[[cytometry_data_source]][, common_markers, with=FALSE]
    setnames(cyto_dat, common_markers, placeholder_common_markers)
    cyto_dat <- Matrix::Matrix(t(as.matrix(cyto_dat)), sparse=TRUE)
    colnames(cyto_dat) <- dat[[cytometry_data_source]][[cell_id_col]]
    
    cyto_dat <- Seurat::CreateSeuratObject(counts = cyto_dat, data = cyto_dat, assay = "cyto")
    cyto_dat <- Seurat::AddMetaData(
        object = cyto_dat,
        metadata = rep("cytometry", ncol(cyto_dat)),
        col.name = "technology"
    )
    
    # Run rPCA
    if (verbose) {
        message('(3/5) Finding integration anchors.')
    }
    
    anchors <- Seurat::FindTransferAnchors(
        reference = citeseq_dat,
        query = cyto_dat,
        reference.assay = "ADT",
        query.assay = "cyto",
        features = placeholder_common_markers,
        normalization.method = 'LogNormalize',
        reduction = 'rpca',
        dims = 1:length(placeholder_common_markers),
        approx.pca = FALSE,
        k.anchor = k_anchor,
        npcs = length(placeholder_common_markers),
    )
    
    # casting the cell type labels to character as seurat doesn't like numeric..
    cell_type_labs <- dat[[citeseq_data_source]][[celltype_label_col]]
    
    if (is.numeric(cell_type_labs)) {
        # TODO probably need to improve this so users do not misunderstand it as
        # their data is modified.
        warning(paste(
            celltype_label_col,
            "in",
            citeseq_data_source,
            "is of type numeric. Converting this to character as rPCA does not accept numeric as cell type label."
        ))
        cell_type_labs <- as.character(cell_type_labs)
    }
    
    if (verbose) {
        message('(4/5) Transferring labels.')
    }
    
    cytometry_labels <- Seurat::TransferData(
        anchorset = anchors,
        refdata = cell_type_labs,
        weight.reduction='rpca.ref'
    )
    
    if (verbose) {
        message('(5/5) Constructing data and adding it to Spectre object.')
    }
    
    cytometry_labels <- as.data.table(cytometry_labels, keep.rownames = TRUE)
    setnames(cytometry_labels, c("rn", "predicted.id"), c("cell_id", "predicted_cell_type"))
    
    
    ### Re-construct Spectre object
    output_to_add <- cytometry_labels[, c("cell_id", "predicted_cell_type")]
    
    cytometry_labels$predicted_cell_type <- NULL
    
    dat = add.new.data(
        spectre_obj = dat, 
        dat = output_to_add, 
        dat_name = output_name,
        metadata = list("rpca_prediction_score" = cytometry_labels)
    )
    
    if (verbose) {
        message(paste(
            'rPCA label transfer finished. A new data and metadata with name',
            output_name,
            "has been added to the Spectre object.")
        )
    }
    
    # check umap. very simple check just to see rpca is correcting the batches.
    # umap_dat = run.umap(dat$cyto_batch_corrected, use.cols = use_cols)
    # make.colour.plot(umap_dat, "UMAP_X", "UMAP_Y", "Batch", randomise.order = FALSE)
    # 
    # umap_pre_cor = run.umap(dat$cyto_batch, use.cols = use_cols)
    # make.colour.plot(umap_pre_cor, "UMAP_X", "UMAP_Y", "Batch", randomise.order = FALSE)
    
    return(dat)
}













