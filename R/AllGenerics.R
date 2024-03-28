#' Do ArcSinh transformation
#'
#' Transform data in selected columns using ArcSinh transformation with a specified co-factor.
#'
#' @seealso \url{https://sydneycytometry.org.au/spectre} for usage instructions and vignettes.
#' @references \url{https://sydneycytometry.org.au/spectre}
#'
#' @param dat NO DEFAULT. 
#' Either a data.table or a Spectre object to apply arc-sinh transformation to.
#' @param use.cols NO DEFAULT. 
#' A vector of character column names to apply arc-sinh transformation to.
#' @param cofactor DEFAULT = 5. Co-factor to use for arcsinh transformation.
#' Supply a data.table with column marker and cofactor if you want to apply different co-factor to different markers.
#' See examples.
#' @param append.cf DEFAULT = FALSE. Appends the co-factor used to the end of 
#' the name of the transformed columns.
#' @param reduce.noise DEFAULT = FALSE. This is an experimental calculation 
#' which should reduce noise from negative values. Use with caution.
#' @param digits DEFAULT = NULL. Number of decimal places as a limit, not used if NULL. 
#' Values beyond will be rounded. 
#' Equal to the number or less (i.e. if 9 is used, but only 5 digits are present, 
#' then 5 digits will be used). Important to control for small floating point 
#' error differences in different OS.
#' @param verbose DEFAULT = TRUE.
#' If TRUE, the function will print progress updates as it executes.
#' 
#' @examples
#' library(data.table)
#' # Assuming dat is a data.table
#' dat <- data.table(NK11=rnorm(10, 2), CD3=rnorm(10,1))
#' # Default co-factor 5 will be used for both NK11 and CD3
#' dat_asinh <- do.asinh(dat, use.cols=c("NK11", "CD3"))
#' dat_asinh
#' 
#' # Apply asinh to only NK11
#' dat_asinh <- do.asinh(dat, use.cols=c("NK11"))
#' dat_asinh
#' 
#' # Apply different co-factor to the markers
#' cofactor_dt <- data.table(marker=c("NK11", "CD3"), cofactor=c(5,10))
#' cofactor_dt
#' dat_asinh <- do.asinh(dat, use.cols=c("NK11", "CD3"), cofactor=cofactor_dt)
#' dat_asinh
#' 
#' 
#' @export
setGeneric("do.asinh", function(dat,
                                use.cols,
                                cofactor = 5,
                                append.cf = FALSE,
                                reduce.noise = FALSE,
                                digits = NULL,
                                verbose = TRUE,
                                ...) {
    standardGeneric("do.asinh")
    
})

#' Run rPCA
#' 
#' Use Seurat rPCA to do batch correction.
#' For more details on rPCA, see https://satijalab.org/seurat/articles/integration_rpca.html.
#'
#' @param dat NO DEFAULT. 
#' A Spectre object containing the data to apply batch correction to.
#' @param dat_name NO DEFAULT. Character. The name of the data in Spectre object to 
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
#' @exportMethod do.asinh
#' @rdname do.asinh
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
#' dat = run.rpca(
#'  dat = dat,
#'  dat_name = "cyto_batch",
#'  output_name = "cyto_batch_corrected",
#'  use_cols = markers,
#'  batch_col = "Batch",
#'  verbose = FALSE,
#'  reference_batch = NULL
#' )
setGeneric("run.rpca", function(dat,
                                dat_name,
                                output_name,
                                use_cols,
                                batch_col,
                                reference_batch = NULL,
                                k_anchor = 5,
                                seed = 42,
                                verbose = TRUE) {
    standardGeneric("run.rpca")
    
})


