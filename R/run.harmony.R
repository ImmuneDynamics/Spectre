#' Run Harmony batch alignment algorithm 
#' 
#' Use Harmony to do batch correction.
#' For more details on Harmony see: 
#' https://portals.broadinstitute.org/harmony/articles/quickstart.html.
#' 
#' @references Korsunsky, I., Millard, N., Fan, J. et al. 
#' Fast, sensitive and accurate integration of single-cell data with Harmony. 
#' Nat Methods 16, 1289–1296 (2019). 
#' https://doi.org/10.1038/s41592-019-0619-0
#' 
#' @references Ashhurst, T.M., Marsh‐Wakefield, F., Putri, G.H. et al.
#' Integration, exploration, and analysis of high‐dimensional single‐cell 
#' cytometry data using Spectre. 
#' Cytometry Part A, 101(3), pp.237-253 (2022).
#' https://doi.org/10.1002/cyto.a.24350
#' 
#' @param dat NO DEFAULT. A data.table.
#' @param cell_id_col Character. The column in `dat` denoting
#' the unique identifier of the cells.
#' @param use_cols NO DEFAULT. 
#' A vector of character column names to apply batch correction to.
#' @param batch_col NO default. The column that denotes the batch or dataset that each cell belongs to
#' @param seed DEFAULT 42. Seed used when harmony runs PCA. 
#' @param do_pca DEFAULT TRUE. Whether to perform PCA on the data.
#' @param npcs DEFAULT 20. If doing PCA on the data, the number of PCs to compute.
#' @param theta DEFAULT NULL. Diversity clustering penalty parameter. Specify for each
#' variable in vars_use. theta=0 does not encourage any
#' diversity. Larger values of theta result in more diverse clusters.
#' @param lambda DEFAULT NULL. Ridge regression penalty parameter. 
#' Specify for each batch.
#' Default lambda=1. Lambda must be strictly positive. Smaller values result
#' in more aggressive correction.
#' @param sigma DEFAULT 0.1.
#' Width of soft kmeans clusters. Default sigma=0.1. Sigma scales
#' the distance from a cell to cluster centroids. Larger values of sigma
#' result in cells assigned to more clusters. Smaller values of sigma make
#' soft kmeans cluster approach hard clustering.
#' @param nclust DEFAULT NULL. Number of clusters in model. 
#' nclust=1 equivalent to simple linear regression.
#' @param tau DEFAULT 0. Protection against overclustering small datasets with 
#' large ones. 'tau' is the expected number of cells per cluster.
#' @param plot_convergence DEFAULT FALSE.
#' Whether to print the convergence plot of the clustering objective function. 
#' TRUE to plot, FALSE to suppress. 
#' This can be useful for debugging.
#' @param return_object DEFAULT FALSE.
#' (Advanced Usage) Whether to return the Harmony object or only the corrected PCA embeddings.
#' If TRUE, this will be stored in the metadata slot.
#' @param verbose DEFAULT FALSE. 
#' Whether to print progress messages. TRUE to print, FALSE to suppress.
#' @param max_iter DEFAULT 10.
#' Maximum number of rounds to run Harmony. 
#' One round of Harmony involves one clustering and one correction step.
#' @param early_stop DEFAULT TRUE.
#' Enable early stopping for harmony. 
#' The harmonization process will stop when the change of 
#' objective function between corrections drops below 1e-4.
#' Only used if 'epsilon.harmony' is set to NULL, 
#' @param ncores DEFAULT 1.
#' Number of processors to be used for math operations when optimized BLAS is available. 
#' If BLAS is not supporting multithreaded then this option has no effect. 
#' By default, ncores=1 which runs as a single-threaded process. 
#' Although Harmony supports multiple cores, 
#' it is not optimized for multithreading. 
#' Increase this number for large datasets iff single-core performance is not 
#' adequate and optimized BLAS is available.
#' @param alpha DEFAULT 0.2. Advanced Harmony parameter.
#' When setting lambda = NULL and use lambda estimation mode, 
#' lambda would be determined by the expected number of cells assuming 
#' independence between batches and clusters. 
#' i.e., lambda = alpha * expected number of cells.
#' Alpha should be 0 < alpha < 1
#' @param block.size DEFAULT 0.05. Advanced Harmony parameter.
#' What proportion of cells to update during clustering. Between 0 to 1. 
#' Larger values may be faster but less accurate.
#' @param max.iter.cluster DEFAULT 20. Advanced Harmony parameter.
#' Maximum number of rounds to run clustering at each round of Harmony.
#' @param epsilon.cluster DEFAULT 0.001. Advanced Harmony parameter.
#' Convergence tolerance for clustering round of Harmony. 
#' Set to -Inf to never stop early.
#' @param epsilon.harmony DEFAULT 0.01. Advanced Harmony parameter.
#' Convergence tolerance for Harmony. 
#' Set to -Inf to never stop early. 
#' When 'epsilon.harmony' is set to not NULL, 
#' then user-supplied values of 'early_stop' is ignored.
#' 
#' @return Returns a data.table with aligned data added in new columns.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#' @author Givanna Putri
#'
#' @import data.table
#'
#' @export
#' 
run.harmony <- function(
        dat,
        cell_id_col,
        use_cols,
        batch_col,
        do_pca = TRUE,
        npcs = 20,
        theta = NULL,
        lambda = 1,
        sigma = 0.1,
        nclust = NULL,
        max_iter = 10,
        early_stop = TRUE,
        ncores = 1,
        plot_convergence = FALSE,
        return_object = FALSE,
        verbose = FALSE,
        seed = 42,
        # advanced parameters for Harmony
        alpha = 0.2,
        tau = 0,
        block.size = 0.05,
        max.iter.cluster = 200,
        epsilon.cluster = 0.001,
        epsilon.harmony = 0.01
) {
    # for testing only
    # download demo batches from immunedynamics data repo then load
    # dat = demo.batches.1
    # dat[, cell_id := paste0("Cell_", seq(nrow(dat)))]
    # use_cols = names(dat)[1:15]
    # cell_id_col = "cell_id"
    # batch_col = "Batch"
    # do_pca = TRUE
    # npcs = 20
    # theta = NULL
    # lambda = 1
    # sigma = 0.1
    # nclust = NULL
    # max_iter = 10
    # early_stop = TRUE
    # ncores = 1
    # plot_convergence = FALSE
    # return_object = TRUE
    # verbose = FALSE
    # seed = 42
    # # advanced parameters for Harmony
    # alpha = 0.2
    # tau = 0
    # block.size = 0.05
    # max.iter.cluster = 200
    # epsilon.cluster = 0.001
    # epsilon.harmony = 0.01
    
    check_packages_installed(c("harmony"))
    
    if (verbose) {
        
        if (do_pca) {
            n_steps <- 4
        } else {
            n_steps <- 3
        }
        cur_step <- 1
        
        message("Running Harmony")
        message(paste0("(", cur_step, "/", n_steps, ") Preparing data."))
    }
    
    if(length(use_cols) <= npcs){
        warning("There are more use_cols than npcs. Setting npcs to number of element in use_cols - 1")
        npcs <- length(use_cols) -1
        
    }
    
    # Harmony doesn't support running PCA anymore. It wants it to supply the PCs coordinate
    if (do_pca) {
        
        if (verbose) {
            cur_step <- cur_step +  1
            message(paste0("(", cur_step, "/", n_steps, ") Running PCA."))
        }
        
        cnt_mtx <- run.pca(
            dat = dat, 
            use.cols = use_cols, 
            scale = TRUE, 
            add.pca.col = TRUE, 
            pca.col.no = npcs, 
            scree.plot = FALSE, 
            variable.contribution = FALSE, 
            plot.individuals = FALSE, 
            plot.variables = FALSE, 
            plot.combined = FALSE)
        # keep just the PCs
        cnt_mtx <- cnt_mtx[, paste0("Dim.", seq(npcs))]
        # rename Dim. to PCs so the column headers are more explicit
        setnames(cnt_mtx, paste0("Dim.", seq(npcs)), paste0("PC_", seq(npcs)))
    } else {
        cnt_mtx <- dat[, use_cols, with = FALSE]
    }
    
    # create metadata containing cell id and batch column
    metadata <- dat[, c(cell_id_col, batch_col), with=FALSE]
    
    if (verbose) {
        if (verbose) {
            cur_step <- cur_step +  1
            message(paste0("(", cur_step, "/", n_steps, ") Running Harmony"))
        }
    }
    
    set.seed(seed)
    
    hrm_res <- harmony::RunHarmony(
        data_mat = cnt_mtx,
        meta_data = metadata,
        vars_use = batch_col,
        theta = theta,
        lambda = lambda,
        sigma = sigma,
        nclust = nclust,
        max_iter = max_iter,
        early_stop = early_stop,
        ncores = ncores,
        plot_convergence = plot_convergence,
        return_object = return_object,
        verbose = verbose,
        seed = seed,
        .options = harmony::harmony_options(
            alpha = alpha,
            tau = tau,
            block.size = block.size,
            max.iter.cluster = max.iter.cluster,
            epsilon.cluster = epsilon.cluster,
            epsilon.harmony = epsilon.harmony
        )
    )
    
    if (verbose) {
        if (verbose) {
            cur_step <- cur_step +  1
            message(paste0("(", cur_step, "/", n_steps, ") Constructing results to return"))
        }
    }
    
    if (return_object) {
        hrm_res_to_return <- data.table(t(hrm_res$Z_corr))
        hrm_res_metadata <- list("harmony_object" = hrm_res)
    } else {
        hrm_res_to_return <- data.table(hrm_res)
        hrm_res_metadata <- list()
    }
    
    if (do_pca) {
        names(hrm_res_to_return) <- paste0("PC_", seq(npcs))
    } else {
        names(hrm_res_to_return) <- use_cols
    }
    
    
    hrm_res_to_return[[cell_id_col]] <- dat[[cell_id_col]]
    hrm_res_to_return[[batch_col]] <- dat[[batch_col]]
    
    # So we can store the parameter as metadata.
    hrm_param <- data.table(
        vars_use = batch_col,
        theta = theta,
        lambda = lambda,
        sigma = sigma,
        nclust = nclust,
        max_iter = max_iter,
        early_stop = early_stop,
        ncores = ncores,
        plot_convergence = plot_convergence,
        return_object = return_object,
        verbose = verbose,
        seed = seed,
        alpha = alpha,
        tau = tau,
        block.size = block.size,
        max.iter.cluster = max.iter.cluster,
        epsilon.cluster = epsilon.cluster,
        epsilon.harmony = epsilon.harmony
    )
    hrm_param_names <- names(hrm_param)
    hrm_param <- transpose(hrm_param)
    setnames(hrm_param, "V1", "value")
    hrm_param[, parameter := hrm_param_names]
    setcolorder(hrm_param, c("parameter", "value"))
    
    # add the parameter into metadata
    hrm_res_metadata$harmony_parameter <- hrm_param
    
    for (meta_name in names(hrm_res_metadata)) {
        # meta_name = names(metadata)[1]
        attr(hrm_res_to_return, meta_name) <- hrm_res_metadata[[meta_name]]
        
    }
    
    return(hrm_res_to_return)
}




