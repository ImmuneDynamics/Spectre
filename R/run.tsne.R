#' run.tsne - ...
#'
#' @param dat data.frame. No default.
#' @param use.cols Vector of numbers, reflecting the columns to use for clustering. No default.
#' @param tsne.seed Numeric. Seed value for reproducibility. Default = 42.
#' @param dims Number of dimensions for output results, either 2 or 3. Default = 2.
#' @param initial_dims Number of dimensions retained in initial PCA step. Default = 50.
#' @param perplexity Default = 30.
#' @param theta Use 0.5 for Barnes-Hut tSNE, 0.0 for exact tSNE (takes longer). Default = 0.5.
#' @param check_duplicates Default = FALSE.
#' @param pca Runs PCA prior to tSNE run. Default = TRUE.
#' @param max_iter Maximum number of iterations. Default = 1000.
#' @param verbose Default = TRUE.
#' @param is_distance Experimental, using X as a distance matrix. Default = FALSE.
#' @param Y_init Recommend NULL for random initialisation. Default = NULL.
#' @param stop_lying_iter Number of iterations of early exaggeration. Default = 250.
#' @param mom_switch_iter Number of iterations before increased momentum of spread. Default = 250.
#' @param momentum Initial momentum of spread. Default = 0.5.
#' @param final_momentum Momentum of spread at 'final_momentum'. Default = 0.8.
#' @param eta Learning rate. Default = 200.
#' @param exaggeration_factor Factor used during early exaggeration. Default = 12.0.
#'
#' This function runs tSNE on a dataframe with cells (rows) vs markers (columns), and returns 'res' with result columns. Uses the Rtsne package. For more information on parameter choices and effects, check out https://distill.pub/2016/misread-tsne/.
#' 
#' @usage run.tsne(dat, use.cols, tsne.seed, dims, initial_dims, perplexity, theta, check_duplicates, pca, max_iter, verbose, is_distance, Y_init, stop_lying_iter, mom_switch_iter, momentum, final_momentum, eta, exaggeration_factor, ...)
#'
#' @export

run.tsne <- function(dat,
                 use.cols,
                 tsne.seed = 42,
                 dims = 2,
                 initial_dims = 50,
                 perplexity = 30,
                 theta = 0.5,
                 check_duplicates = FALSE,
                 pca = TRUE,
                 max_iter = 1000,
                 verbose = TRUE,
                 is_distance = FALSE,
                 Y_init = NULL,
                 stop_lying_iter = 250,
                 mom_switch_iter = 250,
                 momentum = 0.5,
                 final_momentum = 0.8,
                 eta = 200,
                 exaggeration_factor = 12.0
                 ){
  
  ## Check that necessary packages are installed
  if(!is.element('Rtsne', installed.packages()[,1])) stop('Rtsne is required but not installed')
  
  ## Require packages
  require(Rtsne)
  
  ## Run tSNE
  tsne_out <- Rtsne::Rtsne(as.matrix(dat[use.cols]),
        dims = dims,
        initial_dims = initial_dims,
        perplexity = perplexity,
        theta = theta,
        check_duplicates = check_duplicates,
        pca = pca,
        max_iter = max_iter,
        verbose = verbose,
        is_distance = is_distance,
        Y_init = Y_init,
        stop_lying_iter = stop_lying_iter,
        mom_switch_iter = mom_switch_iter,
        momentum = momentum,
        final_momentum = final_momentum,
        eta = eta,
        exaggeration_factor = exaggeration_factor
  )
  
  tsne_out_Y <- tsne_out$Y
  colnames(tsne_out_Y) <- c(paste0("tSNE", "_", tsne.seed, "_", "X"), paste0("tSNE", "_", tsne.seed, "_", "Y"))
  
  tsne.res <- cbind(dat, tsne_out_Y)
  
  assign("tsne.res", tsne.res, envir = globalenv())
}
