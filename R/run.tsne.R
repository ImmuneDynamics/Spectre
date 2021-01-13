#' Run the tSNE algorithm (using Rtsne::Rtsne())
#'
#' Method to run a tSNE dimensionality reduction algorithm.
#' A tSNE (t-distributed stochastic neighbor embedding) plot is a useful means to visualise data.
#' As it is a dimensionality reduction algorithm, some data will be lost.
#' It is good practice to validate any populations (namely through manual gating).
#' Output data will be "tsne.res".
#' Uses the R package "Rtsne" to calculate plots.
#'
#' @param dat NO DEFAULT. data.frame.
#' @param use.cols NO DEFAULT. Vector of numbers, reflecting the columns to use for dimensionality reduction.
#' @param tsne.x.name DEFAULT = "tSNE_X". Character. Name of tSNE x-axis.
#' @param tsne.y.name DEFAULT = "tSNE_Y". Character. Name of tSNE y-axis.
#' @param tsne.seed DEFAULT = 42. Numeric. Seed value for reproducibility.
#' @param dims DEFAULT = 2. Number of dimensions for output results, either 2 or 3.
#' @param initial_dims DEFAULT = 50. Number of dimensions retained in initial PCA step.
#' @param perplexity DEFAULT = 30.
#' @param theta DEFAULT = 0.5. Use 0.5 for Barnes-Hut tSNE, 0.0 for exact tSNE (takes longer).
#' @param check_duplicates DEFAULT = FALSE.
#' @param pca DEFAULT = TRUE. Runs PCA prior to tSNE run.
#' @param max_iter DEFAULT = 1000. Maximum number of iterations.
#' @param verbose DEFAULT = TRUE.
#' @param is_distance DEFAULT = FALSE. Experimental, using X as a distance matrix.
#' @param Y_init DEFAULT = NULL. Recommend NULL for random initialisation.
#' @param stop_lying_iter DEFAULT = 250. Number of iterations of early exaggeration.
#' @param mom_switch_iter DEFAULT = 250. Number of iterations before increased momentum of spread.
#' @param momentum DEFAULT = 0.5. Initial momentum of spread.
#' @param final_momentum DEFAULT = 0.8. Momentum of spread at 'final_momentum'.
#' @param eta DEFAULT = 200. Learning rate.
#' @param exaggeration_factor DEFAULT = 12.0. Factor used during early exaggeration.
#'
#' @usage run.tsne(dat, use.cols, tsne.x.name, tsne.y.name, tsne.seed, dims, initial_dims, perplexity, theta, check_duplicates, pca, max_iter, verbose, is_distance, Y_init, stop_lying_iter, mom_switch_iter, momentum, final_momentum, eta, exaggeration_factor)
#'
#' @examples
#' # Run tSNE on a subset of the  demonstration dataset
#'
#' cell.dat <- do.subsample(Spectre::demo.asinh, 10000) # Subsample the demo dataset to 10000 cells
#' cell.dat <- Spectre::run.tsne(dat = cell.dat,
#'                               use.cols = names(demo.asinh)[c(2:10)])
#'
#' @author Felix Marsh-Wakefield, \email{felix.marsh-wakefield@@sydney.edu.au}
#' @export

run.tsne <- function(dat,
                 use.cols,
                 tsne.x.name = "tSNE_X",
                 tsne.y.name = "tSNE_Y",
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
  
  ## Set the seed
  set.seed(tsne.seed)

  ## Run tSNE
  tsne_out <- Rtsne::Rtsne(as.matrix(dat[, use.cols, with = FALSE]),
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
  colnames(tsne_out_Y) <- c(tsne.x.name, tsne.y.name)

  tsne.res <- cbind(dat, tsne_out_Y)

  assign("tsne.res", tsne.res, envir = globalenv())
}
