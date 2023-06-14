#' DEPRECATED Faster version of run umap.
#'
#' 15/6/2023: This function has been merged into run.umap. 
#' To use, please run run.umap with parameter fast set to TRUE.
#' 
#' This umap function runs the implementation provided by the uwot package.
#' *Note, by default, the calculation will be run in parallel.*
#' To be more specific, `n_threads` and `n_sgd_threads` will be automatically set to maximum number
#' of cores in your computer - 1.
#' However, you can override this by setting both to NULL to run it in serial.
#' See uwot::umap vignette for more information.
#' This umap implementation is ***much much*** faster than the one provided by the
#' umap package.
#' Once this has been long established, we will deprecate the old one.
#'
#' @param dat Data.table containing your cytometry data.
#' @param use.cols Markers to be used to calculate UMAP stored in a vector.
#' @param umap.x.name DEFAULT = "UMAP_X". Character. Name of UMAP x-axis.
#' @param umap.y.name DEFAULT = "UMAP_Y". Character. Name of UMAP y-axis.
#' @param umap.seed DEFAULT = 42. Numeric. Seed value for reproducibility.
#' @param n_threads DEFAULT `detectCores()-1`. Numeric. Number of threads to use (except during stochastic gradient descent). For nearest neighbor search, only applies if `nn_method = "annoy"`. If `n_threads > 1`, then the Annoy index will be temporarily written to disk in the location determined by tempfile.
#' @param n_sgd_threads DEFAULT "auto". Number of threads to use during stochastic gradient descent. If set to > 1, then be aware that if `batch = FALSE`, results will not be reproducible, even if `set.seed` is called with a fixed seed before running. Set to "auto" to use the same value as n_threads.
#' @param batch DEFAULT TRUE. If set to TRUE, then embedding coordinates are updated at the end of each epoch rather than during the epoch. In batch mode, results are reproducible with a fixed random seed even with n_sgd_threads > 1, at the cost of a slightly higher memory use. You may also have to modify learning_rate and increase n_epochs, so whether this provides a speed increase over the single-threaded optimization is likely to be dataset and hardware-dependent.
#' @param ... Parameters to be passed to uwot::umap. See their vignette for more info.
#'
#' @details
#' As of Version 0.1.11, it is now possible to get reproducible results (for a given value of set.seed) when running the optimization step with multiple threads (`n_sgd_threads` greater than 1).
#' You may need to increase n_epochs to get similar levels of convergence.
#' To run in this mode, set batch = TRUE.
#'
#' @import data.table
#' @importFrom parallel detectCores
#' @importFrom uwot umap
#'
#' @author Givanna Putri
#' @export

run.fast.umap <- function(dat,
                          use.cols,
                          umap.x.name = "UMAP_X",
                          umap.y.name = "UMAP_Y",
                          umap.seed = 42,
                          n_threads = detectCores() - 1,
                          n_sgd_threads = "auto",
                          batch = TRUE,
                          ...) {
    # Check input
    if (!"data.frame" %in% class(dat)) {
        stop("dat must be of type data.frame or data.table")
    }

    # Check package dependencies
    check_packages_installed("parallel")
    check_packages_installed("uwot")

    # for testing only
    # dat.copy.sub <- dat.copy[, ..use.cols][sample(.N, 1000)]

    set.seed(umap.seed)
    dat.umap <- uwot::umap(
        X = dat[, use.cols, with = FALSE],
        n_threads = n_threads,
        n_sgd_threads = n_sgd_threads,
        batch = batch,
        ...
    )

    # Preparing data to return
    colnames(dat.umap) <- c(umap.x.name, umap.y.name)
    dat <- cbind(dat, dat.umap)

    return(dat)
}
