#' V2 version of run umap. 
#' 
#' Experimental feature. Use at your own risk!
#' *Note, by default, uwot will be run in parallel by default.*
#' To be more specific, `n_threads` will be automatically set to maximum number
#' of cores in your computer - 1. 
#' However, you can override this by just passing the n_threads parameter.
#' Just a reminder, by default, we do not modify the `n_sgd_threads` parameter.
#' It will be set to 1 unless you explicitly modify it. 
#' Thus by default, stochastic gradient descent will be calculated serially, which
#' is better for result reproducibility as the same seed passed as `umap.seed` will be used.
#' See uwot::umap vignette for more information
#'
#' @param ... Parameters to be passed to uwot::umap. 
#' See their vignette for more info.
#' 
#' @author Givanna Putri
#' @export run.umap.v2

run.umap.v2 <- function(dat,
                        use.cols,
                        umap.x.name = "UMAP_X",
                        umap.y.name = "UMAP_Y",
                        umap.seed = 42,
                        ...) {
    
    
    ## Check that necessary packages are installed
    necessary_packages <- c("uwot", "data.table")
    check_packages_installed(necessary_packages)
    
    # Check input
    if (! "data.frame" %in% class(dat))
        stop("dat must be of type data.frame or data.table")
    
    ## Require packages
    require(uwot)
    require(data.table)
    
    set.seed(umap.seed)
    
    dat.copy <- copy(dat)
    
    # for testing only
    # dat.copy.sub <- dat.copy[, ..use.cols][sample(.N, 1000)]
    
    dat.umap <- uwot::umap(
        X=dat.copy[, ..use.cols],
        ...)
    
    # Preparing data to return
    colnames(dat.umap) <- c(umap.x.name, umap.y.name)
    dat.copy <- cbind(dat.copy, dat.umap)
    
    return(dat.copy)
}
