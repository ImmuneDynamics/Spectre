#' V2 version of run umap. Experimental. Use at your own risk.
#'
#' @param ... Parameters to be passed to uwot::umap. 
#' See their vignette for more info.
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
    
    ## Require packages
    require(uwot)
    require(data.table)
    
    set.seed(umap.seed)
    
    if (! "data.frame" %in% class(dat))
        stop("dat must be of type data.frame or data.table")
    
    dat.copy <- copy(dat)
    
    # for testing only
    # dat.copy.sub <- dat.copy[, ..use.cols][sample(.N, 1000)]
    dat.umap <- uwot::umap(dat.copy[, ..use.cols], ...)
    
    # Preparing data to return
    colnames(dat.umap) <- c(umap.x.name, umap.y.name)
    dat.copy <- cbind(dat.copy, dat.umap)
    
    return(dat.copy)
}
