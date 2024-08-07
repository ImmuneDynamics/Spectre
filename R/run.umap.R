#' Run the UMAP algorithm (using umap::umap())
#'
#' Method to run a UMAP dimensionality reduction algorithm.
#' A UMAP (uniform manifold approximation and projection) plot is a useful means to visualise data.
#' As it is a dimensionality reduction algorithm, some data will be lost.
#' It is good practice to validate any populations (namely through manual gating).
#' For more information on parameter choices, see ?umap::umap.defaults.
#' Uses the R package "umap" to calculate plots and "data.table" to handle data.
#'
#' @param dat NO DEFAULT. Input data.table or data.frame.
#' @param use.cols NO DEFAULT. Vector of column names or numbers for clustering.
#' @param umap.x.name DEFAULT = "UMAP_X". Character. Name of UMAP x-axis.
#' @param umap.y.name DEFAULT = "UMAP_Y". Character. Name of UMAP y-axis.
#' @param umap.seed DEFAULT = 42. Numeric. Seed value for reproducibility.
#' @param neighbours DEFAULT = 15. Numeric. Number of nearest neighbours.
#' @param n_components DEFAULT = 2. Numeric. Number of dimensions for output results.
#' @param metric DEFAULT = "euclidean". Character or function. Determines how distances between data points are computed. Can also be "manhattan".
#' @param n_epochs DEFAULT = 200. Numeric. Number of iterations performed during layout optimisation.
#' @param input DEFAULT = "data". Character. Determines whether primary input argument is a data or distance matrix. Can also be "dist".
#' @param init DEFAULT = "spectral". Character or matrix. Deafult "spectral" computes an initial embedding using eigenvectors of the connectivity graph matrix. Can also use "random" (creates an initial layout based on random coordinates).
#' @param min_dist DEFAULT = 0.1. Numeric. Determines how close points appear in final layout.
#' @param set_op_mix_ratio DEFAULT = 1. Numeric in range [0,1]. Determines who the knn-graph is used to create a fuzzy simplicial graph.
#' @param local_connectivity DEFAULT = 1. Numeric. Used during construction of fuzzy simplicial set.
#' @param bandwidth DEFAULT = 1. Numeric. Used during construction of fuzzy simplicial set.
#' @param alpha DEFAULT = 1. Numeric. Initial value of "learning rate" of layout optimisation.
#' @param gamma DEFAULT = 1. Numeric. Together with alpha, it determines the learning rate of layout optimisation.
#' @param negative_sample_rate DEFAULT = 5. Numeric. Determines how many non-neighbour points are used per point and per iteration during layout optimisation.
#' @param a_gradient DEFAULT = NA. Numeric. Contributes to gradient calculations during layout optimisation. When left at NA, a suitable value will be estimated automatically.
#' @param b_gradient DEFAULT = NA. Numeric. Contributes to gradient calculations during layout optimisation. When left at NA, a suitable value will be estimated automatically.
#' @param spread DEFAULT = 1. Numeric. Used during automatic estimation of a_gradient/b_gradient parameters.
#' @param transform_state DEFAULT = 42. Numeric. Seed for random number generation used during predict().
#' @param knn.repeats DEFAULT = 1. Numeric. Number of times to restart knn search.
#' @param verbose DEFAULT = TRUE. Logical. Determines whether to show progress messages.
#' @param umap_learn_args DEFAULT = NA. Vector. Vector of arguments to python package umap-learn.
#' @param n_threads DEFAULT "auto". Numeric. Number of threads to use (except during stochastic gradient descent). For nearest neighbor search, only applies if `nn_method = "annoy"`. If `n_threads > 1`, then the Annoy index will be temporarily written to disk in the location determined by tempfile. The default "auto" option will automatically set this to the maximum number of threads in the computer - 1.
#' @param n_sgd_threads DEFAULT "auto". Number of threads to use during stochastic gradient descent. If set to > 1, then be aware that if `batch = FALSE`, results will not be reproducible, even if `set.seed` is called with a fixed seed before running. Set to "auto" to use the same value as n_threads.
#' @param batch DEFAULT TRUE. If set to TRUE, then embedding coordinates are updated at the end of each epoch rather than during the epoch. In batch mode, results are reproducible with a fixed random seed even with n_sgd_threads > 1, at the cost of a slightly higher memory use. You may also have to modify learning_rate and increase n_epochs, so whether this provides a speed increase over the single-threaded optimization is likely to be dataset and hardware-dependent.
#' @param fast DEFAULT TRUE Whether to run uwot implementation of UMAP which is much faster.
#'
#' @usage run.umap(dat, use.cols, umap.x.name = "UMAP_X", 
#' umap.y.name = "UMAP_Y", umap.seed = 42, neighbours = 15, 
#' n_components = 2, metric = "euclidean", n_epochs = 200, 
#' input = "data", init = "spectral", min_dist = 0.1, 
#' set_op_mix_ratio = 1, local_connectivity = 1, bandwidth = 1, 
#' alpha = 1, gamma = 1, negative_sample_rate = 5, a_gradient = NA, 
#' b_gradient = NA, spread = 1, transform_state = 42, 
#' knn.repeats = 1, verbose = TRUE, umap_learn_args = NA)
#'
#' @examples
#' # Run UMAP on a subset of the  demonstration dataset
#'
#' cell.dat <- do.subsample(Spectre::demo.clustered, 10000) # Subsample the demo dataset to 10000 cells
#' cell.dat$UMAP_X <- NULL
#' cell.dat$UMAP_Y <- NULL
#' 
#' cell.dat <- Spectre::run.umap(dat = cell.dat,
#'                               use.cols = c("NK11_asinh", "CD3_asinh", 
#'                               "CD45_asinh", "Ly6G_asinh", "CD11b_asinh", 
#'                               "B220_asinh", "CD8a_asinh", "Ly6C_asinh", 
#'                               "CD4_asinh"))
#' @author
#' Thomas Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#' Felix Marsh-Wakefield, \email{felix.marsh-wakefield@@sydney.edu.au}
#' 
#' @import data.table
#' 
#' @export
#' 

run.umap <- function(dat,
                     use.cols,
                     umap.x.name = "UMAP_X",
                     umap.y.name = "UMAP_Y",
                     umap.seed = 42,
                     neighbours = 15,
                     n_components = 2,
                     metric = "euclidean",
                     n_epochs = 200,
                     input = "data",
                     init = "spectral",
                     min_dist = 0.1,
                     set_op_mix_ratio = 1,
                     local_connectivity = 1,
                     bandwidth = 1,
                     alpha = 1,
                     gamma = 1,
                     negative_sample_rate = 5,
                     a_gradient = NA,
                     b_gradient = NA,
                     spread = 1,
                     transform_state = 42,
                     knn.repeats = 1,
                     verbose = TRUE,
                     umap_learn_args = NA,
                     
                     # For Fast UMAP
                     fast = TRUE,
                     n_threads = 'auto',
                     n_sgd_threads = 'auto',
                     batch = TRUE) {
    
    ### Test data
    # dat <- iris
    # umap.seed <- 42
    # use.cols <- c(1:4)
    
    ## Check that necessary packages are installed
    # check_packages_installed(c("data.table"))
    # require(data.table)
    
    if (!"data.frame" %in% class(dat)) {
        stop("dat must be of type data.frame or data.table")
    }
    
    if(fast){
        if (!is.element("parallel", installed.packages()[, 1])){
            message("For 'fast' UMAP, parallel is required but not installed. Switching to slow UMAP")
            fast <- FALSE
        }
    }
    
    
    if(fast == FALSE){
        
        check_packages_installed(c("umap"))
        require(umap)
        
        ###
        custom.config <- umap::umap.defaults
        custom.config$random_state <- umap.seed
        
        custom.config$n_neighbors <- neighbours
        custom.config$n_components <- n_components
        custom.config$metric <- metric
        custom.config$n_epochs <- n_epochs
        custom.config$input <- input
        custom.config$init <- init
        custom.config$min_dist <- min_dist
        custom.config$set_op_mix_ratio <- set_op_mix_ratio
        custom.config$local_connectivity <- local_connectivity
        custom.config$bandwidth <- bandwidth
        custom.config$alpha <- alpha
        custom.config$gamma <- gamma
        custom.config$negative_sample_rate <- negative_sample_rate
        custom.config$a <- a_gradient
        custom.config$b <- b_gradient
        custom.config$spread <- spread
        custom.config$transform_state <- transform_state
        custom.config$knn.repeats <- knn.repeats
        custom.config$verbose <- verbose
        custom.config$umap_learn_args <- umap_learn_args
        
        ###
        
        res <- umap::umap(
            d = dat[, ..use.cols],
            config = custom.config
        )
        
        umap.res <- res$layout
        umap.res <- as.data.frame(umap.res)
        names(umap.res)[names(umap.res) == "V1"] <- umap.x.name
        names(umap.res)[names(umap.res) == "V2"] <- umap.y.name
        
        # assign("umap.res", umap.res, envir = globalenv())
        res <- cbind(dat, umap.res) # Merge UMAP results with data
        return(res)
    }
    
    if (fast){
        
        # Irritating. Can't peeps just settle on using either NA or NULL?
        if (is.na(a_gradient)) {
            a_gradient <- NULL
        }
        if (is.na(b_gradient)) {
            b_gradient <- NULL
        }
        
        set.seed(umap.seed)
        
        # set number of threads
        # the following could easy just check for is not numeric, but whatever. 
        # no need to do anything if n_threads is numeric as it will automatically 
        # be passed on.
        if (n_threads == 'auto' || !is.numeric(n_threads)) {
            n_threads <- parallel::detectCores() - 1
        }
        
        dat.umap <- uwot::umap(
            X = dat[, use.cols, with = FALSE],
            n_threads = n_threads, 
            n_sgd_threads = n_sgd_threads,
            batch = batch,
            n_neighbors = neighbours,
            n_components = n_components,
            metric = metric,
            n_epochs = n_epochs,
            init = init,
            min_dist = min_dist,
            set_op_mix_ratio = set_op_mix_ratio,
            local_connectivity = local_connectivity,
            bandwidth = bandwidth,
            negative_sample_rate = negative_sample_rate,
            a = a_gradient,
            b = b_gradient,
            spread = spread,
            verbose = verbose
        )
        
        # Preparing data to return
        colnames(dat.umap) <- c(umap.x.name, umap.y.name)
        res <- cbind(dat, dat.umap)
        
        return(res)
    }
    
}
