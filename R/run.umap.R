#' Run UMAP
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
#'
#'
#' @examples
#' # Run UMAP on a subset of the  demonstration dataset
#'
#' cell.dat <- do.subsample(Spectre::demo.asinh, 10000) # Subsample the demo dataset to 10000 cells
#' cell.dat <- Spectre::run.umap(
#'   dat = cell.dat,
#'   use.cols = names(demo.asinh)[c(2:10)]
#' )
#' 
#' 
#' @author
#' Thomas Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#' Felix Marsh-Wakefield, \email{felix.marsh-wakefield@@sydney.edu.au}
#' @export
#' 
setGeneric("run.umap", function(
        dat,
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
        fast = FALSE,
        n_threads = 'auto',
        n_sgd_threads = 'auto',
        batch = TRUE,
        ...) {
    standardGeneric("run.umap")
    
})

#' @param data_source Character. The name of the data in Spectre object to 
#' run umap on.
#' Only used if dat is a Spectre object.
#' @param output_name Character. What name should the output
#' data be stored under in the Spectre object.
#' Only used if dat is a Spectre object.
#'
#' @exportMethod run.umap
#' @rdname run.umap
setMethod("run.umap", "Spectre", function(
        dat,
        use.cols,
        data_source,
        output_name,
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
        fast = FALSE,
        n_threads = 'auto',
        n_sgd_threads = 'auto',
        batch = TRUE) {
    
    # for testing
    # dat <- data.table(
    #     cell_id = paste0("cell_", seq(1, 1000)),
    #     marker1 = rnorm(1000, 1),
    #     marker2 = rnorm(1000, 2),
    #     marker3 = rnorm(1000, 3)
    # )
    # 
    # obj <- create.spectre.object(cell_id_col = "cell_id")
    # obj <- add.new.data(obj, dat, "test")
    # dat = obj
    # data_source = "test"
    # output_name = NULL
    # use.cols = c("marker1", "marker2", "marker3")
    # umap.x.name = "UMAP_X"
    # umap.y.name = "UMAP_Y"
    # umap.seed = 42
    # neighbours = 15
    # n_components = 2
    # metric = "euclidean"
    # n_epochs = 200
    # input = "data"
    # init = "spectral"
    # min_dist = 0.1
    # set_op_mix_ratio = 1
    # local_connectivity = 1
    # bandwidth = 1
    # alpha = 1
    # gamma = 1
    # negative_sample_rate = 5
    # a_gradient = NA
    # b_gradient = NA
    # spread = 1
    # transform_state = 42
    # knn.repeats = 1
    # verbose = TRUE
    # umap_learn_args = NA
    # # For Fast UMAP
    # fast = FALSE
    # n_threads = parallel::detectCores() - 1
    # n_sgd_threads = 'auto'
    # batch = TRUE
    
    df <- dat[[data_source]]
    
    umap_res <- run_actual_umap(
        dat = df,
        use.cols = use.cols,
        umap.x.name = umap.x.name,
        umap.y.name = umap.y.name,
        umap.seed = umap.seed,
        neighbours = neighbours,
        n_components = n_components,
        metric = metric,
        n_epochs = n_epochs,
        input = input,
        init = init,
        min_dist = min_dist,
        set_op_mix_ratio = set_op_mix_ratio,
        local_connectivity = local_connectivity,
        bandwidth = bandwidth,
        alpha = alpha,
        gamma = gamma,
        negative_sample_rate = negative_sample_rate,
        a_gradient = a_gradient,
        b_gradient = b_gradient,
        spread = spread,
        transform_state = transform_state,
        knn.repeats = knn.repeats,
        verbose = verbose,
        umap_learn_args = umap_learn_args,
        fast = fast,
        n_threads = n_threads,
        n_sgd_threads = n_sgd_threads,
        batch = batch
    )
    
    # TODO turn this into a function so we can use this in every other function.
    if (is.null(output_name)) {
        warning(paste("Appending umap coordinate to", data_source))
        
        if (umap.x.name %in% names(df)) {
            warning(paste(
                "Column", umap.x.name, "already present in", data_source, ". Overwriting it."
            ))
        }
        df[[umap.x.name]] <- umap_res[[umap.x.name]]
        
        if (umap.y.name %in% names(df)) {
            warning(paste(
                "Column", umap.y.name, "already present in", data_source, ". Overwriting it."
            ))
        }
        df[[umap.y.name]] <- umap_res[[umap.y.name]]
        
        dat <- add.new.data(
            spectre_obj = dat,
            dat = df,
            dat_name = data_source,
            metadata = dat@metadata[[data_source]]
        )
        
    } else {
        # just so the cell id column is first!
        cell_id_col <- dat@cell_id_col
        umap_dat <- data.table(cell_id = df[[cell_id_col]], umap_res)
        setnames(umap_dat, "cell_id", cell_id_col)
        
        dat <- add.new.data(
            spectre_obj = dat,
            dat = umap_dat,
            dat_name = output_name
        )
    }
    
    return(dat)
    
})

#' @rdname run.umap
#' @exportMethod run.umap
setMethod("run.umap", "data.table", function(
        dat,
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
        fast = FALSE,
        n_threads = 'auto',
        n_sgd_threads = 'auto',
        batch = TRUE) {
    
    umap_res <- run_actual_umap(
        dat = dat,
        use.cols = use.cols,
        umap.x.name = umap.x.name,
        umap.y.name = umap.y.name,
        umap.seed = umap.seed,
        neighbours = neighbours,
        n_components = n_components,
        metric = metric,
        n_epochs = n_epochs,
        input = input,
        init = init,
        min_dist = min_dist,
        set_op_mix_ratio = set_op_mix_ratio,
        local_connectivity = local_connectivity,
        bandwidth = bandwidth,
        alpha = alpha,
        gamma = gamma,
        negative_sample_rate = negative_sample_rate,
        a_gradient = a_gradient,
        b_gradient = b_gradient,
        spread = spread,
        transform_state = transform_state,
        knn.repeats = knn.repeats,
        verbose = verbose,
        umap_learn_args = umap_learn_args,
        fast = fast,
        n_threads = n_threads,
        n_sgd_threads = n_sgd_threads,
        batch = batch
    )
    
    if (umap.x.name %in% names(dat)) {
        warning(paste(
            "Column", umap.x.name, "already present in dat. Overwriting it."
        ))
    }
    dat[[umap.x.name]] <- umap_res[[umap.x.name]]
    
    if (umap.y.name %in% names(dat)) {
        warning(paste(
            "Column", umap.y.name, "already present in dat. Overwriting it."
        ))
    }
    dat[[umap.y.name]] <- umap_res[[umap.y.name]]
    
    # dat <- cbind(dat, umap_res)
    
    return(dat)
    
})

#' Internal function which actually run umap.
run_actual_umap <- function(
        dat,
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
        fast = FALSE,
        n_threads = 'auto',
        n_sgd_threads = 'auto',
        batch = TRUE) {
    ### Test data
    # dat <- iris
    # umap.seed <- 42
    # use.cols <- c(1:4)
    
    
    if(isFALSE(fast)){
        
        check_packages_installed(c("umap"))
        
        
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
        # res <- cbind(dat, umap.res) # Merge UMAP results with data
        # return(umap.res)
    } else {
        
        check_packages_installed(c("uwot"))
        
        if (n_threads == 'auto') {
            if ((!is.element("parallel", installed.packages()[, 1]))) {
                warning("Parallel package is not installed! Cannot infer optimal number of threads to use. Setting n_threads to 1.")
                n_threads <- 1
            } else {
                n_threads <- parallel::detectCores() - 1
            }
        }
        
        # Irritating. Can't peeps just settle on using either NA or NULL?
        if (is.na(a_gradient)) {
            a_gradient <- NULL
        }
        if (is.na(b_gradient)) {
            b_gradient <- NULL
        }
        
        set.seed(umap.seed)
        
        umap.res <- uwot::umap(
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
        colnames(umap.res) <- c(umap.x.name, umap.y.name)
        # res <- cbind(dat, umap.res)
    }
    
    umap.res <- data.table(umap.res)
    return(umap.res)
}


