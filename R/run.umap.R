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
#' cell.dat <- Spectre::run.umap(dat = cell.dat,
#'                               use.cols = c("NK11_asinh", "CD3_asinh", 
#'                               "CD45_asinh", "Ly6G_asinh", "CD11b_asinh", 
#'                               "B220_asinh", "CD8a_asinh", "Ly6C_asinh", 
#'                               "CD4_asinh"))
#' @author
#' Thomas Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#' Felix Marsh-Wakefield, \email{felix.marsh-wakefield@@sydney.edu.au}
#' @export

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
                     umap_learn_args = NA
                     )
{

  ### Test data
      # dat <- iris
      # umap.seed <- 42
      # use.cols <- c(1:4)

  ## Check that necessary packages are installed
  if(!is.element('umap', installed.packages()[,1])) stop('umap is required but not installed')
  if(!is.element('data.table', installed.packages()[,1])) stop('data.table is required but not installed')

  ## Require packages
  require(umap)
  require(data.table)

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

  dat.start <- data.table(dat)
  dat.bk <- data.table(dat)

  dat.bk <- dat.bk[, ..use.cols]

  res <- umap::umap(d = dat.bk,
              config = custom.config)

  umap.res <- res$layout
  head(umap.res)

  umap.res <- as.data.frame(umap.res)
  head(umap.res)

  names(umap.res)[names(umap.res) == "V1"] <- umap.x.name
  names(umap.res)[names(umap.res) == "V2"] <- umap.y.name

  #assign("umap.res", umap.res, envir = globalenv())
  res <- cbind(dat.start, umap.res) # Merge UMAP results with data
  return(res)
}
