#' Run the UMAP algorithm (using umap::umap())
#'
#' Method to run a UMAP dimensionality reduction algorithm.
#' A UMAP (uniform manifold approximation and projection) plot is a useful means to visualise data.
#' As it is a dimensionality reduction algorithm, some data will be lost.
#' It is good practice to validate any populations (namely through manual gating).
#' Uses the R package "umap" to calculate plots and "data.table" to handle data.
#'
#' @param dat NO DEFAULT. Input data.table or data.frame.
#' @param use.cols NO DEFAULT. Vector of column names or numbers for clustering.
#' @param umap.x.name DEFAULT = "UMAP_X". Character. Name of UMAP x-axis.
#' @param umap.y.name DEFAULT = "UMAP_Y". Character. Name of UMAP y-axis.
#' @param umap.seed DEFAULT = 42. Numeric. Seed value for reproducibility.
#' @param suffix DEFAULT = NULL
#'
#' @usage run.umap(dat, use.cols, umap.x.name = "UMAP_X", umap.y.name = "UMAP_Y", umap.seed = 42, suffix = NULL)
#'
#' @examples
#' # Run UMAP on a subset of the  demonstration dataset
#'
#' cell.dat <- do.subsample(Spectre::demo.asinh, 10000) # Subsample the demo dataset to 10000 cells
#' cell.dat <- Spectre::run.umap(dat = cell.dat,
#'                               use.cols = names(demo.asinh)[c(2:10)])
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
                     suffix = NULL
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

  dat.start <- dat

  dat <- dat[, use.cols, with = FALSE]

  res <- umap::umap(d = dat,
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
