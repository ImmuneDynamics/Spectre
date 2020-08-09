#' run.umap - ...
#'
#' @param dat NO DEFAULT. Input data.table or data.frame.
#' @param use.cols NO DEFAULT. Vector of column names or numbers for clustering.
#' @param umap.x.name DEFAULT = "UMAP_X". Character. Name of UMAP x-axis.
#' @param umap.y.name DEFAULT = "UMAP_Y". Character. Name of UMAP y-axis.
#' @param umap.seed DEFAULT = 42. Numeric. Seed value for reproducibility.
#' @param suffix DEFAULT = blank.
#'
#' This function runs UMAP on a dataframe with cells (rows) vs markers (columns), and returns 'res' with result columns.
#'
#' @usage run.umap(dat, use.cols, umap.x.name, umap.y.name, umap.seed, suffix, ...)
#'
#' @import data.table
#'
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

  dat <- dat[,..use.cols]

  res <- umap::umap(d = dat,
              condif = custom.config)

  umap.res <- res$layout
  head(umap.res)

  umap.res <- as.data.frame(umap.res)
  head(umap.res)

  if(is.null(suffix) == TRUE){ # if suffix = NULL
    names(umap.res)[names(umap.res) == "V1"] <- umap.x.name
    names(umap.res)[names(umap.res) == "V2"] <- umap.y.name
  }

  if(is.null(suffix) == FALSE){ # if suffix = anything other than NULL
    names(umap.res)[names(umap.res) == "V1"] <- paste0("UMAP_X", "_", suffix)
    names(umap.res)[names(umap.res) == "V2"] <- paste0("UMAP_Y", "_", suffix)
  }

  #assign("umap.res", umap.res, envir = globalenv())
  res <- cbind(dat.start, umap.res) # Merge UMAP results with data
  return(res)
}
