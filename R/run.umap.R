#' run.umap - ...
#'
#' @usage run.umap(x, ...)
#'
#' @param x NO DEFAULT. Input data.table or data.frame.
#' @param use.cols NO DEFAULT. Vector of column names or numbers for clustering.
#' @param umap.seed DEFAULT = 42. Numeric. Seed value for reproducibility.
#' @param suffix DEFAULT = blank.
#'
#'This function runs UMAP on a dataframe with cells (rows) vs markers (columns), and returns 'res' with result columns.
#'
#' @export

run.umap <- function(x,
                     use.cols,
                     umap.seed = 42,
                     suffix = NULL suffix = "test"
                     )
{

  ### Test data
      # x <- iris
      # umap.seed <- 42
      # use.cols <- c(1:4)

  ###
  custom.config <- umap.defaults
  custom.config$random_state <- umap.seed

  res <- umap(d = x[use.cols],
              condif = custom.config)

  umap.res <- res$layout
  head(umap.res)

  umap.res <- as.data.frame(umap.res)
  head(umap.res)

  if(is.null(suffix) == TRUE){ # if suffix = NULL
    names(umap.res)[names(umap.res) == "V1"] <- paste0("UMAP_X")
    names(umap.res)[names(umap.res) == "V2"] <- paste0("UMAP_Y")
  }

  if(is.null(suffix) == FALSE){ # if suffix = anything other than NULL
    names(umap.res)[names(umap.res) == "V1"] <- paste0("UMAP_X", "_", suffix)
    names(umap.res)[names(umap.res) == "V2"] <- paste0("UMAP_Y", "_", suffix)
  }

  assign("umap.res", umap.res, envir = globalenv())
}
