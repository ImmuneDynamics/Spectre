#' run.umap - ...
#'
#' @usage run.umap(x, ...)
#'
#' @param x data.frame. No default.
#' @param use.cols Vector of numbers, reflecting the columns to use for clustering. No default.
#' @param umap.seed Numeric. Seed value for reproducibility. No default.
#'
#'This function runs UMAP on a dataframe with cells (rows) vs markers (columns), and returns 'res' with result columns.
#'
#' @export

run.umap <- function(x,
                 use.cols,
                 umap.seed
                 ){
  #x <- iris
  #umap.seed <- 42
  #use.cols <- c(1:4)

  custom.config <- umap.defaults
  custom.config$random_state <- umap.seed

  res <- umap(d = x[use.cols],
              condif = custom.config)

  umap.res <- res$layout
  head(umap.res)

  umap.res <- as.data.frame(umap.res)
  head(umap.res)

  names(umap.res)[names(umap.res) == "V1"] <- paste0("UMAP", "_", umap.seed, "_", "X")
  names(umap.res)[names(umap.res) == "V2"] <- paste0("UMAP", "_", umap.seed, "_", "Y")

  assign("umap.res", umap.res, envir = globalenv())

}
