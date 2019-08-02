# umap

umap <- function(x,
                 use.cols = c(1:4),
                 umap.seed = 42
                 ){
  #x <- iris
  #umap.seed <- 42
  #use.cols <- c(1:4)
  
  custom.config <- umap.defaults
  custom.config$random_state <- umap.seed
  res <- umap(x[use.cols], custom.config)
  umap.res <- res$layout
  assign("umap.res", umap.res, envir = globalenv())
  
}



