# # Spectre::umap

umap <- function(x,
                 use.cols,
                 umap.seed,
                 ){

  ## Test values
      #x <- iris
      #umap.seed <- 42
      #use.cols <- c(1:4)

  ## Setup parameters
  custom.config <- umap.defaults
  custom.config$random_state <- umap.seed

  ## Run UMAP
  res <- umap(x[use.cols], custom.config)
  umap.res <- res$layout

  ## Modify result 'UMAP column' names
  names(umap.res)[names(umap.res) == "1"] <- paste0("UMAP", "_", umap.seed, "_", "1")
  names(umap.res)[names(umap.res) == "2"] <- paste0("UMAP", "_", umap.seed, "_", "1")

  ## Save UMAP result to global env
  assign("umap.res", umap.res, envir = globalenv())

}

