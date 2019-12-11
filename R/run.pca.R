#' run.pca - ...
#'
#' @usage run.pca(x, use.cols, ...)
#'
#' @param x data.frame. No default.
#' @param use.cols Vector of numbers, reflecting the columns to use for clustering. No default.
#' @param cor a logical value indicating whether the calculation should use the correlation matrix or the covariance matrix. (The correlation matrix can only be used if there are no constant variables.). Default = TRUE.
#' @param scores a logical value indicating whether the score on each principal component should be calculated. Default = TRUE.
#' @param scree.plot option to create scree plots. Note this will require the input of an elbow point during run. Will save generated scree plot. Default = TRUE.
#' @param loading.plot option to create scree plots. Will save generated loading plot. Default = TRUE.
#'
#'This function runs a principal component analysis (PCA) on a dataframe with cells (rows) vs markers (columns), returning chosen figures. Uses the base R package "stats" for PCA, "factoextra" for scree and loading plots.
#'
#' @export

run.pca <- function(x,
                 use.cols,
                 cor = TRUE,
                 scores = TRUE,
                 scree.plot = TRUE,
                 loading.plot = TRUE
                 ){
  
  pca_out <- stats::princomp(as.matrix(x[use.cols]),
                             cor = cor,
                             scores = scores,
                             scree.plot = scree.plot,
                             loading.plot = loading.plot
                             )
  
  if (scree.plot == TRUE) {
    scree_plot <- factoextra::fviz_eig(pca_out, addlabels = TRUE) #creates scree plot; addlabels adds % to plot
    print(scree_plot)
    
    elbow.point <- readline("Type in the elbow point based on the scree plot. Must be positive integer. ")
    elbow.point <- as.numeric(elbow.point)
    
    # Saves scree plot
    ggsave(scree_plot, filename = "Scree plot.pdf")
  }
  
  if (loading.plot == TRUE) {
    loading_plot <- factoextra::fviz_pca_var(pca_out, col.var = "contrib",
                                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
    
    # Saves loadings plot
    ggsave(loading_plot, filename = "Loading plot.pdf")
    
    data.loadings <- unclass(pca_out$loadings) #shows the loadings values for everything!
    fwrite(x = as.data.frame(data.loadings), file = "loadings.csv", row.names = TRUE)
  }
  
}

