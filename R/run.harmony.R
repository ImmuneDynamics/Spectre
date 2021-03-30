#' run.harmony - dun Harmony alignment on a data.table
#'
#' This function allows you to run the 'Harmony' data alignment algorithm on single cell or cytometry data stored in a data.table
#' 
#' @usage run.harmony()
#' 
#' @param dat NO DEFAULT. A data.table with all of the data you wish to align
#' @param align.cols NO default. The columns you wish to align. For cytometry data, this can be the markers themselves or principle components. For single-cell seq data, principle components are recommended.
#' @param batch.col NO default. The column that denotes the batch or dataset that each cell belongs to
#' @param append.name DEFAULT = '_aligned'. Text that will be appended to the new columns containing aligned data
#' @param do_pca DEFAULT = FALSE. Whether to perform PCA on input matrix. 
#' @param npcs If doing PCA on input matrix, number of PCs to compute. 
#' @param theta Diversity clustering penalty parameter. Specify for each
#'  variable in vars_use Default theta=2. theta=0 does not encourage any 
#'  diversity. Larger values of theta result in more diverse clusters. 
#' @param lambda Ridge regression penalty parameter. Specify for each variable
#'  in vars_use. 
#' Default lambda=1. Lambda must be strictly positive. Smaller values result 
#' in more aggressive correction. 
#' @param sigma Width of soft kmeans clusters. Default sigma=0.1. Sigma scales
#'  the distance from a cell to cluster centroids. Larger values of sigma 
#'  result in cells assigned to more clusters. Smaller values of sigma make 
#'  soft kmeans cluster approach hard clustering. 
#' @param nclust Number of clusters in model. nclust=1 equivalent to simple 
#' linear regression. 
#' @param tau Protection against overclustering small datasets with large ones.
#'  tau is the expected number of cells per cluster. 
#' @param block.size What proportion of cells to update during clustering.
#'  Between 0 to 1, default 0.05. Larger values may be faster but less accurate
#' @param max.iter.cluster Maximum number of rounds to run clustering at each 
#' round of Harmony. 
#' @param epsilon.cluster Convergence tolerance for clustering round of 
#' Harmony. Set to -Inf to never stop early. 
#' @param max.iter.harmony Maximum number of rounds to run Harmony. One round
#'  of Harmony involves one clustering and one correction step. 
#' @param epsilon.harmony Convergence tolerance for Harmony. Set to -Inf to
#'  never stop early. 
#' @param plot_convergence Whether to print the convergence plot of the 
#' clustering objective function. TRUE to plot, FALSE to suppress. This can be
#'  useful for debugging. 
#' @param return_object (Advanced Usage) Whether to return the Harmony object 
#' or only the corrected PCA embeddings. 
#' @param verbose DEFAULT = FALSE. Whether to print progress messages. TRUE to print, 
#' FALSE to suppress.
#' @param reference_values (Advanced Usage) Defines reference dataset(s). 
#' Cells that have batch variables values matching reference_values will not 
#' be moved.  
#' @param cluster_prior (Advanced Usage) Provides user defined clusters for 
#' cluster initialization. If the number of provided clusters C is less than K, 
#' Harmony will initialize K-C clusters with kmeans. C cannot exceed K.  
#' 
#' 
#' @return Returns a data.table with aligned data added in new columns.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @examples
#' cell.dat <- run.harmony()
#'
#' @import data.table
#'
#' @export

run.harmony <- function(dat,
                        align.cols,
                        batch.col,
                        append.name = '_aligned',
                        do_pca = FALSE, 
                        npcs = 20, 
                        theta = NULL, 
                        lambda = NULL, 
                        sigma = 0.1, 
                        nclust = NULL, 
                        tau = 0, 
                        block.size = 0.05, 
                        max.iter.harmony = 10, 
                        max.iter.cluster = 200, 
                        epsilon.cluster = 1e-5, 
                        epsilon.harmony = 1e-4, 
                        plot_convergence = FALSE, 
                        return_object = FALSE, 
                        verbose = FALSE, 
                        reference_values = NULL, 
                        cluster_prior = NULL){
  
  ### Packages
  
      if(!is.element('Spectre', installed.packages()[,1])) stop('Spectre is required but not installed')
      if(!is.element('data.table', installed.packages()[,1])) stop('data.table is required but not installed')
      if(!is.element('harmony', installed.packages()[,1])) stop('harmony is required but not installed. You can install harmony by running devtools::install_github("immunogenomics/harmony")')
      
  ### Require packages
  
      require(Spectre)
      require(data.table)
      require(harmony)
  
  ### Data prep
  
      message("run.harmony - preparing data (1/4)")
  
      start.dat <- dat

      dat <- dat[,align.cols, with = FALSE]
      nms <- names(dat)
      dat <- as.matrix(dat)
      
      meta <- data.table()
      meta$CellID <- c(1:nrow(dat))
      meta$CellID <- as.character(meta$CellID)
      meta <- cbind(meta, start.dat[,batch.col, with = FALSE])
      meta <- as_tibble(meta)
      
  ### Run harmony
      
      message("run.harmony - running harmony (2/4)")
      
      hrm.res <- harmony::HarmonyMatrix(data_mat = dat, 
                                        meta_data = meta, 
                                        vars_use = batch.col, 
                                        do_pca = do_pca, 
                                        npcs = npcs, 
                                        theta = theta, 
                                        lambda = lambda, 
                                        sigma = sigma, 
                                        nclust = nclust, 
                                        tau = tau, 
                                        block.size = block.size, 
                                        max.iter.harmony = max.iter.harmony, 
                                        max.iter.cluster = max.iter.cluster, 
                                        epsilon.cluster = epsilon.cluster, 
                                        epsilon.harmony = epsilon.harmony, 
                                        plot_convergence = plot_convergence, 
                                        return_object = return_object, 
                                        verbose = verbose, 
                                        reference_values = reference_values, 
                                        cluster_prior = cluster_prior)
      
  ### Final preparation and return
      
      message("run.harmony - harmony complete, finalising data (3/4)")
      
      hrm.res <- as.data.table(hrm.res)
      names(hrm.res) <- paste0(names(hrm.res), append.name)
      hrm.res
      
      final.res <- cbind(start.dat, hrm.res)
      message("run.harmony - returning data (4/4)")
      
      return(final.res)
}
