#' run.phenograph - Run phenograph clustering
#'
#' This function allows you to perform phenograph clustering on a data.table. Uses 'Rphenograph' from https://github.com/JinmiaoChenLab/Rphenograph
#'
#' @usage run.phenograph()
#'
#' @param ref.dat NO DEFAULT. A data.table consisting of the 'refernece' data you will use to train the alignment algorithm

#' @return Returns a data.table with clustering results added
#'
#' @author
#' Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#' Felix Marsh-Wakefield, \email{felix.marsh-wakefield@@sydney.edu.au}
#' 
#' @references Ashhurst, T. M., et al. (2019). \url{https://www.ncbi.nlm.nih.gov/pubmed/31077106}
#'
#' @examples
#' \dontrun{
#' cell.dat <- Spectre::demo.clustered
#' cell.dat <- Spectre::do.subsample(cell.dat, 5000)
#' 
#' use.cols <- c("NK11_asinh", "CD3_asinh", "CD45_asinh", 
#' "Ly6G_asinh", "CD11b_asinh", "B220_asinh", "CD8a_asinh", 
#' "Ly6C_asinh", "CD4_asinh")
#' 
#' cell.dat <- run.phenograph(dat = cell.dat, use.cols = use.cols)
#' }
#' @import data.table
#'
#' @export

run.phenograph <- function(dat,
                           use.cols,
                           k = 45,
                           clust.name = "Phenograph_cluster") {

  ### Testing data

  # dat <- Spectre::demo.asinh
  # dat <- Spectre::do.subsample(dat, 5000)
  # use.cols <- names(dat)[c(11:19)]
  # k <- 45
  # clust.name <- 'Phenograph_cluster'

  ### Data setup

  dat.start <- dat
  dat <- dat[, use.cols, with = FALSE]
  dat <- as.matrix(dat)

  ### Clustering

  res <- Rphenograph::Rphenograph(dat, k = k)
  res <- factor(igraph::membership(res[[2]]))

  ### Attach result

  res <- as.data.table(res)
  names(res) <- clust.name
  dat.start <- cbind(dat.start, res)

  ### Return

  message("Phenograph clustering complete")
  return(dat.start)
}
