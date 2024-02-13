#' run.phenograph - Run phenograph clustering
#'
#' This function allows you to perform phenograph clustering on a data.table
#'
#' @usage run.phenograph()
#'
#' @param ref.dat NO DEFAULT. A data.table consisting of the 'refernece' data you will use to train the alignment algorithm

#' @return Returns a data.table with clustering results added
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references Ashhurst, T. M., et al. (2019). \url{https://www.ncbi.nlm.nih.gov/pubmed/31077106}
#'
#' @examples
#'
#' cell.dat <- Spectre::demo.asinh
#' cell.dat <- Spectre::do.subsample(cell.dat, 5000)
#'
#' use.cols <- names(cell.dat)[c(11:19)]
#'
#' cell.dat <- run.phenograph(dat = cell.dat, use.cols = use.cols)
#'
#'
#' @export

run.phenograph <- function(dat,
                           use.cols,
                           k = 45,
                           clust.name = "Phenograph_cluster") {

  ### Packages

  # require: Rphenograph, data.table, Rcpp
  # "JinmiaoChenLab/Rphenograph"
  check_packages_installed(c("Rphenograph", "Rcpp"))
  
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

  res <- Rphenograph(dat, k = k)
  res <- factor(membership(res[[2]]))

  ### Attach result

  res <- as.data.table(res)
  names(res) <- clust.name
  dat.start <- cbind(dat.start, res)

  ### Return

  message("Phenograph clustering complete")
  return(dat.start)
}
