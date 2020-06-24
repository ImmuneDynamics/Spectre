#' Run the ChronoClust algorithm
#' 
#' Method to run ChronoClust clustering algorithm.
#' For this to work properly, a Python session must already be embedded within the currently running R session.
#' This is because ChronoClust is originally written in Python 3. 
#' Spectre includes a function \code{\link{run.prepare.chronoclust}} that uses Reticulate to prepare and embed Python session within R session.
#' 
#' @seealso \url{https://github.com/ghar1821/Chronoclust} for ChronoClust's Python implementation.
#' @seealso \code{\link{run.prepare.chronoclust}} for how to embed Python session within R session.
#' @import data.table
#' @import reticulate
#' @import Spectre
#' 
#' @param dat NO DEFAULT. Data.frame. Data to be clustered.
#' @param timepoint.col NO DEFAULT. Column name which represents the time point of each cell (data point) in dat.
#' @param use.cols NO DEFAULT. Vector of column names to use for clustering.
#' @param config DEFAULT = NULL. A named numeric list. List containing the value of ChronoClust's hyper-parameter. The name of the element must correspond to one of ChronoClust's parameter name such as epsilon, upsilon, etc. The numeric value must correspond to the value assigned for the corresponding parameter. 
#' \emph{Only include parameters that you want to override.} Those you prefer to set to default value need not be included in the list.
#' @param clust.name DEFAULT = "ChronoClust_cluster". Character. Name of the resulting 'cluster'
#'
#'
#'@usage
#'run.chronoclust(dat, timepoint.col, use.cols, config=NULL, clust.name = "ChronoClust_cluster")
#' 
#' 
#' @examples
#' # Read data
#' data.list <- Spectre::read.files(file.loc = PrimaryDirectory, file.type = ".csv", do.embed.file.names = TRUE)
#' cell.dat <- Spectre::do.merge.files(dat = data.list)
#' 
#' # Specify clustering column
#' ColumnNames <- as.matrix(unname(colnames(cell.dat)))
#' ClusteringColNos <- c(1:3)
#' ClusteringCols <- ColumnNames[ClusteringColNos]
#' 
#' # Specify time point column
#' timepoint.col <- 'day'
#' 
#' # Prepare Python session
#' run.prepare.chronoclust(environment_name = "chronoclust-R", 
#'                         create_environment = TRUE, 
#'                         install_dependencies = TRUE)
#' 
#' 
#' # Specify the parameters that need to be overriden
#' config <- list(beta = 0.2, 
#'                delta = 0.05, 
#'                epsilon = 0.03, 
#'                mu = 0.01)
#' 
#' cell.dat <- Spectre::run.chronoclust(dat=cell.dat, 
#'                                      timepoint.col=timepoint.col,
#'                                      use.cols=ClusteringCols,
#'                                      config=config)
#'
#' @author Givanna Putri, \email{ghar1821@@uni.sydney.edu.au}
#' @export

run.chronoclust <- function(dat,
                            timepoint.col,
                            use.cols,
                            config=NULL,
                            clust.name = "ChronoClust_cluster") {

  if(!is.element('reticulate', installed.packages()[,1])) stop('reticulate is required but not installed')
  if(!is.element('Spectre', installed.packages()[,1])) stop('Spectre is required but not installed')
  if(!is.element('data.table', installed.packages()[,1])) stop('data.table is required but not installed')
  
  require(reticulate)
  require(Spectre)
  require(data.table)
  
  ## Backup the data first
  dat.bk <- data.table::data.table(dat)
  
  ## Store this and remember to change working directory after the function finishes
  current.work.dir <- getwd()
  
  ## Create directory to store input file for ChronoClust. It requires each time point to be stored in separate csv file.
  ## Each csv file must also contain markers to be used for clustering.
  input.cc.dir <- paste(current.work.dir, 'input_chronoclust', sep = '/')
  dir.create(input.cc.dir, showWarnings = FALSE)
  setwd(input.cc.dir)
  
  ## Split the data into different csv files based on time point and columns for clustering.
  timepoints <- unique(dat.bk[[timepoint.col]])
  for (time.point in timepoints) {
    dat.subset <- dat.bk[dat.bk[[timepoint.col]] == time.point,]
    dat.subset <- dat.subset[, use.cols, with = FALSE]
    
    Spectre::write.files(dat.subset, time.point)
  }
  
  ## Store the files in a vector
  input.cc.files <- list.files(input.cc.dir, ".csv")
  input.cc.files <- paste(input.cc.dir, input.cc.files, sep="/")
  
  ## Setup Chronoclust
  # Create output directory
  output.cc.dir <- paste(current.work.dir, 'output_chronoclust', sep = '/')
  dir.create(output.cc.dir, showWarnings = FALSE)
  chronoclust <- import("chronoclust")

  # Run ChronoClust
  if (is.null(config)) {
    chronoclust$app$run(data=input.cc.files, 
                        output_directory=output.cc.dir)
  } else {
    all.chronoclust.param <- list(beta= 0.2,
                                  delta= 0.05,
                                  epsilon= 0.03,
                                  lambda= 2,
                                  k= 15,
                                  mu= 0.005,
                                  pi= 0,
                                  omicron= 0.00001,
                                  upsilon= 2)
    
    # work out which parameter has been specified by user
    user.def.config <- names(config)
    
    # work out which parameter has not been specified by user
    missing.chronoclust.param <- setdiff(names(all.chronoclust.param), user.def.config)
    
    # copy the config and modify this
    config.copy <- config
    for (missing.param in missing.chronoclust.param) {
      config.copy[[missing.param]] <- all.chronoclust.param[[missing.param]]
    }
    
    # run chronoclust
    chronoclust$app$run(data=input.cc.files, 
                        output_directory=output.cc.dir, 
                        param_beta=config.copy[['beta']], 
                        param_delta=config.copy[['delta']],
                        param_epsilon=config.copy[['epsilon']], 
                        param_lambda=config.copy[['lambda']], 
                        param_k=config.copy[['k']], 
                        param_mu=config.copy[['mu']],
                        param_pi=config.copy[['pi']], 
                        param_omicron=config.copy[['omicron']], 
                        param_upsilon=config.copy[['upsilon']])
  }
  
  # Read the result files and merge that into the dat as cluster
  setwd(output.cc.dir)
  timepoints.from.0 <-c(0: (length(timepoints)-1))
  clusters <- lapply(timepoints.from.0, function(tp) {
    cluster.dat <- read.csv(paste0("cluster_points_D", tp, ".csv"))
    return(cluster.dat$cluster_id)
  })
  names(clusters) <- timepoints
  
  # Prepare to append as column
  # First, convert the list into vector
  cluster.col <- unlist(clusters, use.names=FALSE)
  dat.bk[,clust.name] <- cluster.col

  ## Clean up
  # Delete the input file directory and the output directory
  unlink(input.cc.dir, recursive = TRUE)
  unlink(output.cc.dir, recursive = TRUE)
  # Set the working directory back to where it was
  setwd(current.work.dir)
  
  return(dat.bk)
  
}
