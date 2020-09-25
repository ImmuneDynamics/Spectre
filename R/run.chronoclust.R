#' Run the ChronoClust algorithm
#'
#' Method to run ChronoClust clustering algorithm.
#' For this to work properly, a Python session must already be embedded within the currently running R session.
#' This is because ChronoClust is originally written in Python 3.
#' Spectre includes a function \code{\link{run.prepare.chronoclust}} that uses Reticulate to prepare and embed Python session within R session.
#'
#' @seealso \url{https://github.com/ghar1821/Chronoclust} for ChronoClust's Python implementation.
#' @seealso \code{\link{run.prepare.chronoclust}} for how to embed Python session within R session.
#'
#' @param dat NO DEFAULT. Data.frame. Data to be clustered.
#' @param timepoint.col NO DEFAULT. Column name which represents the time point of each cell (data point) in dat.
#' @param timepoints NO DEFAULT. The time points (in order).
#' @param use.cols NO DEFAULT. Vector of column names to use for clustering.
#' @param config DEFAULT = NULL. A named numeric list. List containing the value of ChronoClust's hyper-parameter. The name of the element must correspond to one of ChronoClust's parameter name such as epsilon, upsilon, etc. The numeric value must correspond to the value assigned for the corresponding parameter.
#' \emph{Only include parameters that you want to override.} Those you prefer to set to default value need not be included in the list.
#' @param clust.name DEFAULT = "ChronoClust_cluster". Character. Name of the resulting 'cluster'
#' @param clean.up DEFAULT = FALSE. Whether to remove the files chronoclust produces
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
#' cluster.cols <- names(cell.dat)[c(1:8,10:20)]
#'
#' # Specify time point column
#' timepoint.col <- "Group"
#' timepoints <- c("Mock", "WNV-01", "WNV-02", "WNV-03", "WNV-04", "WNV-05")
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
#'                                      timepoints=timepoints,
#'                                      use.cols=cluster.cols,
#'                                      config=config)
#'
#' @author Givanna Putri, \email{givanna.haryonoputri@@sydney.edu.au}
#' @export

run.chronoclust <- function(dat,
                            timepoint.col,
                            timepoints,
                            use.cols,
                            config = NULL,
                            clust.name = "ChronoClust_cluster",
                            clean.up = FALSE) {
  
  if(!is.element('reticulate', installed.packages()[,1])) stop('reticulate is required but not installed')
  if(!is.element('Spectre', installed.packages()[,1])) stop('Spectre is required but not installed')
  if(!is.element('data.table', installed.packages()[,1])) stop('data.table is required but not installed')
  
  require(reticulate)
  require(Spectre)
  require(data.table)
  
  # testing
  # timepoint.col = 'Group'
  # timepoints = c('Mock', 'WNV-01', 'WNV-02', 'WNV-03', 'WNV-04', 'WNV-05')
  # use.cols = c(1:8,10:20)
  
  message("Preparing data")
  
  ## Backup the data first
  dat.bk <- data.table(dat)
  
  ## Order the data based on the timepoints given
  dat.bk <- dat.bk[order(match(dat.bk[[timepoint.col]], as.character(timepoints)))]
  
  
  ## Store this and remember to change working directory after the function finishes
  current.work.dir <- getwd()
  
  ## Create directory to store input file for ChronoClust. It requires each time point to be stored in separate csv file.
  ## Each csv file must also contain markers to be used for clustering.
  input.cc.dir <- paste(current.work.dir, 'input_chronoclust', sep = '/')
  dir.create(input.cc.dir, showWarnings = FALSE)
  setwd(input.cc.dir)
  
  ## Split the data into different csv files based on time point and columns for clustering.
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
  
  # Run ChronoClust
  # TODO update the python version to reflect this
  default.cc.params <- list(beta= 0.2,
                            delta= 0.05,
                            epsilon= 0.03,
                            lambda= 2,
                            k= 15,
                            mu= 0.005,
                            pi= 0,
                            omicron= 0.00001,
                            upsilon= 2)
  if (is.null(config)) {
    config.copy <- default.cc.params
  } else {
    config.copy <- config
  }
  
  # work out which parameter has been specified by user
  user.def.config <- names(config.copy)
  
  # work out which parameter has not been specified by user
  missing.chronoclust.param <- setdiff(names(default.cc.params), user.def.config)
  
  # copy the config and modify this
  for (missing.param in missing.chronoclust.param) {
    config.copy[[missing.param]] <- default.cc.params[[missing.param]]
  }
  
  message("Running ChronoClust")
  # run chronoclust
  chronoclust <- import("chronoclust")
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
  
  
  setwd(output.cc.dir)
  message("Inferring lineage ID and association ID")
  
  cell_dat_list = list()
  result_dat <- fread("result.csv")
  timepoints.from.0 <- c(0: (length(timepoints)-1))
  for (i in timepoints.from.0) {
    time.point <- timepoints[i+1]
    dat.subset <- dat.bk[dat.bk[[timepoint.col]] == time.point,]
    
    # read cluster points
    cluster.dat <- fread(paste0("cluster_points_D", i, ".csv"))
    dat.subset[[paste0(clust.name, '_lineage')]] <- cluster.dat$cluster_id
    
    result_dat_sub <- result_dat[result_dat$timepoint == i,]
    
    res <- result_dat_sub$tracking_by_association
    names(res) <- result_dat_sub$tracking_by_lineage
    res['None'] <- 'None'
    
    assoc <- sapply(cluster.dat$cluster_id, function(cl) {
      res[cl]
    })
    
    dat.subset[[paste0(clust.name, '_association')]] <- as.vector(assoc)
    cell_dat_list[[time.point]] <- dat.subset
  }
  
  message("Stitching results back to data frame")
  cell_dat <- rbindlist(cell_dat_list)
  
  if (clean.up) {
    message("Cleaning up")
    ## Clean up
    # Delete the input file directory and the output directory
    unlink(input.cc.dir, recursive = TRUE)
    unlink(output.cc.dir, recursive = TRUE)
  }
  
  # Set the working directory back to where it was
  setwd(current.work.dir)
  
  return(cell_dat)
  
}
