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
#' 
#' @param data.files NO DEFAULT. A character vector. Path and name of the csv files storing the time series dataset. \strong{1 csv file per time point}.
#' @param output.dir NO DEFAULT. Character. Path to the directory where the clustering results will be stored.
#' @param config A named numeric list. List containing the value of ChronoClust's hyper-parameter. The name of the element must correspond to one of ChronoClust's parameter name such as epsilon, upsilon, etc. The numeric value must correspond to the value assigned for the corresponding parameter. 
#' \emph{Only include parameters that you want to override.} Those you prefer to set to default value need not be included in the list.
#'
#'
#'@usage
#'run.chronoclust(data.files, output.dir, config = NULL)
#' 
#' 
#' @examples
#' # Prepare Python session
#' run.prepare.chronoclust(environment_name = "chronoclust-R", 
#'                         create_environment = TRUE, 
#'                         install_dependencies = TRUE)
#' 
#' # Specify the data files location and output directory
#' data.files <- c("cc_example/day1.csv", "cc_example/day2.csv")
#' output.dir <- "cc_example/output"
#' 
#' # Specify the parameters that need to be overriden
#' config <- list(beta = 0.2, 
#'                delta = 0.05, 
#'                epsilon = 0.03, 
#'                mu = 0.01)
#' 
#' run.chronoclust(data.files = data.files, 
#'                 output.dir = output.dir, 
#'                 config = config)
#'
#' @author Givanna Putri, \email{ghar1821@@uni.sydney.edu.au}
#' @export

run.chronoclust <- function(dat,
                            timepoint.col,
                            use.cols,
                            config=NULL) {

  if(!is.element('reticulate', installed.packages()[,1])) stop('reticulate is required but not installed')
  if(!is.element('Spectre', installed.packages()[,1])) stop('Spectre is required but not installed')
  
  require(reticulate)
  require(Spectre)
  
  ## Store this and remember to change working directory after the function finishes
  current.work.dir <- getwd()
  
  ## Create directory to store input file for ChronoClust. It requires each time point to be stored in separate csv file.
  ## Each csv file must also contain markers to be used for clustering.
  input.cc.dir <- 'input_chronoclust'
  dir.create(input.cc.dir, showWarnings = FALSE)
  setwd(input.cc.dir)
  
  ## Split the data into different csv files based on time point and columns for clustering.
  Spectre::do.split.data(dat, timepoint.col, use.cols)
  
  ## Store the files in a vector
  input.cc.files <- list.files(input.cc.dir, ".csv")
  input.cc.files <- paste(input.cc.dir, input.cc.files, sep="/")
  
  
  ## Setup Chronoclust
  # Create output directory
  output.cc.dir <- 'output_chronoclust'
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
  
  # Delete the input file directory and the output directory
  unlink(input.cc.dir, recursive = TRUE)
  unlink(output.cc.dir, recursive = TRUE)
  # Set the working directory back to where it was
  setwd(current.work.dir)
  
}
