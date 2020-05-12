#' run.chronoclust - ...
#' @usage run.chronoclust(x, ...)
#'
#' @param data.files NO DEFAULT. A vector containing the location of csv file containing the time series datasets. 
#' Make sure that each 1 csv file is data for only 1 time point. 
#' Thus if you have 3 time points, this vector have 3 elements in it.
#' @param output.dir NO DEFAULT. The directory for Chronoclust to write its result to.
#' @param config DEFAULT NULL. Specify the parameter for ChronoClust as a list. If none is given, ChronoClust will
#' run on a default config. Only specify the parameter that you do not wish to use the default value for.
#'
#'
#' This function run a time-series clustering and cluster tracking algorithm Chronoclust.
#' ChronoClust was originally written in Python 3. It'll be run using Reticulate.
#' Please run run.prepare.chronoclust before running this function to
#' properly setup the python environment for running chronoclust.
#'
#' @export

run.chronoclust <- function(data.files,
                            output.dir,
                            config=NULL) {

  # Setup Chronoclust
  require(reticulate)
  chronoclust <- import("chronoclust")

  if (is.null(config)) {
    chronoclust$app$run(data=data.files, 
                        output_directory=output.dir)
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
    chronoclust$app$run(data=data.files, 
                        output_directory=output.dir, 
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
  
}
