#' run.chronoclust - ...
#' @usage run.chronoclust(x, ...)
#'
#' @param data.files NO DEFAULT. A vector containing the location of csv file containing the time series datasets. 
#' Make sure that each 1 csv file is data for only 1 time point. 
#' Thus if you have 3 time points, this vector have 3 elements in it.
#' @param output.dir NO DEFAULT. The directory for Chronoclust to write its result to.
#' @param conda.env.name NO DEFAULT. The conda virtual environment name where Chronoclust wil run.
#' @param conda.env.location DEFAULT NULL. Specify where anaconda is installed. If nothing is given, 
#' then reticulate will look at where it's normally installed /Users/<yourname>/anaconda.
#' @param config DEFAULT NULL. Specify the parameter for ChronoClust as a list. If none is given, ChronoClust will
#' run on a default config. Only specify the parameter that you do not wish to use the default value for.
#'
#'
#' This function run a time-series clustering and cluster tracking algorithm Chronoclust.
#' ChronoClust was originally written in Python 3. It'll be run using Reticulate.
#' Only anaconda/miniconda virtual environment is supported by this function. This environment need to be setup before running the function.
#' Please visit the various guide in the following link to install anaconda/miniconda and setup an environment:
#' https://docs.conda.io/projects/conda/en/latest/user-guide/index.html
#' ChronoClust requires the following python package to run: numpy, pandas, scikit-learn, tqdm, numba.
#' 
#' ChronoClust package is available for download through Github: https://github.com/ghar1821/chronoclust.
#' To install ChronoClust in your environment, please visit the Github page and follow the instruction:
#' "How do I use Chronoclust".
#' 
#' Please make sure that your environment is running Python 3.7.
#' ChronoClust is not supported for Python 2.
#'
#' @export

run.chronoclust <- function(data.files,
                            output.dir,
                            conda.env.name,
                            conda.env.location=NULL,
                            config=NULL) {
  # activate the conda env
  if (is.null(conda.env.location)) {
    use_condaenv(conda.env.name)
  } else {
    use_condaenv(conda.env.name, conda=conda.env.location)
  }


  # Setup Chronoclust
  chronoclust <- import("chronoclust")  # Import Chronoclust API

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
