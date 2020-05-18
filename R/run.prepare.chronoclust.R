#' run.prepare.chronoclust - This function will prepare a Python environment to run ChronoClust.
#' You need to install anaconda package management software for this to work.
#' If Reticulate is kind, it might prompt you to install it when running this function if you haven't got it installed.
#' However, please do not count on it.
#' If it does not prompt you, and you encounter an error such as no anaconda installation is found,
#' please visit https://www.anaconda.com/products/individual to manually download and install anaconda.
#' 
#' For advanced users, if your anaconda is installed in custom location,
#' you can specify this location as environment_path parameter.
#' @usage run.prepare.chronoclust(environment_name, environment_path, create_environment, install_dependencies)
#'
#' @param environment_name NO DEFAULT. Character. The name of the anaconda environment where chronoclust will execute.
#' @param environment_path DEFAULT NULL. Character. For expert users only. If you have custom
#' anaconda installation i.e. not installed to default location, you can specify
#' the location here.
#' @param create_environment DEFAULT TRUE. Logical. Whether to newly create the environment or not. If FALSE, make sure the
#' environment exists in your anaconda installation.
#' @param install_dependencies DEFAULT TRUE. Logical. Whether to install chronoclust and all its dependencies into the environment.
#' If TRUE, it will install all the dependencies, and restart R session.
#' When R session is restarted, the variables you have previously assigned
#' will remain intact, but you will need to reload all the libraries again.
#' 
#' @author Givanna Putri, \email{ghar1821@@uni.sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}.
#'
#' @examples
#' run.prepare.chronoclust(environment_name = "chronoclust-R", create_environment = TRUE, install_dependencies = TRUE)
#'
#' @export

run.prepare.chronoclust <- function(environment_name,
                                            environment_path=NULL,
                                            create_environment=TRUE,
                                            install_dependencies=TRUE) {
  if(!is.element('reticulate', installed.packages()[,1])) stop('reticulate is required but not installed')
  require(reticulate)
  
  if(create_environment) {
    conda_create(environment_name)
  } else {
    if (is.null(environment_path)) {
      use_condaenv(environment_name, required = TRUE)
    } else {
      use_condaenv(environment_name, conda=environment_path, required = TRUE)
    }
    
  }
  
  if (install_dependencies) {
    # install chronoclust dependencies
    conda_install(envname = environment_name, 
                  python_version = "3.7.0",
                  packages = c("numpy", 
                               "pandas", 
                               "scipy", 
                               "scikit-learn", 
                               "tqdm", 
                               "numba"))
    conda_install(envname = environment_name,
                  packages = c("Chronoclust"),
                  pip=TRUE)
    
    # install chronoclust
    #download.file(url = "https://github.com/ghar1821/Chronoclust/archive/master.zip",
    #              destfile = "chronoclust-master.zip")
    #unzip(zipfile = "chronoclust-master.zip")
    #WorkingDirectory <- getwd()
    #setwd(paste(WorkingDirectory, "chronoclust-master", sep = '/'))
    #py_run_string(paste0("from setuptools import sandbox; sandbox.run_setup('",getwd(),"/setup.py', ['install'])"))
    #chronoclust <- import_from_path("chronoclust")
    #setwd(WorkingDirectory)
    .rs.restartR()
  }
}
