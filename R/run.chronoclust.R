#' run.chronoclust - ...
#' @usage run.chronoclust(x, ...)
#' 
#' @param config.xml.file.location NO DEFAULT. Location of xml file describing the parameter values for Chronoclust to operate on.
#' @param input.xml.file.location NO DEFAULT. Location of xml file describing the location of the dataset (for each timepoint) for Chronoclust to operate on.
#' @param output.dir.location NO DEFAULT. The directory for Chronoclust to write its result to.
#' @param virtualenv.name DEFAULTS to r-chronoclust. The name of the virtual environment to run Chronoclust on.
#' @param create.virtualenv DEFAULTS to TRUE. Whether to create the virtual environment to run Chronoclust on.
#' @param delete.virtualenv DEFAULTS to TRUE. Whether to delete the virtual environment to run Chronoclust on after it finishes.
#' 
#' This function run a time-series clustering and cluster tracking algorithm Chronoclust.
#' Chronoclust is originally written in python. Hence it requires Reticulate to run.
#' Reticulate executes Chronoclust on a virtual environment. 
#' 
#' @export

run.chronoclust <- function(config.xml.file.location,
                            input.xml.file.location,
                            output.dir.location,
                            virtualenv.name = 'r-chronoclust',
                            create.virtualenv = TRUE,
                            delete.virtualenv = TRUE) {
  # Create new virtual env if needed
  if (create.virtualenv) {
    virtualenv_create(virtualenv.name)  
  }
  # install chronoclust
  virtualenv_install(virtualenv.name, "chronoclust")
  
  # Setup Chronoclust
  chronoclust <- import("chronoclust")  # Import Chronoclust API
  
  chronoclust$chronoclust$run(config_xml = config.xml.file.location, 
                              input_xml = input.xml.file.location, 
                              log_dir = output.dir.location, 
                              output_dir = output.dir.location)
  # Clean up if needed (may user want to reuse?)
  if (delete.virtualenv) {
    virtualenv_remove(virtualenv.name)  
  }
  
}