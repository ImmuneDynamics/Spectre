#' run.chronoclust - ...
#' @usage run.chronoclust(x, ...)
#' 
#' @param 
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