#' run.chronoclust - ...
#' @usage run.chronoclust(x, ...)
#'
#' @param config.xml.file.location NO DEFAULT. Location of xml file describing the parameter values for Chronoclust to operate on.
#' @param input.xml.file.location NO DEFAULT. Location of xml file describing the location of the dataset (for each timepoint) for Chronoclust to operate on.
#' @param output.dir.location NO DEFAULT. The directory for Chronoclust to write its result to.
#' @param conda.env.name NO DEFAULT. The conda virtual environment name where Chronoclust wil run.
#' @param conda.env.location DEFAULT NULL. Specify where anaconda is installed. If nothing is given, then reticulate will look at where it's normally installed /Users/<yourname>/anaconda
#'
#'
#' This function run a time-series clustering and cluster tracking algorithm Chronoclust.
#' ChronoClust was originally written in Python 3. It'll be run using Reticulate.
#' Only anaconda/miniconda virtual environment is supported by this function. This environment need to be setup before running the function.
#' Please visit the various guide in the following link to install anaconda/miniconda and setup an environment:
#' https://docs.conda.io/projects/conda/en/latest/user-guide/index.html
#' ChronoClust requires the following python package to run: numpy, pandas, scikit-learn, tqdm, numba.
#' You can either install them separately after setting up your anaconda environment
#' or just simply import them from the prepared environment script provided in our github repository (worked_example/Chronoclust_workflow)
#'
#' Please make sure that your environment is running Python 3 (preferably 3.7).
#' ChronoClust is not supported for Python 2.
#'
#' @export

run.chronoclust <- function(config.xml.file.location,
                            input.xml.file.location,
                            output.dir.location,
                            conda.env.name,
                            conda.env.location=NULL) {
  # activate the conda env
  if (is.null(conda.env.location)) {
    use_condaenv(conda.env.name)
  } else {
    use_condaenv(conda.env.name, conda=conda.env.location)
  }


  # Setup Chronoclust
  chronoclust <- import("chronoclust")  # Import Chronoclust API

  chronoclust$main$main$run(config_xml = config.xml.file.location,
                              input_xml = input.xml.file.location,
                              log_dir = output.dir.location,
                              output_dir = output.dir.location)


}
