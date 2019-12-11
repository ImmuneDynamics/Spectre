#' run.minmax.normaliser - ...
#' 
#' @usage run.minmax.normaliser(x, ...)
#' 
#' @param data NO DEFAULT. A dataframe containing cells (rows) vs features/markers (columns) to be normalised.
#' 
#' Normalise the data to range between 0 and 1.
#' 
#' @return Normalised data
#' 
#' @export
run.minmax.normaliser <- function(data) {
  norm.data <- (data - min(data))/(max(data)-min(data))
  return(norm.data)
}