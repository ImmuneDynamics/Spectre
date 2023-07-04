#' spectre object
#'
#' @export

setClass("spectre",representation(data='list',
                                  analysis="list",
                                  other='list',
                                  cellid='vector'))
