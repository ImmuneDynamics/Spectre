#' List switch or inversion
#'
#' @usage do.list.switch(dat)
#'
#' @param dat NO DEFAULT. A list of vectors that you wish to switch/invert
#'
#' @return Returns a data.table of the switched/inverted li
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}.
#'
#' @examples
#' lst <- list("A" = c(1,2,3), "B" = c(4,5,6), "C" = c(7,8,9))
#' lst
#'
#' res <- do.list.switch(lst)
#' res
#'
#' @import data.table
#'
#' @export

do.list.switch <- function(dat){

  res.list <- list()

  for(i in c(1:length(dat))){
    # i <- 1
    a <- names(dat)[[i]]

    temp <- dat[[a]]
    temp.dt <- data.table('Values' = temp,
                          'New' = rep(a, length(temp)))

    res.list[[i]] <- temp.dt
  }
  res.dat <- rbindlist(res.list)
  return(res.dat)
}
