#' List switch or inversion
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
