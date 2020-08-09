#' Summarise the column names and row numbers for elements of a list
#'
#' @import data.table
#'
#' @export

do.list.summary <- function(dat)

{

  res.list <- list()

  ncol.check <- list()
  colName.check <- list()
  nrow.check <- list()

  for(i in c(1:(length(dat)))){ncol.check[[i]] <- length(names(dat[[i]]))} # creates a list of the number of columns in each sample
  for(i in c(1:(length(dat)))){colName.check[[i]] <- names(dat[[i]])}
  name.table <- data.frame(matrix(unlist(colName.check), nrow = length(dat), byrow = T))
  for(i in c(1:(length(dat)))){nrow.check[[i]] <- nrow(dat[[i]])}

  ncol.check <- as.matrix(ncol.check)
  nrow.check <- as.matrix(nrow.check)

  res.list$name.table <- name.table
  res.list$ncol.check <- ncol.check
  res.list$nrow.check <- nrow.check

  return(res.list)
}

