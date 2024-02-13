#' Combine two columns
#'
#' Combine the values of two columns into a single column, with an option to remove any 'NA' values.
#'
#' @usage do.combine.cols(dat, col1, col2, na.rm)
#'
#' @param dat NO DEFAULT. Data.table.
#' @param col1 NO DEFAULT. Name of the column you want to combine into.
#' @param col2 NO DEFAULT. Name of the column you want to combine from. This column will be deleted after combination into col1.
#' @param na.rm DEFAULT = TRUE. If NA values are present in the columns, they will be removed.
#'
#' @return Returns a data.table with adjusted columns
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @export

do.combine.cols <- function(dat,
                            col1,
                            col2,
                            na.rm = TRUE) {
  # require: data.table
  

  ### NA checks

  na.checks <- rowSums(is.na(dat[, c(col1, col2), with = FALSE]))

  if (any(na.checks < 1)) {
    stop("Less that one 'NA' per row")
  }

  ### Check class of col1

  cls <- class(dat[[col1]])

  ### Column combination

  if (isTRUE(na.rm)) {
    dat[[col1]][is.na(dat[[col1]])] <- "PLACEHOLDER"
    dat[[col2]][is.na(dat[[col2]])] <- "PLACEHOLDER"
    dat[[col1]] <- paste0(dat[[col1]], dat[[col2]])
    dat[[col1]] <- gsub("PLACEHOLDER", "", dat[[col1]])
    dat[[col2]] <- NULL
  } else {
    dat[[col1]] <- paste0(dat[[col1]], dat[[col2]])
    dat[[col2]] <- NULL
  }

  ### Class adjust

  class(dat[[col1]]) <- cls

  ### Return

  return(dat)
}
