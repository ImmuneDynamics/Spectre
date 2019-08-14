#' annotate.files - Add column of new values or keywords to samples within a list of sample dataframes
#'
#' This function allows you to add a column of new values or keywords to the rows (cells) of each sample dataframe within a list.
#'
#' @param x The list of sample dataframes. No default.
#' @param col.name Name of the new column. No Default.
#' @param d List of values or keywords to add to each. No default.
#' @param assignments A list designating which new values or keywords are added to each sample. No default.
#' @param add.num Logical, add a 'numbered' version of the new values or keywords. Defaults to FALSE.
#'
#' @return a list (called data.list) of dataframes -- one file per dataframe -- with the a column of the new values or keywords added.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references \url{https://sydneycytometry.org.au/spectre}
#'
#' @usage See \url{https://sydneycytometry.org.au/spectre} for usage instructions and vignettes.
#'
#' @examples
#' annotate.files()
#'
#' @export

annotate.files <- function(x,
                          col.name,
                          d,
                          assignments,
                          add.num = FALSE)
  {

  ## For testing
      #x = data.list
      #col.name = "Sample"
      #d = sample.names
      #assignments = sample.assign
      #add.num = TRUE

  ## Add keywords to dataset
      num.of.d <- c(1:length(d))
      for(a in num.of.d){
        for(i in c(assignments[[a]])){
          x[[i]][[col.name]] <- d[[a]]

          if(add.num == TRUE){
            n <- paste0(col.name, "_Num")
            x[[i]][[n]] <- a
          }

        }
      }

  assign("data.list", x, envir = globalenv())
}

