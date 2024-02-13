#' do.reorder - custom row re-ordering of a data.table
#'
#' This function allows for a user-defined re-ordering of data.table rows
#'
#' @param dat NO DEFAULT. A data.table
#' @param use.col DEFAULT = NULL. The column to use for re-ordering
#' @param new.order DEFAULT = NULL. A vector of row values for the specified column, in the desired order
#'
#' @usage do.reorder(dat, use.col, new.order)
#'
#' @author
#' Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#' Givanna Putri
#'
#'
#' @export

do.reorder <- function(dat, use.col, new.order) {
    ### Packages

    # require: data.table

    ### Test

    # dat <- Spectre::demo.clustered use.col <- 'FileName' as.matrix(unique(dat[[use.col]])) new.order <-
    # unique(dat[[use.col]])[c(3,4,1,2,5:12)]

    # New implementation that might be quicker?  check the values in new.order is in the column to reorder dat

    if (!setequal(new.order, unique(dat[[use.col]])))
        stop("The values provided in 'new.order' do not reflect all unique entries in 'use.col' column")


    dat[[use.col]] <- factor(dat[[use.col]], levels = new.order)
    data.table::setorderv(dat, use.col)
    
    return(dat)



    ### Checks

    # if (!all(sort(new.order) == sort(unique(dat[[use.col]])))) { 
    #   stop("The values provided in 'new.order' do not reflect all unique entries in dat[[use.col]]") 
    # } 
    # ### Add new column new.order.
    # nums <- c(1:length(new.order)) 
    # new.order <- data.table(new.order, new.order.nums) 
    # dat <- do.add.cols(dat, use.col, new.order, 'new.order', show.status = FALSE)
    # setorderv(dat, 'new.order.nums') ### Cleanup dat[['new.order.nums']] <- NULL 
    # ### Return 
    # return(dat)
}

