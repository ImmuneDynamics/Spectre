#' do.stats
#' 
#' The function has been depreciated. Please use 'create.stats' instead
#' 
#' @usage do.stats(dat, use.cols, sample.col, grp.col, comparisons)
#'
#' @import data.table
#'
#' @export do.stats

do.stats <- function(dat,
                     use.cols,
                     sample.col,
                     grp.col,
                     
                     comparisons, # make a tool to create a factorial comparison design -- for now just specify manually
                     
                     variance.test = "kruskal.test", ## add ANOVA
                     pairwise.test = "wilcox.text", ## Add t-test
                     corrections = "fdr"){
    .Deprecated("create.stats")
    # stop("The 'do.stats' function has been depreciated. Please use 'create.stats' instead")
  
}