#' Statistical analysis
#'
#' @import data.table
#' 
#' @usage NULL
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
  
 warning("'do.stats' function has been depreciated. Please use 'create.stats' instead")
  
}