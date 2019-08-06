# colour.plot.loop

colour.plot.loop <- function(d,
                             x.axis,
                             y.axis,
                             col.axis,
                             diff.by, # name of column to differentiate by
                             path, # path of output directory
                             title,
                             colours,
                             dot.size){

  ## input = dataframe, embedded names to differentiate by



  ## filter to numeric columns only (others use factor.plots)
      numeric.only <- sapply(d, is.numeric)
      d <- numeric.only

  ## set X and Y limits
      Xmax <- max(d[[x.axis]])
      Ymax <- max(d[[y.axis]])

      Xmin <- min(d[[x.axis]])
      Ymin <- min(d[[y.axis]])

  ##

      subset(x, x[[samp.col]] == nam)




}


