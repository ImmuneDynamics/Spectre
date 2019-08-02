# Spectre::subsample

subsample <- function(x,
                      method, # random, per.sample
                      samp.col, # column than determines sample names
                      targets, # c(1000, 1500, ...)
                      seed){ # 42

  ## Test data
  #x <- cell.dat
  #method <- "per.sample"
  #samp.col <- "FileName"
  #targets <- c(rep(1000, 12))
  #seed <- 42

  ## IF random # WORKS
  if(method == "random"){
    set.seed(seed)
    subsample.res <- x[sample(1:nrow(x), targets), ]
    subsample.res <- as.data.frame(subsample.res)
    assign("subsample.res", subsample.res, envir = globalenv())
  }

  ## IF per.sample
  if(method == "per.sample"){

    # Create list of unique sample names
    sample.list <- unique(x[samp.col])
    sample.list <-sample.list[,1]
    sample.list

    # Create res data.frame
    subsample.res <- data.frame()

    # Loop
    for (i in c(1:length(sample.list))) {
      nam <- sample.list[i]
      nsub <- targets[i]
      data.temp <- subset(x, x[[samp.col]] == nam) # works
      nrow(data.temp)
      set.seed(seed)
      data.temp <- data.temp[sample(1:nrow(data.temp), nsub), ]
      nrow(data.temp)
      subsample.res <- rbind(subsample.res, data.temp)
    }
    dim(subsample.res)
    assign("subsample.res", subsample.res, envir = globalenv())
  }
}
