#' do.subsample - ...
#'
#' @usage do.subsample(x, method, samp.col, targets, seed)
#'
#' @param x. Input dataframe with cells (rows) vs markers (columns). No default.
#' @param method. Character. Can be 'random' (downsampling single number from merged data), 'per.sample' (specifying number of cells each sample contributes) or 'min.per.sample' (each sample contributes the same amount of data based on sample with lowest count). No default
#' @param samp.col Character. Name of the column that reflects sample names.
#' @param targets List of downsample targets (if random, then absolute target; if per.sample, then targets per sample, must be in the same order as the unique sample names appear in the data frame or data table; can leave blank if min.per.sample). No default.
#' @param seed Numeric. Seed for reproducibility. Default = 42.
#'
#' @return Returns a subsampled data.table
#'
#' This function facilitates downsampling of a dataframe.
#'
#' @export

do.subsample <- function(x,
                         method, # random, per.sample
                         samp.col, # column than determines sample names
                         targets, # c(1000, 1500, ...)
                         seed = 42){

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
    #assign("subsample.res", subsample.res, envir = globalenv())
    return(subsample.res)
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
    #assign("subsample.res", subsample.res, envir = globalenv())
    return(subsample.res)

  }

  ## IF min.per.sample
  if(method == "min.per.sample") {
    # Create list of unique sample names
    sample.list <- unique(x[samp.col])
    sample.list <-sample.list[,1]
    sample.list

    #nrow.check = list()
    #for(i in c(1:(length(DataList)))){nrow.check[[i]] <- nrow(DataList[[i]])}
    #DownSampleTargets <- c(rep(nrow.check[[which.min(nrow.check)]], each=length(unique(AllSampleNos))))

    #min(data.frame(table(x[[samp.col]]))$Freq) #calculates count of each parameter (samp.col) in data (x), selecting the minimum number
    # Sets downsample target to be the same for each sample, based on whichever has the smallest number of cells
    targets <- c(rep(min(data.frame(table(x[[samp.col]]))$Freq), each=length(unique(x[[samp.col]]))))

    # Create res data.frame
    subsample.res <- data.frame()

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
    #assign("subsample.res", subsample.res, envir = globalenv())
    return(subsample.res)
  }

}
