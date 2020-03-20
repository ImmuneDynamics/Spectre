#' do.subsample - ...
#'
#' @usage do.subsample(dat, method, samp.col, targets, seed)
#'
#' @param dat NO DEFAULT. Input dataframe with cells (rows) vs markers (columns).
#' @param method NO DEFAULT. Character. Can be 'random' (downsampling single number from merged data), 'per.sample' (specifying number of cells each sample contributes) or 'min.per.sample' (each sample contributes the same amount of data based on sample with lowest count).
#' @param samp.col NO DEFAULT. Character. Name of the column that reflects sample names.
#' @param targets NO DEFAULT. List of downsample targets (if random, then absolute target; if per.sample, then targets per sample, must be in the same order as the unique sample names appear in the data frame or data table; can leave blank if min.per.sample).
#' @param seed DEFAULT = 42. Numeric. Seed for reproducibility.
#'
#' @return Returns a subsampled data.table
#'
#' This function facilitates downsampling of a dataframe.
#'
#' @export

do.subsample <- function(dat,
                         method, # random, per.sample
                         samp.col, # column than determines sample names
                         targets, # c(1000, 1500, ...)
                         seed = 42){

  ## Test data
      #dat <- cell.dat
      #method <- "per.sample"
      #samp.col <- "FileName"
      #targets <- c(rep(1000, 12))
      #seed <- 42

  ## IF random # WORKS
  if(method == "random"){
    set.seed(seed)
    subsample.res <- dat[sample(1:nrow(dat), targets), ]
    subsample.res <- as.data.frame(subsample.res)
    #assign("subsample.res", subsample.res, envir = globalenv())
    return(subsample.res)
  }

  ## IF per.sample
  if(method == "per.sample"){

    # Create list of unique sample names
    sample.list <- unique(dat[samp.col])
    sample.list <-sample.list[,1]
    sample.list

    # Create res data.frame
    subsample.res <- data.frame()

    # Loop
    for (i in c(1:length(sample.list))) {
      nam <- sample.list[i]
      nsub <- targets[i]
      data.temp <- subset(dat, dat[[samp.col]] == nam) # works
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
    sample.list <- unique(dat[samp.col])
    sample.list <-sample.list[,1]
    sample.list

    #nrow.check = list()
    #for(i in c(1:(length(DataList)))){nrow.check[[i]] <- nrow(DataList[[i]])}
    #DownSampleTargets <- c(rep(nrow.check[[which.min(nrow.check)]], each=length(unique(AllSampleNos))))

    #min(data.frame(table(dat[[samp.col]]))$Freq) #calculates count of each parameter (samp.col) in data (dat), selecting the minimum number
    # Sets downsample target to be the same for each sample, based on whichever has the smallest number of cells
    targets <- c(rep(min(data.frame(table(dat[[samp.col]]))$Freq), each=length(unique(dat[[samp.col]]))))

    # Create res data.frame
    subsample.res <- data.frame()

    for (i in c(1:length(sample.list))) {
      nam <- sample.list[i]
      nsub <- targets[i]
      data.temp <- subset(dat, dat[[samp.col]] == nam) # works
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
