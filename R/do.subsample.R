#' Subsample data
#'
#' Method to subsample data.
#' Can subsample by randomly selecting a desired number of cells from all samples (method = "random"), subsample by specifying the exact number of cells for each sample (method= "per.sample"), or by subsampling the same number of cells from each sample based on the sample with the lowest count (method = "min.per.sample").
#' Useful to decrease total cells for generating dimensionality reduction plots (tSNE/UMAP).
#'
#' @param dat NO DEFAULT. Input dataframe with cells (rows) vs markers (columns).
#' @param targets NO DEFAULT. List of downsample targets. If divide.by is specified, then must be a vector of subsample targets in the same order as the unique divide.by entries.
#' @param divide.by DEFAULT = NULL. Character. Name of the column that reflects groupings of cells (sample names, group names etc) if you want to subsample by each.
#' @param min.per DEFAULT = FALSE. If TRUE, and samp.col is specified, each sample contributes the same amount of data based on sample with lowest count.
#' @param seed DEFAULT = 42. Numeric. Seed for reproducibility.
#'
#' @usage do.subsample(dat, targets, samp.col, min.per, seed)
#'
#' @examples
#' Spectre::do.subsample(dat = Spectre::demo.start,
#'                       targets = 10000)
#'
#' @author
#' Thomas Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#' Felix Marsh-Wakefield, \email{felix.marsh-wakefield@@sydney.edu.au}
#'
#' @export

do.subsample <- function(dat,
                         targets,

                         divide.by = NULL,
                         min.per = FALSE,
                         seed = 42){

  ## Check that necessary packages are installed
  if(!is.element('Spectre', installed.packages()[,1])) stop('Spectre is required but not installed')
  if(!is.element('data.table', installed.packages()[,1])) stop('data.table is required but not installed')

  ## Require packages
  require(Spectre)
  require(data.table)

  ## IF random
  if(is.null(divide.by)){
    set.seed(seed)
    subsample.res <- dat[sample(1:nrow(dat), targets), ]
    subsample.res <- as.data.table(subsample.res)
    return(subsample.res)
  }

  ## IF per.sample
  if(!is.null(divide.by)){

    if(min.per == FALSE){
      # Create list of unique sample names
      sample.list <- unique(dat[,divide.by,with = FALSE])
      sample.list <- sample.list[[1]]

      # Create res data.frame
      subsample.res <- data.frame()

      # Loop
      for (i in c(1:length(sample.list))) {
        nam <- sample.list[i]
        nsub <- targets[i]

        data.temp <- dat[dat[[divide.by]] == nam,]

        nrow(data.temp)
        set.seed(seed)
        data.temp <- data.temp[sample(1:nrow(data.temp), nsub), ]
        nrow(data.temp)
        subsample.res <- rbind(subsample.res, data.temp)
      }

      dim(subsample.res)
      #assign("subsample.res", subsample.res, envir = globalenv())
      subsample.res <- as.data.table(subsample.res)
      return(subsample.res)
    }

    ## IF min.per.sample
    if(min.per == TRUE){
      # Create list of unique sample names
      sample.list <- unique(dat[divide.by])
      sample.list <- sample.list[,1]
      sample.list

      #nrow.check = list()
      #for(i in c(1:(length(DataList)))){nrow.check[[i]] <- nrow(DataList[[i]])}
      #DownSampleTargets <- c(rep(nrow.check[[which.min(nrow.check)]], each=length(unique(AllSampleNos))))

      #min(data.frame(table(dat[[divide.by]]))$Freq) #calculates count of each parameter (divide.by) in data (dat), selecting the minimum number
      # Sets downsample target to be the same for each sample, based on whichever has the smallest number of cells
      targets <- c(rep(min(data.frame(table(dat[[divide.by]]))$Freq), each=length(unique(dat[[divide.by]]))))

      # Create res data.frame
      subsample.res <- data.frame()

      for (i in c(1:length(sample.list))) {
        nam <- sample.list[i]
        nsub <- targets[i]
        data.temp <- subset(dat, dat[[divide.by]] == nam) # works
        nrow(data.temp)
        set.seed(seed)
        data.temp <- data.temp[sample(1:nrow(data.temp), nsub), ]
        nrow(data.temp)
        subsample.res <- rbind(subsample.res, data.temp)
      }
      dim(subsample.res)
      #assign("subsample.res", subsample.res, envir = globalenv())
      subsample.res <- as.data.table(subsample.res)
      return(subsample.res)
    }
  }

}
