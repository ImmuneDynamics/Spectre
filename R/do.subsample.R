#' Subsample data
#'
#' Method to subsample data.
#' Can subsample by randomly selecting a desired number of cells from all samples (DEFAULT), subsample by specifying the exact number of cells for each sample (specify divide.by), or by subsampling the same number of cells from each sample based on the sample with the lowest count (specify divide.by and min.per).
#' Useful to decrease total cells for generating dimensionality reduction plots (tSNE/UMAP).
#'
#' @param dat NO DEFAULT. Input dataframe with cells (rows) vs markers (columns).
#' @param targets NO DEFAULT. Vector of downsample targets. If divide.by is specified, then must be a vector of subsample targets in the same order as the unique divide.by entries (e.g. unique(dat[[divide.by]])). Can also provide as a data.table or data.frame where the first column is the unique entries in the divide.by argument (i.e. unique(dat[[divide.by]])), and the second column should be the targets. In this case, does not have to be in the order they appear in the dataset, but the 'divide.by' argument must be set.
#' @param divide.by DEFAULT = NULL. Character. Name of the column that reflects groupings of cells (sample names, group names etc) if you want to subsample by each.
#' @param min.per DEFAULT = FALSE. If TRUE, and divide.by is specified, each sample contributes the same amount of data based on sample with lowest count.
#' @param seed DEFAULT = 42. Numeric. Seed for reproducibility.
#'
#' @usage do.subsample(dat, targets, samp.col, min.per, seed)
#'
#' @examples
#' # Subsample 10,000 cells randomly from the total dataset
#' sub.dat <- Spectre::do.subsample(
#'   dat = Spectre::demo.start,
#'   targets = 10000
#' )
#'
#' # Subsample based on the sample with the smallest number of cells
#' sub.dat.sample <- Spectre::do.subsample(
#'   dat = Spectre::demo.start,
#'   divide.by = "FileName",
#'   min.per = TRUE
#' )
#'
#' @author
#' Thomas Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#' Felix Marsh-Wakefield, \email{felix.marsh-wakefield@@sydney.edu.au}
#'
#' @import data.table
#'
#' @export do.subsample

do.subsample <- function(dat,
                         targets,
                         divide.by = NULL,
                         min.per = FALSE,
                         seed = 42) {

  # require: data.table
  
  ## Test

  # dat <- Spectre::demo.clustered
  # targets <- data.table('Group' = c('WNV', 'Mock'),
  #                       'Targets' = c(10000, 1000))
  # divide.by = 'Group'
  # min.per = FALSE
  # seed = 42

  ## Setup
  dat <- as.data.table(dat)

  ## IF random
  if (is.null(divide.by)) {
    set.seed(seed)
    subsample.res <- dat[sample(1:nrow(dat), targets), ]
    subsample.res <- as.data.table(subsample.res)
    return(subsample.res)
  }

  ## IF per.sample
  if (!is.null(divide.by)) {
    if (min.per == FALSE) {

      # Create list of unique sample names
      if (is.vector(targets)) {
        sample.list <- unique(dat[, divide.by, with = FALSE])
        sample.list <- sample.list[[1]]
      }

      if (is.data.frame(targets)) {
        sample.list <- targets[[1]]
        targets <- targets[[2]]
      }

      # Create res data.frame
      subsample.res <- data.frame()

      # Loop
      for (i in c(1:length(sample.list))) {
        nam <- sample.list[i]
        nsub <- targets[i]

        data.temp <- dat[dat[[divide.by]] == nam, ]

        nrow(data.temp)
        set.seed(seed)
        data.temp <- data.temp[sample(1:nrow(data.temp), nsub), ]
        nrow(data.temp)
        subsample.res <- rbind(subsample.res, data.temp)
      }

      dim(subsample.res)
      # assign("subsample.res", subsample.res, envir = globalenv())
      subsample.res <- as.data.table(subsample.res)
      return(subsample.res)
    }

    ## IF min.per.sample
    if (min.per == TRUE) {
      # Create list of unique sample names
      sample.list <- unique(dat[, divide.by, with = FALSE])
      sample.list <- sample.list[[1]]
      sample.list

      # nrow.check = list()
      # for(i in c(1:(length(DataList)))){nrow.check[[i]] <- nrow(DataList[[i]])}
      # DownSampleTargets <- c(rep(nrow.check[[which.min(nrow.check)]], each=length(unique(AllSampleNos))))

      # min(data.frame(table(dat[[divide.by]]))$Freq) #calculates count of each parameter (divide.by) in data (dat), selecting the minimum number
      # Sets downsample target to be the same for each sample, based on whichever has the smallest number of cells
      targets <- c(rep(min(data.frame(table(dat[[divide.by]]))$Freq), each = length(unique(dat[[divide.by]]))))

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
      # assign("subsample.res", subsample.res, envir = globalenv())
      subsample.res <- as.data.table(subsample.res)
      return(subsample.res)
    }
  }
}
