#' prep.cytonorm - Prepare reference data into a FlowSOM object
#'
#' This function allows you to prepare reference data ahead of performing batch alignment.
#'
#' @usage prep.cytonorm()
#'
#' @param dat NO DEFAULT. A data.table consisting of the 'refernece' data you will use to train the alignment algorithm
#' @param cellular.cols NO DEFAULT. A vector of column names from the data.table that contain 'cellular' markers
#' @param cluster.cols NO DEFAULT. A vector of column names from the data.table that contain markers you wish to use for clusteirng
#' @param batch.col NO DEFAULT. Name of the column that contains batch names
#' @param sample.col DEFAULT = NULL. Name of the column that contains sample names
#' @param dir DEFAULT = getwd(). Sets the working directory to operate from. Because this function involves some reading/writing of files, it's best to set this to somewhere static in case the active working directory moves to a subfolder, and then doesn't return because the function runs into an error.
#' @param xdim DEFAULT = 5. Size of X-axis of FlowSOM grid.
#' @param ydim DEFAULT = 5. Size of Y-axis of FlowSOM grid.
#' @param meta.k DEFAULT = 10. Number of metaclusters. If set to 1, will map all cells to a single metacluster
#' @param seed DEFAULT = 42. Seed for reproducibility.
#' @param mem.ctrl DEFAULT = TRUE. Allows the function to clear held memory on occasion.
#'
#' @return Returns an object which represents the alignment model. In this preparation stage, it contains the FlowSOM object containing the reference data. The 'train.align' function can then be used to calculate the conversions between batches for each metacluser.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references Ashhurst, T. M., et al. (2019). \url{https://www.ncbi.nlm.nih.gov/pubmed/31077106}
#'
#' @examples
#' align.model <- prep.cytonorm()
#'
#' @import data.table
#'
#' @export

prep.cytonorm <- function(dat,
                          cellular.cols,
                          cluster.cols,
                          batch.col,
                          sample.col = NULL,
                          dir = getwd(),
                          xdim = 5,
                          ydim = 5,
                          meta.k = 10, # can be 1, or >3
                          seed = 42,
                          mem.ctrl = TRUE) {

  ### Check that necessary packages are installed
  if (!is.element("Spectre", installed.packages()[, 1])) stop("Spectre is required but not installed")
  if (!is.element("data.table", installed.packages()[, 1])) stop("data.table is required but not installed")
  if (!is.element("CytoNorm", installed.packages()[, 1])) stop("CytoNorm is required but not installed")
  if (!is.element("flowCore", installed.packages()[, 1])) stop("flowCore is required but not installed")
  if (!is.element("Biobase", installed.packages()[, 1])) stop("Biobase is required but not installed")

  ### Require packages
  require(data.table)
  require(CytoNorm)
  # require(flowCore)
  require(Biobase)

  ### Directory setup

  setwd(dir)
  starting.dir <- getwd()
  message("Working directory is '", starting.dir, "'")

  ### Initial checks

  ## Are multiple samples present in each batch ?

  ## Does 'dat' have more than 0 rows

  if (nrow(dat) == 0) {
    setwd(dir)
    stop("Error -- your 'dat' has no rows, please check")
  }

  ## Are all 'cluster.cols' present in 'cellular.cols'


  ## meta.k value
  if (meta.k == 2) {
    setwd(dir)
    stop("Error -- cannot create '2' metaclusters -- please set meta.k to '1' or a value >3")
  }

  all.cols <- unique(c(cellular.cols, cluster.cols))
  all.cols

  value <- dat[, all.cols, with = FALSE]

  if (isFALSE(all(sapply(value, is.numeric)))) {
    message("It appears that one column in your dataset is non numeric")
    print(sapply(value, is.numeric))
    setwd(dir)
    stop("do.asinh stopped")
  }

  rm(value)

  if (mem.ctrl == TRUE) {
    gc()
  }

  ### Initial preparation

  if (!is.null(sample.col)) {
    dat <- as.data.table(dat)
    dat <- dat[, c(sample.col, batch.col, all.cols), with = FALSE]
  }

  if (is.null(sample.col)) {
    dat <- as.data.table(dat)
    dat <- dat[, c(batch.col, all.cols), with = FALSE]
  }

  ### Metacluster settings
  # Essentially if meta.k = 1, then run with 5 metaclusters and replace the metacluster values at the end.
  # If meta.k is set to >1, then run as per normal

  if (meta.k == 1) {
    meta.k <- 5
    one.clust <- TRUE
  } else {
    one.clust <- FALSE
  }

  ### Create temp folder

  ## Set working directory
  setwd(starting.dir)

  ## Remove any previous 'tmp' folder
  unlink("tmp-cytonorm-fsom", recursive = TRUE)

  ## Create new 'tmp' folder
  dir.create("tmp-cytonorm-fsom", showWarnings = FALSE)
  setwd("tmp-cytonorm-fsom")

  ## Set working directory
  setwd(starting.dir)

  ### Convert sample.col into a numerical string

  if (!is.null(sample.col)) {
    samps <- unique(dat[[sample.col]])
    samps.num <- c(1:length(samps))

    smp.tb <- cbind(samps, samps.num)
    x <- as.data.table(dat[[sample.col]])
    names(x) <- sample.col

    x <- do.add.cols(x, sample.col, smp.tb, "samps")
    x <- x$samps.num

    dat[[sample.col]] <- x
    dat[[sample.col]] <- as.numeric(dat[[sample.col]])

    rm(x)
    rm(samps)
    rm(samps.num)
  }


  ### Write files (one per batch)

  message("Step 1/4 - Splitting files for use with original FlowSOM function")

  dat.list <- unique(dat[[batch.col]])
  dat.list

  setwd(starting.dir)
  dir.create("tmp-cytonorm-fsom", showWarnings = FALSE)
  setwd("tmp-cytonorm-fsom")

  for (i in c(1:length(dat.list))) {
    # i <- 1
    a <- dat.list[[i]]
    temp <- dat[dat[[batch.col]] == a, ]

    if (!is.null(sample.col)) {
      temp <- temp[, c(sample.col, all.cols), with = FALSE]
    }

    if (is.null(sample.col)) {
      temp <- temp[, c(all.cols), with = FALSE]
    }

    write.files(temp,
      file.prefix = a,
      write.csv = FALSE,
      write.fcs = TRUE
    )

    rm(i)
    rm(a)
  }

  rm(dat.list)

  ## This is key -- here we get a list (sorted by reading file from disk)
  files <- list.files(getwd(), ".fcs")
  files

  file.nums <- c(1:length(files))

  setwd(starting.dir)

  if (mem.ctrl == TRUE) {
    gc()
  }

  ### Run FlowSOM

  message("Step 2/4 - Running FlowSOM")

  setwd(starting.dir)
  dir.create("tmp-cytonorm-fsom", showWarnings = FALSE)
  setwd("tmp-cytonorm-fsom")

  # message(paste0(" -- files are :", list.files(getwd(), '.fcs')))

  fsom <- prepareFlowSOM(files,
    colsToUse = cluster.cols,
    nCells = NULL, # Any subsampling is done prior to using this function
    FlowSOM.params = list(
      xdim = xdim,
      ydim = ydim,
      nClus = meta.k,
      scale = FALSE
    ), # Any transformations are done prior to using this function
    seed = seed
  )

  setwd(starting.dir)
  unlink("tmp-cytonorm-fsom", recursive = TRUE)

  if (nrow(fsom$data) != nrow(dat)) {
    stop("Error - the numer of rows (cells) is different in the starting dataset and the FlowSOM prepared dataset")
  }

  if (mem.ctrl == TRUE) {
    gc()
  }

  ### Pull results together

  message("Step 3/4 - Preparing FlowSOM object and results data.table")

  if (one.clust == TRUE) {
    length(fsom$metaclustering)
    fsom$metaclustering <- rep(1, length(fsom$metaclustering))
  }

  A <- fsom$data
  B <- fsom$map$mapping[, 1]
  C <- fsom$metaclustering[fsom$map$mapping[, 1]]

  dt <- as.data.table(A)
  dt <- cbind(A, prep.fsom.cluster = B, prep.fsom.metacluster = C)
  dt <- as.data.table(dt)

  files <- gsub(".fcs", "", files)

  res <- named.list(fsom, dt, cellular.cols, cluster.cols, files, file.nums)

  rm(fsom)
  rm(dt)

  if (mem.ctrl == TRUE) {
    gc()
  }

  ### Return

  message("Step 4/4 - FlowSOM preparation for alignment model complete")
  setwd(starting.dir)

  return(res)
}
