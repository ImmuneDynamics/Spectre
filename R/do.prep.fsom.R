#' do.prep.fsom - Function to prepare a FlowSOM (fsom) object for alignment with CytoNorm
#'
#' This function allows you to prepare a FlowSOM (fsom) object for alignment with CytoNorm.
#'
#' @usage do.prep.fsom(dat, use.cols, sample.col, xdim, ydim, nClus, nCells, scale, seed)
#'
#' @param dat NO DEFAULT. data.table of the dataset you wish to run in FlowSOM.
#' @param use.cols NO DEFAULT. Vector of character etnries. Columns to use to calculate FlowSOM.
#' @param sample.col NO DEFAULT. Character, name of the column denoting samples
#' @param xdim DEFAULT = 5. Size of X-axis of FlowSOM grid.
#' @param xdim DEFAULT = 5. Size of Y-axis of FlowSOM grid.
#' @param nClus DEFAULT = 10. Number of metaclusters. If set to 1, will map all cells to a single metacluster
#' @param nCells DEFAULT = NULL. Number of cells to downsample to, otherwise 'NULL' will use all available cells
#' @param scale DEFAULT = FALSE. Do features need to be scaled?
#' @param seed DEFAULT = 42. Seed for reproducibility.
#'
#' @return Returns a FlowSOM object.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references Ashhurst, T. M., et al. (2019). \url{https://www.ncbi.nlm.nih.gov/pubmed/31077106}
#'
#' @examples
#' ref.fsom <- do.prep.fsom()
#'
#' @export

do.prep.fsom <- function(dat,
                         use.cols,
                         sample.col,
                         batch.col,
                         xdim = 5,
                         ydim = 5,
                         nClus = 10,
                         nCells = NULL,
                         scale = FALSE,
                         seed = 42){

  ### Check that necessary packages are installed
      if(!is.element('Spectre', installed.packages()[,1])) stop('Spectre is required but not installed')
      if(!is.element('data.table', installed.packages()[,1])) stop('data.table is required but not installed')
      if(!is.element('CytoNorm', installed.packages()[,1])) stop('CytoNorm is required but not installed')
      if(!is.element('flowCore', installed.packages()[,1])) stop('flowCore is required but not installed')
      if(!is.element('Biobase', installed.packages()[,1])) stop('Biobase is required but not installed')

  ### Require packages
      require(Spectre)
      require(data.table)
      require(CytoNorm)
      require(flowCore)
      require(Biobase)

  ### Test data

      # dat <- Spectre::demo.clustered
      # dat$FlowSOM_cluster <- NULL
      # dat$FlowSOM_metacluster <- NULL
      #
      # as.matrix(names(dat))
      # use.cols <- names(dat)[c(11:19)]
      #
      # sample.col <- "Sample"
      # batch.col <- "Batch"
      #
      # xdim = 5
      # ydim = 5
      # nClus = 1
      # nCells = NULL
      # scale = FALSE
      # seed = 1

  ### Single cluster preferences

      if(nClus == 1){
        nClus <- 5
        one.clust <- TRUE
      } else {
        one.clust <- FALSE
      }

  ### Create 'temp' folder
  starting.dir <- getwd()

  setwd(starting.dir)
  dir.create("tmp-cytonorm-fsom")
  setwd("tmp-cytonorm-fsom")

  ### Write files

  message("Step 1/4 - Splitting files for use with original FlowSOM function")

  dat.list <- unique(dat[[sample.col]])
  dat.list

  # dat.list <- dat.list[order(dat.list)]
  # dat.list

  for(i in c(1:length(dat.list))){
    # i <- 1

    a <- dat.list[[i]]

    temp <- dat[dat[[sample.col]] == a,]

    write.files(temp,
                file.prefix = a,
                write.csv = FALSE,
                write.fcs = TRUE)
  }

  files <- list.files(getwd(), ".fcs")
  files

  ## Get batch for each file

  batches <- vector()

  for(a in files){
    #a <- "export_TA203-2.fcs"
    a <- gsub(".fcs", "", a)
    temp <- dat[dat[[sample.col]] == a,]
    res <- temp[[batch.col]][1]

    ## Should add a test to check for multiple batch entries per sample

    batches <- cbind(batches, res)
  }

  batches <- as.vector(batches)
  batches

  ### Run FlowSOM on ref data (--> save as ref.ff, also save as ref.dat)

  message("Step 2/4 - Running FlowSOM")

  fsom <- prepareFlowSOM(files,
                         colsToUse = use.cols,
                         nCells = nCells, #########################
                         FlowSOM.params = list(xdim = xdim,
                                               ydim = ydim,
                                               nClus = nClus,
                                               scale = scale),
                         seed = seed)

  if(nrow(fsom$FlowSOM$data) != nrow(dat)){
    stop("Error - the numer of rows (cells) is different in the starting dataset and the FlowSOM prepared dataset")
  }

  #nrow(fsom$FlowSOM$data)
  #nrow(dat)

  setwd(starting.dir)

  ### Data, clusters, and metaclusters
  # fsom$FlowSOM$data # data
  # fsom$FlowSOM$map$mapping[,1] # 1* clusters
  # unique(fsom$FlowSOM$map$mapping[,1]) # unique 1* clusters
  # fsom$metaclustering # metaclusters

  ### Mapping
  # fsom$metaclustering[fsom$FlowSOM$map$mapping[,1]]

  ### Merge outputs into Spectre format

  message("Step 3/4 - Preparing FlowSOM object and results data.table")

  if(one.clust == TRUE){
    length(fsom$metaclustering)
    fsom$metaclustering <- rep(1, length(fsom$metaclustering))
  }

  A <- fsom$FlowSOM$data
  B <- fsom$FlowSOM$map$mapping[,1]
  C <- fsom$metaclustering[fsom$FlowSOM$map$mapping[,1]]

  if(nrow(A) != length(B)){
    stop("Error - the numer of rows (cells) is different in the data, clusters, and/or metaclusters")
  }

  if(nrow(A) != length(C)){
    stop("Error - the numer of rows (cells) is different in the data, clusters, and/or metaclusters")
  }

  if(length(B) != length(C)){
    stop("Error - the numer of rows (cells) is different in the data, clusters, and/or metaclusters")
  }

  # nrow(A)
  # length(B)
  # length(C)

  ### Merge

  fsom.dt <- as.data.table(A)
  fsom.dt <- cbind(A, prep.fsom.cluster = B, prep.fsom.metacluster = C)
  fsom.dt <- as.data.table(fsom.dt)
  # str(fsom.dt)

  ### Results
  setwd(starting.dir)
  unlink("tmp-cytonorm-fsom", recursive = TRUE)

  files <- gsub(".fcs", "", files)

  fsom$files <- files
  fsom$filenums <- unique(fsom.dt$File)
  fsom$batches <- batches
  fsom$features <- use.cols

  res <- named.list(fsom, fsom.dt)

  message("Step 4/4 - FlowSOM complete")
  return(res)
}

