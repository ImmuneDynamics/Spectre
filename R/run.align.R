#' run.align - Run the alignment model on a target data.table
#'
#' This function allows you to prepare reference data ahead of performing batch alignment, using on of the CytofBatchAdjust functions
#'
#' @usage run.align()
#'
#' @param ref.dat NO DEFAULT. A data.table consisting of the 'refernece' data you will use to train the alignment algorithm
#' @param target.dat NO DEFAULT.A data.table consisting of the 'target' data you will use to align the data
#' @param batch.col NO DEFAULT. Column name of the data.table that contains batch labels
#' @param align.cols NO DEFAULT. Vector of column names to align.
#' @param method DEFAULT = 'quantile'. Can be 'quantile', 'SD', or a percentile indicated as '50p' (50th percentile) or '95p' (95th percentile) etc.
#' @param append.name DEFAULT = "_aligned". This will be appended to the column names containing the aligned data
#' @param dir DEFAULT = getwd(). Sets the working directory to operate from. Because this function involves some reading/writing of files, it's best to set this to somewhere static in case the active working directory moves to a subfolder, and then doesn't return because the function runs into an error.
#' @param mem.ctrl DEFAULT = TRUE. Allows the function to clear held memory on occasion.
#'
#' @return Returns a data.table with aligned data added in new columns.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references Ashhurst, T. M., et al. (2019). \url{https://www.ncbi.nlm.nih.gov/pubmed/31077106}
#'
#' @examples
#' cell.dat <- run.align()
#'
#'
#' @export

run.align <- function(ref.dat,
                      target.dat,
                      batch.col,
                      align.cols,
                      cluster.col = NULL,
                      goal = NULL,
                      method = "quantile",
                      append.name = "_aligned",
                      dir = getwd(),
                      mem.ctrl = TRUE) {

  ### Notes

  # 95p | SD | quantile
  # quantile: quantile normalization
  # SD: scaling to reference batch standard deviation
  # 50p: scaling to reference batch 50th percentile (median)
  # 95p: scaling to reference batch 95th percentile
  # Batches may be scaled to an user-defined percentile by specifying any number (1-100) followed by the letter 'p'. For example method="80p" would scale channels to the 80th percentile of the reference batch.

  # Make batch 1 the average of the other batches?
  # Else set a goal -- make that the first batch...

  ### Demo

  # setwd("/Users/thomasa/Desktop/")
  # dir.create("Align test")
  # setwd("Align test")
  # dir <- getwd()
  #
  # ref.dat <- Spectre::demo.batches.1
  # set.seed(42)
  # nr <- sample(nrow(ref.dat))
  # ref.dat <- ref.dat[nr,]
  #
  # target.dat <- Spectre::demo.batches.1
  # set.seed(42)
  # nr <- sample(nrow(target.dat))
  # target.dat <- target.dat[nr,]
  #
  # batch.col <- 'Batch'
  # align.cols <- names(Spectre::demo.batches.1)[c(1:15)]
  # method <- "95p"
  # cluster.col <- NULL
  # goal <- NULL

  # cluster.col <- 'Population'
  # goal <- 'B'

  ### Packages

  # require data.table, flowCore
  
  ### Setup

  # TODO: what is this? why in the beginning?
  unlink("Temp -- cytofbatchadjust", recursive = TRUE)

  start.dat <- target.dat

  target.dat$ADJUST_ORDER <- c(1:nrow(target.dat))

  # ref.dat$ADJUST_ORDER <- c(1:nrow(ref.dat))
  # target.dat$ADJUST_ORDER <- c(1:nrow(target.dat))

  ### Check batch entries

  message("Alignment - setting up")

  REF.BATCHES <- sort(unique(ref.dat[[batch.col]]))
  TRG.BATCHES <- sort(unique(target.dat[[batch.col]]))

  for (i in TRG.BATCHES) {
    # i <- TRG.BATCHES[[1]]

    if (!any(REF.BATCHES == i)) {
      stop("One of the batches in your target data is not represented in the reference data")
    }
  }

  rm(REF.BATCHES)
  rm(TRG.BATCHES)

  all.batches <- sort(unique(ref.dat[[batch.col]]))
  all.batches

  ## Sort to put goal batch at front (if 'goal' is set)

  if (!is.null(goal)) {
    all.batches <- all.batches[c(which(all.batches == goal), which(all.batches != goal))]
  }

  all.batch.nums <- c(1:length(all.batches))
  all.batch.nums

  batch.tb <- cbind(all.batches, all.batch.nums)
  batch.tb <- as.data.table(batch.tb)
  names(batch.tb) <- c(batch.col, "ADJUST_BATCH_NUMBERS")
  batch.tb

  ###

  ref.dat <- do.add.cols(ref.dat, batch.col, batch.tb, batch.col, show.status = FALSE)
  target.dat <- do.add.cols(target.dat, batch.col, batch.tb, batch.col, show.status = FALSE)

  batch.col <- "ADJUST_BATCH_NUMBERS"

  ref.dat[[batch.col]] <- as.numeric(ref.dat[[batch.col]])
  target.dat[[batch.col]] <- as.numeric(target.dat[[batch.col]])

  ### Sort

  setorderv(ref.dat, "ADJUST_BATCH_NUMBERS")
  setorderv(target.dat, "ADJUST_BATCH_NUMBERS")

  target.brcd <- target.dat[, "ADJUST_ORDER", with = FALSE]
  target.brcd

  gc()

  # ### REFERENCE DAT - Add per batch barcode
  #
  #     ref.dat$PER_BATCH_BARCODE <- 'EMPTY'
  #
  #     for(i in unique(ref.dat[['ADJUST_BATCH_NUMBERS']])){
  #       # i <- unique(ref.dat$ADJUST_BATCH_NUMBERS)[[1]]
  #       ref.dat[ref.dat[['ADJUST_BATCH_NUMBERS']] == i,]$PER_BATCH_BARCODE <- c(1:nrow(ref.dat[ref.dat[['ADJUST_BATCH_NUMBERS']] == i,]))
  #     }
  #
  #     ref.dat$PER_BATCH_BARCODE <- as.numeric(ref.dat$PER_BATCH_BARCODE)
  #
  # ### TARGET DAT - Add per batch barcode
  #
  #     target.dat$PER_BATCH_BARCODE <- 'EMPTY'
  #
  #     for(i in unique(target.dat[['ADJUST_BATCH_NUMBERS']])){
  #       # i <- unique(target.dat$ADJUST_BATCH_NUMBERS)[[1]]
  #       target.dat[target.dat[['ADJUST_BATCH_NUMBERS']] == i,]$PER_BATCH_BARCODE <- c(1:nrow(target.dat[target.dat[['ADJUST_BATCH_NUMBERS']] == i,]))
  #     }
  #
  #     target.dat$PER_BATCH_BARCODE <- as.numeric(target.dat$PER_BATCH_BARCODE)
  #
  # # ref.dat[['ADJUST_ORDER']] <- as.numeric(ref.dat[['ADJUST_ORDER']])
  # # target.dat[['ADJUST_ORDER']] <- as.numeric(target.dat[['ADJUST_ORDER']])
  #
  # gc()

  ### Setup a temp directory

  setwd(dir)
  getwd()

  unlink("Temp -- cytofbatchadjust", recursive = TRUE)
  dir.create("Temp -- cytofbatchadjust")
  setwd("Temp -- cytofbatchadjust")

  temp.dir <- getwd()

  ########################################################################################################
  ################################# Whole dataset (i.e. no clusters) #####################################
  ########################################################################################################

  if (is.null(cluster.col)) {

    ### Write REFERENCE and TARGET FCS files to disk

    message("Alignment - preparing reference and target files")

    ##
    ref.dat

    for (i in unique(ref.dat[[batch.col]])) {
      # i <- unique(ref.dat[[batch.col]])[[1]]

      temp <- ref.dat[ref.dat[[batch.col]] == i, ]
      temp <- temp[, c(align.cols), with = FALSE]
      temp$PER_BATCH_BARCODE <- c(1:nrow(temp))

      # temp$`PER_BATCH_BARCODE`
      # str(temp)

      ## Write FCS
      metadata <- data.frame(name = dimnames(temp)[[2]], desc = paste(dimnames(temp)[[2]]))

      ## Create FCS file metadata - ranges, min, and max settings
      metadata$range <- apply(apply(temp, 2, range), 2, diff)
      metadata$minRange <- apply(temp, 2, min)
      metadata$maxRange <- apply(temp, 2, max)

      ## Create flowframe with tSNE data
      temp <- new("flowFrame", exprs = as.matrix(temp), parameters = Biobase::AnnotatedDataFrame(metadata))

      ## Save flowframe as .fcs file -- save data (with new tSNE parameters) as FCS
      write.FCS(temp, paste0("BATCH_", i, "_", "Ref.fcs"))

      rm(temp)
      rm(metadata)
      rm(i)
      gc()
    }

    rm(ref.dat)
    gc()

    ##
    target.dat

    for (i in unique(target.dat[[batch.col]])) {
      # i <- unique(target.dat[[batch.col]])[[1]]

      temp <- target.dat[target.dat[[batch.col]] == i, ]
      temp <- temp[, c(align.cols), with = FALSE]
      temp$PER_BATCH_BARCODE <- c(1:nrow(temp))

      ## Write FCS
      metadata <- data.frame(name = dimnames(temp)[[2]], desc = paste(dimnames(temp)[[2]]))

      ## Create FCS file metadata - ranges, min, and max settings
      metadata$range <- apply(apply(temp, 2, range), 2, diff)
      metadata$minRange <- apply(temp, 2, min)
      metadata$maxRange <- apply(temp, 2, max)

      ## Create flowframe with tSNE data
      temp <- new("flowFrame", exprs = as.matrix(temp), parameters = Biobase::AnnotatedDataFrame(metadata))

      ## Save flowframe as .fcs file -- save data (with new tSNE parameters) as FCS
      write.FCS(temp, paste0("BATCH_", i, "_.fcs"))

      rm(temp)
      rm(metadata)
      rm(i)
      gc()
    }

    # target.dat <- as.data.table(target.dat[,c(batch.col, 'PER_BATCH_BARCODE'),with = FALSE])
    # names(target.dat) <- 'PER_BATCH_BARCODE'

    gc()

    ### Write .txt. file of 'channels to adjust

    write(align.cols, "ChannelsToAdjust.txt")

    ### Run

    message("Alignment - applying alignment process")

    unlink("OUTPUT", recursive = TRUE)
    # dir.create("OUTPUT")

    BatchAdjust(
      basedir = ".",
      outdir = "OUTPUT",
      channelsFile = "ChannelsToAdjust.txt",
      batchKeyword = "BATCH_",
      anchorKeyword = "Ref",
      method = method,
      transformation = FALSE,
      addExt = NULL,
      plotDiagnostics = FALSE
    )

    ### Read in output

    setwd(dir)
    setwd("Temp -- cytofbatchadjust")
    setwd("OUTPUT")

    output.files <- list.files(getwd(), ".fcs")
    output.files <- output.files[!grepl("_Ref.fcs", output.files)]

    res.list <- list()

    for (i in output.files) {
      # i <- output.files[[1]]

      temp <- exprs(flowCore::read.FCS(i,
        transformation = FALSE,
        # truncate_max_range = FALSE
      ))

      temp <- temp[1:nrow(temp), 1:ncol(temp)]
      temp <- as.data.table(temp)

      i <- gsub("BATCH_", "", i)
      i <- gsub("_.fcs", "", i)
      i <- gsub("_", "", i)

      temp[["BATCH_CHECK"]] <- i
      temp[["BATCH_CHECK"]] <- as.numeric(temp[["BATCH_CHECK"]])

      # setorderv(temp, 'PER_BATCH_BARCODE')

      res.list[[i]] <- temp

      rm(i)
      rm(temp)
    }

    res <- rbindlist(res.list, fill = TRUE)

    setorderv(res, "PER_BATCH_BARCODE")
    setorderv(res, "BATCH_CHECK")

    res <- cbind(target.brcd, res)

    setorderv(res, "ADJUST_ORDER")

    ### Checks

    # if(any(target.dat$ADJUST_BATCH_NUMBERS != res$BATCH_CHECK)){
    #   message('WARNING: Row barcode for adjusted data did not match original target data')
    # }
    #
    # if(any(target.dat$PER_BATCH_BARCODE != res$PER_BATCH_BARCODE)){
    #   message('WARNING: Row barcode for adjusted data did not match original target data')
    # }
    # if(any(target.dat$ADJUST_ORDER != res$ADJUST_ORDER)){
    #   message('WARNING: Row barcode for adjusted data did not match original target data')
    # }

    res <- res[, align.cols, with = FALSE]

    names(res) <- paste0(names(res), append.name)

    ### Finalise

    message("Alignment - finalising")

    setwd(dir)
    unlink("Temp -- cytofbatchadjust", recursive = TRUE)

    start.dat <- cbind(start.dat, res)
    return(start.dat)
  }

  ########################################################################################################
  ########################################### CLUSTERS ###################################################
  ########################################################################################################

  if (!is.null(cluster.col)) {
    stop("The use of 'cluster.col' is not active yet")
  }

  # if(!is.null(cluster.col)){
  #
  #   message('Alignment - preparing reference and target files BY CLUSTER')
  #   all.res <- list()
  #
  #   ### Loop
  #
  #   for(a in sort(unique(ref.dat[[cluster.col]]))){
  #     # a <- sort(unique(ref.dat[[cluster.col]]))[[1]]
  #
  #     ### Prep
  #
  #         message(paste0(" -- running ", "'", a, "'"))
  #
  #         setwd(temp.dir)
  #         dir.create(a)
  #         setwd(a)
  #
  #         clust.dir <- getwd()
  #
  #         ref.clust <- ref.dat[ref.dat[[cluster.col]] == a,]
  #         target.clust <- target.dat[target.dat[[cluster.col]] == a,]
  #
  #     ### Write REFERENCe files to disk
  #
  #         for(i in unique(ref.clust[[batch.col]])){
  #           # i <- unique(ref.dat[[batch.col]])[[1]]
  #
  #           temp <- ref.clust[ref.clust[[batch.col]] == i,]
  #           temp <- temp[, c(align.cols), with = FALSE]
  #           temp$PER_BATCH_BARCODE <- c(1:nrow(temp))
  #
  #           ## Write FCS
  #           metadata <- data.frame(name=dimnames(temp)[[2]], desc=paste(dimnames(temp)[[2]]))
  #
  #           ## Create FCS file metadata - ranges, min, and max settings
  #           metadata$range <- apply(apply(temp,2,range),2,diff)
  #           metadata$minRange <- apply(temp,2,min)
  #           metadata$maxRange <- apply(temp,2,max)
  #
  #           ## Create flowframe with tSNE data
  #           temp <- new("flowFrame", exprs=as.matrix(temp), parameters=Biobase::AnnotatedDataFrame(metadata))
  #
  #           ## Save flowframe as .fcs file -- save data (with new tSNE parameters) as FCS
  #           write.FCS(temp, paste0("BATCH_", i, "_", "Ref.fcs"))
  #
  #           rm(temp)
  #           rm(metadata)
  #           rm(i)
  #           gc()
  #         }
  #
  #         gc()
  #
  #     ### Write TARGET data to disk
  #
  #         for(i in unique(target.clust[[batch.col]])){
  #           # i <- unique(target.dat[[batch.col]])[[1]]
  #
  #           temp <- target.clust[target.clust[[batch.col]] == i,]
  #           temp <- temp[, c(align.cols), with = FALSE]
  #           temp$PER_BATCH_BARCODE <- c(1:nrow(temp))
  #
  #           ## Write FCS
  #           metadata <- data.frame(name=dimnames(temp)[[2]], desc=paste(dimnames(temp)[[2]]))
  #
  #           ## Create FCS file metadata - ranges, min, and max settings
  #           metadata$range <- apply(apply(temp,2,range),2,diff)
  #           metadata$minRange <- apply(temp,2,min)
  #           metadata$maxRange <- apply(temp,2,max)
  #
  #           ## Create flowframe with tSNE data
  #           temp <- new("flowFrame", exprs=as.matrix(temp), parameters=Biobase::AnnotatedDataFrame(metadata))
  #
  #           ## Save flowframe as .fcs file -- save data (with new tSNE parameters) as FCS
  #           write.FCS(temp, paste0("BATCH_", i, "_.fcs"))
  #
  #           rm(temp)
  #           rm(metadata)
  #           rm(i)
  #           gc()
  #         }
  #
  #         # target.dat <- as.data.table(target.dat[,c(batch.col, 'PER_BATCH_BARCODE'),with = FALSE])
  #         # names(target.dat) <- 'PER_BATCH_BARCODE'
  #
  #         gc()
  #
  #     ### Write .txt. file of 'channels to adjust
  #
  #         write(align.cols, 'ChannelsToAdjust.txt')
  #
  #     ### Run
  #
  #         message('Alignment - applying alignment process')
  #
  #         unlink("OUTPUT", recursive = TRUE)
  #         #dir.create("OUTPUT")
  #
  #         BatchAdjust(
  #           basedir=".",
  #           outdir="OUTPUT",
  #           channelsFile = "ChannelsToAdjust.txt",
  #           batchKeyword="BATCH_",
  #           anchorKeyword = "Ref",
  #           method=method,
  #           transformation=FALSE,
  #           addExt=NULL,
  #           plotDiagnostics=FALSE)
  #
  #     ### Read in output
  #
  #         setwd(clust.dir)
  #         setwd("OUTPUT")
  #
  #         output.files <- list.files(getwd(), '.fcs')
  #         output.files <- output.files[!grepl('_Ref.fcs', output.files)]
  #
  #         res.list <- list()
  #
  #         for(i in output.files){
  #           # i <- output.files[[1]]
  #
  #           temp <- exprs(flowCore::read.FCS(i,
  #                                            transformation = FALSE,
  #                                            #truncate_max_range = FALSE
  #           ))
  #
  #           temp <- temp[1:nrow(temp),1:ncol(temp)]
  #           temp <- as.data.table(temp)
  #
  #           i <- gsub("BATCH_", "", i)
  #           i <- gsub("_.fcs", "", i)
  #           i <- gsub("_", "", i)
  #
  #           temp[['BATCH_CHECK']] <- i
  #           temp[['BATCH_CHECK']] <- as.numeric(temp[['BATCH_CHECK']])
  #
  #           # setorderv(temp, 'PER_BATCH_BARCODE')
  #
  #           res.list[[i]] <- temp
  #
  #           rm(i)
  #           rm(temp)
  #         }
  #
  #         res <- rbindlist(res.list, fill = TRUE)
  #
  #         setorderv(res, 'PER_BATCH_BARCODE')
  #         setorderv(res, 'BATCH_CHECK')
  #
  #         all.res[[a]] <- res
  #       }
  #
  #   ###
  #
  #
  #           res <- cbind(target.brcd, res)
  #
  #           setorderv(res, 'ADJUST_ORDER')
  #
  #           ### Checks
  #
  #           # if(any(target.dat$ADJUST_BATCH_NUMBERS != res$BATCH_CHECK)){
  #           #   message('WARNING: Row barcode for adjusted data did not match original target data')
  #           # }
  #           #
  #           # if(any(target.dat$PER_BATCH_BARCODE != res$PER_BATCH_BARCODE)){
  #           #   message('WARNING: Row barcode for adjusted data did not match original target data')
  #           # }
  #           # if(any(target.dat$ADJUST_ORDER != res$ADJUST_ORDER)){
  #           #   message('WARNING: Row barcode for adjusted data did not match original target data')
  #           # }
  #
  #           res <- res[, align.cols, with = FALSE]
  #           names(res) <- paste0(names(res), append.name)
  #
  #
  #
  #
  #   ### Finalise
  #
  #
  #
  #
  #       message('Alignment - finalising')
  #
  #       setwd(dir)
  #       unlink("Temp -- cytofbatchadjust", recursive = TRUE)
  #
  #
  #       #
  #       #
  #       #
  #
  #
  #
  #
  #
  #       start.dat <- cbind(start.dat, res)
  #       return(start.dat)
  #
  # }
}
