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
#' @import data.table
#'
#' @export

run.align <- function(ref.dat,
                       target.dat,
                       batch.col,
                       align.cols,
                       
                       method = "quantile", # 'harmony'
                       append.name = "_aligned",
                       dir = getwd(),
                       mem.ctrl = TRUE
){
  
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
      
      # ref.dat <- ref.dat
      # target.dat <- target.dat
      # dir <- working.dir
      # 
      # batch.col <- 'Batches'
      # align.cols <- names(dat)[c(11:16,18)] 
      # method <- "95p"
      
#######################################################################################################
########################################### Harmony methods ###########################################
#######################################################################################################

  if(method == 'harmony'){

  }
  
#######################################################################################################
###################################### CytofBatchAdjust methods #######################################
#######################################################################################################
  
  ### Setup
      
      unlink("Temp -- cytofbatchadjust", recursive = TRUE)
      
      start.dat <- target.dat
      
      ref.dat$ADJUST_ORDER <- c(1:nrow(ref.dat))
      target.dat$ADJUST_ORDER <- c(1:nrow(target.dat))
  
  ### Adjust batches to match requirements
      
      ## Check batch entries
      
      REF.BATCHES <- sort(unique(ref.dat[[batch.col]]))
      TRG.BATCHES <- sort(unique(target.dat[[batch.col]]))
      
      for(i in TRG.BATCHES){
        # i <- TRG.BATCHES[[1]]
        
        if(!any(REF.BATCHES == i)){
          stop("One of the batches in your target data is not represented in the reference data")
        }
      }
      
      rm(REF.BATCHES)
      rm(TRG.BATCHES)
      
      all.batches <- sort(unique(ref.dat[[batch.col]]))
      all.batches
      
      all.batch.nums <- c(1:length(all.batches))
      all.batch.nums
      
      batch.tb <- cbind(all.batches, all.batch.nums)
      batch.tb <- as.data.table(batch.tb)
      names(batch.tb) <- c(batch.col, "ADJUST_BATCH_NUMBERS")
      batch.tb
      
      ## 
      
      ref.dat <- do.add.cols(ref.dat, batch.col, batch.tb, batch.col)
      target.dat <- do.add.cols(target.dat, batch.col, batch.tb, batch.col)
      
      batch.col <- 'ADJUST_BATCH_NUMBERS'
      
      ref.dat[[batch.col]] <- as.numeric(ref.dat[[batch.col]])
      target.dat[[batch.col]] <- as.numeric(target.dat[[batch.col]])
      
  ### Setup a temp directory
  
      setwd(dir)
      getwd()
      
      unlink("Temp -- cytofbatchadjust", recursive = TRUE)
      dir.create("Temp -- cytofbatchadjust")
      setwd("Temp -- cytofbatchadjust")
      
      temp.dir <- getwd()  
      
  ### Write REFERENCE and TARGET FCS files to disk
      
      ##
      ref.dat
      
      for(i in all.batch.nums){
        # i <- all.batch.nums[[2]]
        
        temp <- ref.dat[ref.dat[[batch.col]] == i,]
        temp <- temp[, c('ADJUST_ORDER', align.cols), with = FALSE]
        
        ## Write FCS
        metadata <- data.frame(name=dimnames(temp)[[2]], desc=paste(dimnames(temp)[[2]]))
        
        ## Create FCS file metadata - ranges, min, and max settings
        metadata$range <- apply(apply(temp,2,range),2,diff)
        metadata$minRange <- apply(temp,2,min)
        metadata$maxRange <- apply(temp,2,max)
        
        ## Create flowframe with tSNE data
        temp.ff <- new("flowFrame", exprs=as.matrix(temp), parameters=Biobase::AnnotatedDataFrame(metadata))
        
        ## Save flowframe as .fcs file -- save data (with new tSNE parameters) as FCS
        write.FCS(temp.ff, paste0("BATCH_", i, "_", "Ref.fcs"))
        
        rm(temp)
        rm(temp.ff)
        rm(metadata)
        rm(i)
      }
      
      ## 
      target.dat
      
      for(i in unique(target.dat[[batch.col]])){
        # i <- unique(target.dat[[batch.col]])[[1]]
        
        temp <- target.dat[target.dat[[batch.col]] == i,]
        temp <- temp[, c('ADJUST_ORDER', align.cols), with = FALSE]
        
        ## Write FCS
        metadata <- data.frame(name=dimnames(temp)[[2]], desc=paste(dimnames(temp)[[2]]))
        
        ## Create FCS file metadata - ranges, min, and max settings
        metadata$range <- apply(apply(temp,2,range),2,diff)
        metadata$minRange <- apply(temp,2,min)
        metadata$maxRange <- apply(temp,2,max)
        
        ## Create flowframe with tSNE data
        temp.ff <- new("flowFrame", exprs=as.matrix(temp), parameters=Biobase::AnnotatedDataFrame(metadata))
        
        ## Save flowframe as .fcs file -- save data (with new tSNE parameters) as FCS
        write.FCS(temp.ff, paste0("BATCH_", i, "_.fcs"))
        
        rm(temp)
        rm(temp.ff)
        rm(metadata)
        rm(i)
      }
      
  ### Write .txt. file of 'channels to adjust
  
      write(align.cols, 'ChannelsToAdjust.txt')
  
  ### Run
      
      unlink("OUTPUT", recursive = TRUE)
      #dir.create("OUTPUT")
      
      BatchAdjust(
        basedir=".",
        outdir="OUTPUT",
        channelsFile = "ChannelsToAdjust.txt",
        batchKeyword="BATCH_",
        anchorKeyword = "Ref",
        method=method,
        transformation=FALSE,
        addExt=NULL,
        plotDiagnostics=FALSE)
      
      
  ### Read in output
      
      setwd(dir)
      setwd("Temp -- cytofbatchadjust")
      setwd("OUTPUT")
      
      output.files <- list.files(getwd(), '.fcs')
      output.files <- output.files[!grepl('_Ref.fcs', output.files)]
      
      res.list <- list()
      
      for(i in output.files){
        # i <- output.files[[1]]
        
        temp <- exprs(flowCore::read.FCS(i, transformation = FALSE))
        temp <- temp[1:nrow(temp),1:ncol(temp)]
        temp <- as.data.table(temp)
        
        res.list[[i]] <- temp
        
        rm(i)
        rm(temp)
      }
      
      res <- rbindlist(res.list, fill = TRUE)
      setorderv(res, 'ADJUST_ORDER')
      
  ### Checks
      
      if(any(target.dat$ADJUST_ORDER != res$ADJUST_ORDER)){
        stop('Result rows could not be matched with original target data')
      }
      
      res <- res[, align.cols, with = FALSE]
      
      names(res) <- paste0(names(res), append.name)
  
  ### Finalise
      
      setwd(dir)
      unlink("Temp -- cytofbatchadjust", recursive = TRUE)
      
      start.dat <- cbind(start.dat, res)
      return(start.dat)        
}  


