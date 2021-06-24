#' train.cytonorm - Prepare reference data into a FlowSOM object
#'
#' This function allows you to learn the conversions to use for batch alignment.
#' 
#' @usage train.align()
#'
#' @param dat NO DEFAULT. A data.table consisting of the 'refernece' data you will use to train the alignment algorithm
#' @param cellular.cols NO DEFAULT. A vector of column names from the data.table that contain the markers to be aligned
#' @param method DEFAULT = 'cytonorm'. In future additional methods for batch alignment will be added.
#' @param cytonorm.goal DEFAULT = 'mean'. Target values for alignment. Can be 'mean' for the average of all batches, or a specific batch can be selected.
#' @param cytonorm.nQ DEFAULT = 101. Number of quantiles.
#' @param dir DEFAULT = getwd(). Sets the working directory to operate from. Because this function involves some reading/writing of files, it's best to set this to somewhere static in case the active working directory moves to a subfolder, and then doesn't return because the function runs into an error.
#' @param mem.ctrl DEFAULT = TRUE. Allows the function to clear held memory on occasion.
#'
#' @return Returns the alignment model object, with an additional 'conversions' element added.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references Ashhurst, T. M., et al. (2019). \url{https://www.ncbi.nlm.nih.gov/pubmed/31077106}
#'
#' @examples
#' align.model <- train.cytonorm()
#'
#' @import data.table
#'
#' @export

train.cytonorm <- function(model,
                          align.cols,
                          
                          method = 'cytonorm', # cytonorm or quantile 
                          
                          ## CytoNorm settings
                          cytonorm.goal = "mean",
                          cytonorm.nQ = 101,
                          
                          ## Quantile settings,
                          quantile.min = 0.001,
                          quantile.max = 0.009,
                          
                          ## General settings
                          dir = getwd(),
                          mem.ctrl = TRUE){
  
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
      # require(flowCore)
      require(Biobase)
  
  ### Check the alignment model object
      message("Training alignment - setup")
      
      model$fsom$data
  
  ### Directory setup
      
      setwd(dir)
      starting.dir<- getwd()
      message("Working directory is '", starting.dir, "'")
  
#########################################################################################           
######################################## Quantiles ######################################
#########################################################################################      
      
if(method == 'quantile'){

  stop("Quantile not yet supported")
  
}
      
#########################################################################################           
######################################## CytoNorm #######################################
#########################################################################################      
     
if(method == 'cytonorm'){
      
  ### Some settings
      
      outputDir = "./tmp"
      normMethod.train = QuantileNorm.train
      normParams = list(nQ = cytonorm.nQ, goal = cytonorm.goal) ### Can we do more?? Default = 101 Can also be a batch #
      
      transformList = NULL # in the normal fsom, the default is NULL, but not here -- might be a bug report for Sofie
      clean = TRUE
      plot = FALSE
      verbose = TRUE
      
      clusterRes <- list()
  
  ### File preparation
      message("Training alignment - file (batch) preparation")
      
      file.dat <- model$fsom$data
      file.dat <- as.data.table(file.dat)
      
      train.files <- unique(file.dat$File)
      train.files
      
      model$files
      model$file.nums
      
      if(!all(model$file.nums == train.files)){
        setwd(starting.dir)
        stop("Error -- file mismatch has occured")
      }
      
  ### Loop - each training file as flowframe
      message("Training alignment - file (batch) and metacluster splitting")
      
      setwd(starting.dir)
      dir.create("tmp-train/")
      setwd("tmp-train/")
      getwd()
      
      for(file in train.files){
        # file <- train.files[[1]]
        
        message(" -- running File '", file, "'")
        
        ## Extract file as new flowFrame
        file.dt <- file.dat[file.dat[["File"]] == file,]
        
        metadata <- data.frame(name=dimnames(file.dt)[[2]], desc=dimnames(file.dt)[[2]])
        metadata$range <- apply(apply(file.dt,2,range),2,diff)
        metadata$minRange <- apply(file.dt,2,min)
        metadata$maxRange <- apply(file.dt,2,max)
        
        ## Create flowframe
        file.ff <- new("flowFrame", exprs=as.matrix(file.dt), parameters=Biobase::AnnotatedDataFrame(metadata))
        file.ff
        
        fsom_file <- FlowSOM::NewData(model$fsom, file.ff)
        fsom_file
        
        cellClusterIDs <- model$fsom$metaclustering[fsom_file$map$mapping[,1]]
        
        ## Split individual metaclustersclusters
        
        for (cluster in unique(model$fsom$metaclustering)) {
          # cluster <- unique(model$fsom$metaclustering)[[1]]
          
          if (sum(FlowSOM::GetMetaclusters(fsom_file,
                                           model$fsom$metaclustering) == cluster) > 0) {
            suppressWarnings(
              write.FCS(
                file.ff[cellClusterIDs == cluster,],
                file=file.path(getwd(), #####
                               #outputDir,
                               paste0(gsub("[:/]","_",file),
                                      "_fsom",cluster,".fcs"))))
              # flowCore::write.FCS(
              #   file.ff[cellClusterIDs == cluster,],
              #   file=file.path(getwd(), #####
              #                  #outputDir,
              #                  paste0(gsub("[:/]","_",file),
              #                         "_fsom",cluster,".fcs"))))
          }
        }
      }
  
  ### Learn conversions
      message("Training alignment - learning conversions")
      
      clusterRes <- list()
      
      for (cluster in unique(model$fsom$metaclustering)) {
        if(verbose) message("Processing cluster ",cluster)
        
        normParams_tmp <- c(normParams,
                            list(files = file.path(getwd(),
                                                   #outputDir,
                                                   paste0(gsub("[:/]", "_", train.files),
                                                          "_fsom", cluster, ".fcs")),
                                 labels = as.character(model$files),
                                 channels = align.cols,
                                 transformList = NULL,
                                 verbose = verbose,
                                 plot = plot))
        
        normParams_tmp <- normParams_tmp[unique(names(normParams_tmp))]
        
        if(is.list(normParams[["goal"]])){
          normParams_tmp[["goal"]] <- normParams[["goal"]][[cluster]]
        }
        
        clusterRes[[cluster]] <- do.call(normMethod.train, normParams_tmp)
      }
  
  ### Cleanup
      message("Training alignment - cleanup")
      
      if(clean){
        for(cluster in unique(model$fsom$metaclustering)){
          tmp_files <- file.path(getwd(),
                                 #outputDir,
                                 paste0(gsub("[:/]", "_", train.files),
                                        "_fsom", cluster, ".fcs"))
          
          file.remove(tmp_files[file.exists(tmp_files)])
        }
        setwd(starting.dir)
        unlink("tmp-train", recursive = TRUE)
      }
  
  ### Add conversions to model
  
      model$align.cols <- align.cols
      model$method <- method
      model$conversions <- clusterRes
      
      if(is.null(model$fsom) |
         is.null(model$method) | 
         is.null(model$conversions)){
        stop("The 'conversions' have not been added correctly")
      }
      
      if(mem.ctrl == TRUE){
        gc()
      }
      
}
 
#########################################################################################           
######################################## Wrap up ########################################
#########################################################################################      

  ### Return modified FlowSOM object
      message("Training alignment - training complete")
      return(model)
      
}

