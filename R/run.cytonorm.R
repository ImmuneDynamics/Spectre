#' run.cytonorm - Run the alignment model on a target data.table
#'
#' This function allows you to prepare reference data ahead of performing batch alignment. 
#' 
#' @usage run.align()
#'
#' @param dat NO DEFAULT. A data.table consisting of the 'refernece' data you will use to train the alignment algorithm
#' @param model NO DEFAULT. A batch alignment conversion model object created by the prep.align() and train.align() functions.
#' @param batch.col NO DEFAULT. Column name of the data.table that contains batch labels
#' @param append.name DEFAULT = "_aligned". This will be appended to the column names containing the aligned data
#' @param dir DEFAULT = getwd(). Sets the working directory to operate from. Because this function involves some reading/writing of files, it's best to set this to somewhere static in case the active working directory moves to a subfolder, and then doesn't return because the function runs into an error.
#' @param mem.ctrl DEFAULT = TRUE. Allows the function to clear held memory on occasion.

#' @return Returns a data.table with aligned data added in new columns.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references Ashhurst, T. M., et al. (2019). \url{https://www.ncbi.nlm.nih.gov/pubmed/31077106}
#'
#' @examples
#' cell.dat <- run.cytonorm()
#'
#' @import data.table
#'
#' @export

run.cytonorm <- function(dat,
                         model,
                         batch.col,
                         append.name = "_aligned",
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
  require(flowCore)
  require(Biobase)
  
  ### Checks
  
  if(is.null(model$fsom) |
     is.null(model$dt) |
     is.null(model$cellular.cols) |
     is.null(model$cluster.cols) |
     is.null(model$files) |
     is.null(model$file.nums)){
    stop("The model object has not been setup correctly")
  }
  
  if(is.null(model$fsom) |
     is.null(model$method) | 
     is.null(model$conversions)){
    stop("The 'conversions' have not been added to the model correctly")
  }
  
  
  ### What type of conversion model
  
  method <- model$method
  
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
    
    ### Preferences
    
    transformList = NULL
    transformList.reverse = NULL
    normMethod.normalize = QuantileNorm.normalize
    #outputDir = "Target_Normalized"
    prefix = "TargetNorm_"
    clean = TRUE
    verbose = TRUE
    
    ### Directories
    
    setwd(dir)
    starting.dir<- getwd()
    message("Working directory is '", starting.dir, "'")
    
    ### Setup data
    
    starting.dat <- dat
    
    starting.dat$ALGN_BARCODE <- c(1:nrow(starting.dat))
    starting.dat$ALGN_BARCODE <- as.numeric(starting.dat$ALGN_BARCODE)
    
    dat <- starting.dat
    
    cellular.cols <- model$cellular.cols
    cluster.cols <- model$cluster.cols
    align.cols <- model$align.cols
    
    clust.align <- unique(cluster.cols, align.cols)
    all.cols <- unique(cellular.cols, clust.align)
    
    all.cols <- c(all.cols, 'ALGN_BARCODE')
    dat <- dat[,c(batch.col, all.cols), with = FALSE]
    
    ### Split files 
    
    message("Alignment -- splitting files (batches)")
    
    setwd(starting.dir)
    
    unlink("tmp-train", recursive = TRUE)
    dir.create("tmp-target")
    setwd("tmp-target")
    unlink("EachCluster") # in case there is anything leftover
    
    dat.list <- unique(dat[[batch.col]])
    
    for(i in c(1:length(dat.list))){
      # i <- 1
      
      a <- dat.list[[i]]
      temp <- dat[dat[[batch.col]] == a,]
      temp <- temp[,c(all.cols),with = FALSE]
      
      write.files(temp,
                  file.prefix = a,
                  write.csv = FALSE,
                  write.fcs = TRUE)
      
      rm(i)
      rm(a)
    }
    
    rm(dat.list)
    
    target.files <- list.files(getwd(), ".fcs")
    target.file.nums <- c(1:length(target.files))
    target.labels <- gsub(".fcs", "", target.files)
    
    ### Process each file (batch)
    
    message("Alignment -- processing files (batches)")
    
    cellClusterIDs <- list()
    meta <- list()
    cluster_files <- list()
    
    for(file in target.files){
      # file <- target.files[[1]]
      if(verbose) message("-- Splitting ",file)
      
      setwd(dir)
      dir.create("tmp-target")
      setwd("tmp-target")
      
      ff <- flowCore::read.FCS(file)
      
      setwd(dir)
      dir.create("tmp-target", showWarnings = FALSE)
      setwd("tmp-target")
      dir.create("EachCluster", , showWarnings = FALSE)
      setwd("EachCluster")
      
      fsom_file <- FlowSOM::NewData(model$fsom$FlowSOM,ff)
      
      cellClusterIDs[[file]] <- model$fsom$metaclustering[fsom_file$map$mapping[,1]]
      
      for(cluster in unique(model$fsom$metaclustering)){
        if (sum(FlowSOM::GetMetaclusters(fsom_file,
                                         model$fsom$metaclustering) == cluster) > 0) {
          f <- file.path(getwd(),
                         #outputDir,
                         paste0(gsub("[:/]","_",file),
                                "_fsom", cluster, ".fcs"))
          suppressWarnings(
            flowCore::write.FCS(ff[cellClusterIDs[[file]] == cluster],
                                file = f)
          )
        }
      }
      rm(ff)
    }
    
    ### Apply normalization on each cluster (modifies the cluster FCS file on disk)
    message("Alignment -- aligning data")
    
    for(cluster in unique(model$fsom$metaclustering)){
      
      if(verbose) message("-- Processing cluster ",cluster)
      files_tmp <- file.path(getwd(),
                             #outputDir,
                             paste0(gsub("[:/]", "_", target.files), "_fsom", cluster, ".fcs")
      )
      
      labels_tmp <- target.labels[file.exists(files_tmp)]
      files_tmp <- files_tmp[file.exists(files_tmp)] #### different order
      normMethod.normalize(model = model$conversions[[cluster]],
                           files = files_tmp,
                           labels = labels_tmp,
                           outputDir = file.path(getwd()),
                           prefix = "Norm_",
                           transformList = NULL,
                           transformList.reverse = NULL,
                           removeOriginal = TRUE,
                           verbose = verbose)
    }
    
    ### Combine clusters into one final fcs file
    message("Alignment -- merging aligned data for each metacluster")
    
    cluster.files <- list.files(getwd(), ".fcs")
    norm.dat <- list()
    
    for(i in cluster.files){
      # i <- cluster.files[[1]]
      temp <- flowCore::read.FCS(i)
      dt <- as.data.table(temp@exprs)
      
      mc <- sub(".*.fcs_", "", i)
      mc <- sub(".fcs", "", mc)
      mc <- sub('fsom', "", mc)
      
      nme <- sub(".fcs_.*", "", i)
      nme <- sub("Norm_", "", nme)
      
      align.mc.nme <- paste0("Alignment_MC", append.name)
      
      dt[[batch.col]] <- nme
      dt[[align.mc.nme]] <- mc
      norm.dat[[i]] <- dt
    }
    
    ## Remove any previous 'tmp' folder
    setwd(starting.dir)
    unlink("tmp-target", recursive = TRUE)
    
    names(norm.dat)
    
    norm.dt <- rbindlist(norm.dat, fill = TRUE)
    norm.dt <- setorderv(norm.dt, 'ALGN_BARCODE')
    
    # Before
    starting.dat
    
    # After
    norm.dt
    
    if(!all(starting.dat$ALGN_BARCODE == norm.dt$ALGN_BARCODE)){
      setwd(dir)
      message("WARNING -- the pre- and post-alignment row barcodes do not align")
    }
    
    if(!all(starting.dat[[batch.col]] == norm.dt[[batch.col]])){
      setwd(dir)
      message("WARNING -- the pre- and post-alignment batch labels")
    }
    
    align.mc <- norm.dt[[align.mc.nme]]
    align.mc <- as.data.table(align.mc)
    names(align.mc) <- align.mc.nme
    align.mc[[align.mc.nme]] <- as.numeric(align.mc[[align.mc.nme]])
    
    norm.dt <- norm.dt[,c(align.cols),with = FALSE]
    names(norm.dt) <- paste0(names(norm.dt), append.name)
    
    ### Creating combined data.table
    message("Alignment -- creating final data.table")
    
    starting.dat$ALGN_BARCODE <- NULL
    final.dat <- cbind(starting.dat, norm.dt, align.mc)
    
  }
  
  #########################################################################################           
  ######################################## WRAP UP ######################################## 
  #########################################################################################           
  
  ### Return
  message("Alignment -- complete")
  return(final.dat)
  
}


