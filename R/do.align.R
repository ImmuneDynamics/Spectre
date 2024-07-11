#' do.align - Function to align multiple batches of a dataset.
#'
#' This function allows you to align multiple batches of a dataset, using either a "Quantile" approach, or using "CytoNorm"
#'
#' @usage do.align(ref.dat, target.dat, batch.col, align.cols, method, goal, 
#' nQ, Qmin, Qmax, write.ref.fcs, write.target.fcs)
#'
#' @param ref.dat NO DEFAULT. If using method = "CytoNorm", this must be a FlowSOM object created by Spectre::do.prep.fsom. If using "Quantiles", then this must be
#' @param target.dat NO DEFAULT. A data.table of data you wish to align
#' @param batch.col NO DEFAULT. Character, column that denotes batches
#' @param align.cols NO DEFAULT. Character, a vector of columns that you wish to align.
#' @param method DEFAULT = "CytoNorm". Character, can be "CytoNorm" or "Quantile".
#' @param goal DEFAULT = "Mean". For method = "CytoNorm". Character, the goal to align to. Can either be "mean", or can be the name of one of the batches (e.g. "Batch1").
#' @param nQ DEFAULT = 101. For method = "CytoNorm". Numeric, the number of quanitles to use.
#' @param Qmin DEFAULT = 0.01. For method = "Quantiles". Numeric, the minimum valus for the lower threshold.
#' @param Qmax DEFAULT = 0.99. For method = "Quantiles". Numeric, the maximum value for the upper threshold.
#' @param write.ref.fcs DEFAULT = TRUE. Logical, do you want a quick copy of the reference sample FCS files.
#' @param write.target.fcs DEFAULT = TRUE. Logical, do you want a quick copy of the ALIGNED target sample FCS files.
#' @param mem.ctrl DEFAULT = TRUE. Runs gc() (garbage collection) after a number of steps to free up memory that hasn't been released quickly enough.
#'
#' @return Returns a data.table with selected columns (align.cols) replaced with aligned data.
#'
#' @author Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#'
#' @references Ashhurst, T. M., et al. (2019). \url{https://www.ncbi.nlm.nih.gov/pubmed/31077106}
#'
#' @usage do.align(ref.dat, target.dat, batch.col, align.cols, method = "CytoNorm", 
#' goal = "mean", nQ = 101, Qmin = 0.01, Qmax = 0.99, write.ref.fcs = TRUE,
#' write.target.fcs = TRUE, mem.ctrl = TRUE)
#'
#' @export do.align

do.align <- function(ref.dat,
                     target.dat,

                     ## Columns
                     #sample.col, # used to divide 'ref' and 'target' data into 'samples'
                     # @param sample.col NO DEFAULT. Character, column that denotes samples

                     batch.col, # Column denoting batches
                     align.cols, # Channels to align

                     ## Method
                     method = "CytoNorm", # Can also be "Quantile"

                     ## CytoNorm settings
                     goal = "mean",
                     nQ = 101,

                     ## Quantile settings
                     Qmin = 0.01,
                     Qmax = 0.99,

                     ## Writing result files
                     write.ref.fcs = TRUE,
                     write.target.fcs = TRUE,
                     mem.ctrl = TRUE
){

  message("The 'do.align' function has been depreciated. Please use 'run.align' instead. Use '? run.align' for more information.")
  
  
#   ###########################################################################################
#   ### Package check
#   ###########################################################################################
# 
#   ### Check that necessary packages are installed
#   if(!is.element('Spectre', installed.packages()[,1])) stop('Spectre is required but not installed')
#   if(!is.element('CytoNorm', installed.packages()[,1])) stop('CytoNorm is required but not installed')
# 
#   ### Require packages
#   require(Spectre)
#   Spectre::package.load()
#   require(data.table)
#   require(CytoNorm)
# 
#   ###########################################################################################
#   ### Test data
#   ###########################################################################################
# 
#   # ref.dat = ref.fsom
#   # target.dat = target.dat
#   #
#   # ## Columns
#   # sample.col <- batch.col
#   # batch.col = batch.col
#   # align.cols = align.cols
#   #
#   # ## CytoNorm settings
#   # method = "CytoNorm"
#   # goal = "mean"
#   # nQ = 101
#   #
#   # ## Writing result files
#   # write.ref.fcs = TRUE
#   # write.target.fcs = TRUE
# 
#   ###########################################################################################
#   ### If using quantile method
#   ###########################################################################################
# 
#   if(method == "Quantile"){
# 
#     stop("do.align -- The Quantile method is not active yet, please use the CytoNorm method")
# 
#     # ## Normalise data (prior to concatenation)
#     # ## --> norm by total dataframe
#     # DataList.norm <- DataList
#     # reset.quant <- function(x){ #Source: https://stackoverflow.com/questions/13339685/how-to-replace-outliers-with-the-5th-and-95th-percentile-values-in-r
#     #   quantiles <- quantile( x, c(.0, .995 ) )
#     #   x[ x < quantiles[1] ] <- quantiles[1]
#     #   x[ x > quantiles[2] ] <- quantiles[2]
#     #   x
#     # } #Sets anything less/greater than 0th/99.5th percentile to be equal to 0th/99.5th percentile
#     # norm.fun <- function(x) {99*(x - min(x, na.rm=TRUE))/(max(x,na.rm=TRUE) -min(x, na.rm=TRUE))+1} #normalises data so min is 1, max is 100
#     # #norm.quant.fun <- function(x) {(x - quantile(x, 0.05, na.rm=TRUE))/(quantile(x, 0.95, na.rm=TRUE) -quantile(x, 0.05, na.rm=TRUE))} #Normalises based on 5th/95th percentile, rather than min/max
#     # is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan)) #"NaN" means "not a number" (happens because all values are the same, so normalising gives back "0/0": NaN)
#     # for (i in AllSampleNos) {
#     #   #DataList.norm[[i]][,1:(length(DataList.norm[[i]])-4)] <- as.data.frame(lapply(DataList.norm[[i]][,1:(length(DataList.norm[[i]])-4)], reset.quant))
#     #   DataList.norm[[i]][,1:(length(DataList.norm[[i]])-4)] <- as.data.frame(lapply(DataList.norm[[i]][,1:(length(DataList.norm[[i]])-4)], norm.fun)) #Ignores last 4 columns which include names
#     #   # DataList.norm[[i]][,1:(length(DataList.norm[[i]])-4)] <- as.data.frame(lapply(DataList.norm[[i]][,1:(length(DataList.norm[[i]])-4)], norm.quant.fun)) #Ignores last 4 columns which include names
#     #   DataList.norm[[i]][is.nan.data.frame(DataList.norm[[i]])] <- 0 #Replaces NaN with 0
#     # }
# 
#   }
# 
#   ###########################################################################################
#   ### Setup for CytoNorm
#   ###########################################################################################
# 
#   if(method == "CytoNorm"){
# 
#     message("do.align -- Step 1/5. Setup started")
# 
#     ### Initial stuff
# 
#     sample.col <- batch.col
#     starting.dir <- getwd()
# 
#     labels <- align.cols
#     channels <- align.cols
#     channel.nums <- as.numeric(as.factor(channels))
# 
#     non.channels <- names(target.dat[, -channel.nums, with=FALSE])
# 
#     ### Merge outputs into Spectre format
#     A <- ref.dat$fsom$FlowSOM$data
#     B <- ref.dat$fsom$FlowSOM$map$mapping[,1]
#     C <- ref.dat$fsom$metaclustering[ref.dat$fsom$FlowSOM$map$mapping[,1]]
# 
#     #nrow(A)
#     #length(B)
#     #length(C)
# 
#     if(nrow(A) != length(B)){
#       stop("Spectre: there was a problem with your 'ref.dat' file")
#     }
# 
#     if(nrow(A) != length(C)){
#       stop("Spectre: there was a problem with your 'ref.dat' file")
#     }
# 
#     if(length(B) != length(C)){
#       stop("Spectre: there was a problem with your 'ref.dat' file")
#     }
# 
#     ### Merge
# 
#     ref.dt <- as.data.table(A)
#     ref.dt <- cbind(A, FlowSOM_cluster = B, FlowSOM_metacluster = C)
#     ref.dt <- as.data.table(ref.dt)
#     #ref.dt
#     #str(ref.dt)
# 
#     test <- ref.dt[,channels, with = FALSE]
#     test2 <- target.dat[,c(channels), with = FALSE]
# 
#     if(!exists("test")){
#       warning("It seems the column names you have specified are not all present in the 'reference' dataset. Please check your entries and try again.")
#     }
# 
#     if(!exists("test2")){
#       warning("It seems the column names you have specified are not all present in the 'target' dataset. Please check your entries and try again.")
#     }
# 
#     rm(test)
#     rm(test2)
# 
#     if(mem.ctrl == TRUE){
#       gc()
#     }
# 
#     #### Cell barcodes for target.dat
#     target.dat[["temp-pre-alignment-barcode"]] <- c(1:nrow(target.dat))
# 
#     ### Setup fsom
#     fsom <- ref.dat$fsom
# 
#     ### Reference samples
#     fsom$files
#     fsom$filenums
#     fsom$batches
# 
#     if(is.unsorted(fsom$filenums) == TRUE){
#       stop("Error -- your files (check $fsom$filenums and $fsom$files) are not in order, and may cause a batch mismatch")
#     }
# 
#     if(all(unique(ref.dt[["File"]]) == unique(fsom$filenums)) == FALSE){
#       stop("Error -- your reference files appear to have been ordered incorrectly in the prep.fsom() function, and may not match their corresponding batches")
#     }
# 
#     if(length(fsom$files) != length(unique(fsom$batches))){
#       stop("Error -- you have a different number of 'files' compared to 'batches'. Please check your prep.fsom() arguments and try again.")
#     }
# 
#     # ref.list <- unique(ref.dt[["File"]]) ##### Causes a mis-match based on ordering
# 
#     ref.list <- fsom$files # ordering based on the fsom.prep file -- disk file order, which is required here
#     ref.list
# 
#     # ref.labels <- vector()
#     #
#     # for(i in ref.list){
#     #   # i <- 1
#     #   ref.labels[i] <- fsom$batches[[i]]
#     # }
#     #
#     # ref.labels <- as.numeric(factor(ref.labels))
# 
#     ref.labels <- fsom$batches # ordering based on the fsom.prep file -- disk file order, which is required here
# 
# 
#     fsom$files
#     ref.list
# 
#     fsom$batches
#     ref.labels
# 
#     ### Target samples
# 
#     target.list <- unique(target.dat[[sample.col]])
#     target.list <- as.character(target.list)
# 
#     # Check if target list already contains "___"
#     # if('___' %in% target.list){
#     #   stop('Error -- it seems your file names contain three underscores in at least one name, which is not compatible with this function (as we use three underscores to add file order information to the file names. Please rename the files and re-run the function.')
#     # }
# 
#     ### ### ### ### ###
#     target.list <- target.list[order(target.list)] ######################## Order it so that it matches the file order when they are made into FCS files
#     target.list
# 
#     starting.target.list <- target.list
# 
#     nr.before <- nrow(target.dat)
# 
#     for(i in c(1:length(target.list))){
#       # i <- 1
# 
#       a <- target.list[[i]]
# 
#       if(nchar(i) == 1){
#         o <- paste0(0, 0, 0)
#       }
# 
#       if(nchar(i) == 2){
#         o <- paste0(0, 0)
#       }
# 
#       if(nchar(i) == 3){
#         o <- paste0(0)
#       }
# 
#       if(nchar(i) == 4){
#         o <- paste0(i)
#       }
# 
#       if(nchar(i) == 5){
#         stop("Error -- the current function configuration is not designed for managing over 9999 batches. How did you even do that? Get in touch at thomas.ashhurst@sydneycytometry.org.au for help.")
#       }
# 
#       temp <- target.list[[i]]
#       temp <- paste0("TARGETFILE", o, i, "___", temp)
#       target.list[[i]] <- temp
#     }
# 
#     d <- data.table("Start" = starting.target.list, "ORDEREDFILE" = target.list)
# 
#     target.dat <- do.embed.columns(target.dat, sample.col, add.dat = d, add.by = "Start", rmv.ext = FALSE)
# 
#     # target.dat[target.dat[[sample.col]] == a, "ORDEREDFILE"] <- temp
#     #  target.dat[target.dat[[sample.col]] == temp,sample.col, with = FALSE]
# 
#     start.sample.col <- sample.col
#     sample.col <- "ORDEREDFILE"
# 
#     nr.after <- nrow(target.dat)
# 
#     if(nr.before != nr.after){
#       stop("Error -- something has gone wrong in renaming the file/batches in the target dataset")
#     }
# 
#     target.labels <- vector()
# 
#     for(i in c(1:length(target.list))){
#       # i <- 1
#       nme <- target.list[i]
#       temp <- target.dat[target.dat[[sample.col]] == nme,]
#       tmp <- temp[[batch.col]][1]
#       target.labels[i] <- as.character(tmp)
#     }
# 
#     target.labels
# 
#     if(all(gsub(".*___", "", target.list) == target.labels) == FALSE){
#       stop("Error -- your list of target batches has mixed batch labels.")
#     }
# 
#     ###########################################################################################
#     ### TRAIN THE MODEL -- Create quantile conversion model (using ref samples)
#     ###########################################################################################
# 
#     message("do.align -- Step 2/5. Training started")
# 
#     ### Some preferences
# 
#     transformList = NULL # in the normal fsom, the default is NULL, but not here -- might be a bug report for Sofie
#     outputDir = "./tmp"
# 
#     normMethod.train = QuantileNorm.train
#     normParams = list(nQ = nQ, goal = goal) ### Can we do more?? Default = 101 Can also be a batch #
#     clean = TRUE
#     plot = FALSE
#     verbose = TRUE
#     #labels <- as.character(dat.training$Batch)
# 
#     clusterRes <- list()
# 
#     ### Convert 'reference' samples into 'temp' FCS files
# 
#     setwd(starting.dir)
#     dir.create("tmp-train")
#     setwd("tmp-train")
# 
#     ref.list
#     ref.labels
#     ref.fsom$fsom$filenums
# 
#     ordr.check <- list()
# 
#     for(i in ref.fsom$fsom$filenums){
# 
#       a <- ref.list[[i]]
#       ordr.check[[i]] <- a
# 
#       if(nchar(i) == 1){
#         o <- paste0(0, 0, 0, i)
#       }
# 
#       if(nchar(i) == 2){
#         o <- paste0(0, 0, i)
#       }
# 
#       if(nchar(i) == 3){
#         o <- paste0(0, i)
#       }
# 
#       if(nchar(i) == 4){
#         o <- paste0(i)
#       }
# 
#       if(nchar(i) == 5){
#         stop("Error -- the current function configuration is not designed for managing over 9999 batches. How did you even do that? Get in touch at thomas.ashhurst@sydneycytometry.org.au for help.")
#       }
# 
#       temp <- ref.dt[ref.dt[["File"]] == i,]
# 
#       write.files(temp,
#                   file.prefix = paste0("File_", o, "_Batch_", a),
#                   write.csv = FALSE,
#                   write.fcs = TRUE)
# 
#       # write.files(dat = ref.dt[ref.dt[["File"]] == i,],
#       #             file.prefix = i,
#       #             write.csv = FALSE,
#       #             write.fcs = TRUE)
#     }
# 
#     train.files <- list.files(getwd(), ".fcs")
# 
#     train.files
#     ref.list
# 
#     if(length(train.files) != length(ref.list)){
#       stop("Error -- somehow the number of files generated is larger than your batches.")
#     }
# 
# 
#     ### Split files by clusters
# 
#     dirCreated = FALSE
#     if (!dir.exists(outputDir)) {
#       dirCreated = dir.create(outputDir)
#     }
# 
#     # Split files by clusters
# 
#     for(file in train.files){
#       if(verbose) message("Splitting ",file)
# 
#       ff <- flowCore::read.FCS(file)
# 
#       if (!is.null(transformList)) {
#         ff <- flowCore::transform(ff,transformList)
#       }
# 
#       # Map the file to the FlowSOM clustering
#       fsom_file <- FlowSOM::NewData(fsom$FlowSOM, ff)
# 
#       # Get the metacluster label for every cell
#       cellClusterIDs <- fsom$metaclustering[fsom_file$map$mapping[,1]]
# 
#       for (cluster in unique(fsom$metaclustering)) {
#         if (sum(FlowSOM::GetMetaclusters(fsom_file,
#                                          fsom$metaclustering) == cluster) > 0) {
#           suppressWarnings(
#             flowCore::write.FCS(
#               ff[cellClusterIDs == cluster,],
#               file=file.path(outputDir,
#                              paste0(gsub("[:/]","_",file),
#                                     "_fsom",cluster,".fcs"))))
#         }
#       }
#     }
# 
# 
#     ### Learn quantiles for each cluster
#     clusterRes <- list()
# 
#     for (cluster in unique(fsom$metaclustering)) {
#       if(verbose) message("Processing cluster ",cluster)
#       if (plot) {
#         grDevices::pdf(file.path(outputDir,
#                                  paste0("CytoNorm_norm_Cluster",
#                                         cluster, ".pdf")),
#                        height = 3*(2*length(train.files)+2),
#                        width = 3*(length(channels)+1))
#       }
# 
#       normParams_tmp <- c(normParams,
#                           list(files = file.path(outputDir,
#                                                  paste0(gsub("[:/]", "_", train.files),
#                                                         "_fsom", cluster, ".fcs")),
#                                labels = as.character(ref.labels),
#                                channels = channels,
#                                transformList = NULL,
#                                verbose = verbose,
#                                plot = plot))
# 
#       normParams_tmp <- normParams_tmp[unique(names(normParams_tmp))]
# 
#       if(is.list(normParams[["goal"]])){
#         normParams_tmp[["goal"]] <- normParams[["goal"]][[cluster]]
#       }
# 
#       clusterRes[[cluster]] <- do.call(normMethod.train,
#                                        normParams_tmp)
# 
#       if (plot) { grDevices::dev.off() }
#     }
# 
# 
#     ### Cleanup
# 
#     if(clean){
#       for(cluster in unique(fsom$metaclustering)){
#         tmp_files <- file.path(outputDir,
#                                paste0(gsub("[:/]", "_", train.files),
#                                       "_fsom", cluster, ".fcs"))
# 
#         file.remove(tmp_files[file.exists(tmp_files)])
#       }
#       if(dirCreated & !plot){
#         unlink(outputDir, recursive=TRUE)
#       }
#     }
# 
#     ### Establish model
# 
#     model <- named.list(fsom, clusterRes) ###### mismatch
#     #model
# 
#     if(mem.ctrl == TRUE){
#       gc()
#     }
# 
#     ###########################################################################################
#     ### APPLY THE MODEL -- Align datasets (apply to all samples)
#     ###########################################################################################
# 
#     message("do.align -- Step 3/5. Alignment started")
# 
#     ### Write 'training' list samples to FCS files
# 
#     getwd()
# 
#     setwd(starting.dir)
#     dir.create("tmp-target")
#     setwd("tmp-target")
# 
#     for(i in target.list){
#       write.files(dat = target.dat[target.dat[[sample.col]] == i,],
#                   file.prefix = i,
#                   write.csv = FALSE,
#                   write.fcs = TRUE)
#     }
# 
#     target.files <- list.files(getwd(), ".fcs")
#     # trg.test <- gsub(".fcs", "", target.files)
# 
#     target.files <- target.files[order(target.files)]
#     # trg.test <- trg.test[order(trg.test)]
# 
#     # if(all(target.list == trg.test) == FALSE){
#     #   stop("Error -- the target files (batches) have become dis-ordered when writing to/reading from disk.")
#     # }
# 
#     if(all(gsub(".*___", "", target.list) == target.labels) == FALSE){
#       stop("Error -- the target files (batches) and corresponding batches don't match.")
#     }
# 
# 
#     # QUICKTEST <- list()
#     #
#     # for(file in target.files) { # Loop to read files into the list
#     #   tempdata <- exprs(flowCore::read.FCS(file, transformation = FALSE))
#     #   tempdata <- tempdata[1:nrow(tempdata),1:ncol(tempdata)]
#     #   tempdata <- as.data.table(tempdata)
#     #   file <- gsub(".fcs", "", file)
#     #   QUICKTEST[[file]] <- tempdata
#     #   rm(tempdata)
#     # }
#     #
#     # QUICKTEST[[1]]
#     # QUICKTEST[[2]]
#     # QUICKTEST[[3]]
#     # QUICKTEST[[4]]
#     # QUICKTEST[[5]]
#     # QUICKTEST[[6]]
#     # QUICKTEST[[7]]
#     # QUICKTEST[[8]]
#     # QUICKTEST[[9]]
#     # QUICKTEST[[10]]
# 
# 
#     ## Align data
#     # CytoNorm.normalize(model = model,
#     #                    files = target.files,
#     #                    labels = target.labels,
#     #                    transformList = NULL,
#     #                    transformList.reverse = NULL,
#     #                    normMethod.normalize = QuantileNorm.normalize,
#     #                    outputDir = "Target_Normalized",
#     #                    prefix = "TaretNorm_",
#     #                    clean = TRUE,
#     #                    verbose = TRUE)
# 
#     ### Run model
# 
#     model = model
#     files = target.files
#     labels = target.labels
#     transformList = NULL
#     transformList.reverse = NULL
#     normMethod.normalize = QuantileNorm.normalize
#     outputDir = "Target_Normalized"
#     prefix = "TargetNorm_"
#     clean = TRUE
#     verbose = TRUE
# 
#     if(is.null(model$fsom) |
#        is.null(model$clusterRes)){
#       stop("The 'model' paramter should be the result of using the
#                trainQuantiles function.")
#     }
# 
#     if(length(labels) != length(files)){
#       stop("Input parameters 'labels' and 'files' should have the same length")
#     }
# 
#     # Create output directory
#     if(!dir.exists(outputDir)){
#       dir.create(outputDir)
#     }
# 
#     fsom <- model$fsom
#     clusterRes <- model$clusterRes
# 
#     ### Split files by clusters
# 
#     cellClusterIDs <- list()
#     meta <- list()
#     cluster_files <- list()
# 
# 
#     for(file in files){
#       # file <- "01_Air_01.fcs"
#       if(verbose) message("Splitting ",file)
#       ff <- flowCore::read.FCS(file)
# 
#       if(!is.null(transformList)){
#         ff <- flowCore::transform(ff, transformList)
#         # meta[[file]] <- list()
#         # meta[[file]][["description_original"]] <- ff@description
#         # meta[[file]][["parameters_original"]] <- ff@parameters
#       }
# 
#       fsom_file <- FlowSOM::NewData(fsom$FlowSOM,ff)
#       cellClusterIDs[[file]] <- fsom$metaclustering[fsom_file$map$mapping[,1]]
# 
#       for(cluster in unique(fsom$metaclustering)){
#         if (sum(FlowSOM::GetMetaclusters(fsom_file,
#                                          fsom$metaclustering) == cluster) > 0) {
#           f <- file.path(outputDir,
#                          paste0(gsub("[:/]","_",file),
#                                 "_fsom", cluster, ".fcs"))
#           suppressWarnings(
#             flowCore::write.FCS(ff[cellClusterIDs[[file]] == cluster],
#                                 file = f)
#           )
#         }
#       }
#       rm(ff)
#     }
# 
#     ### Apply normalization on each cluster
#     for(cluster in unique(fsom$metaclustering)){
#       if(verbose) message("Processing cluster ",cluster)
#       files_tmp <- file.path(outputDir,
#                              paste0(gsub("[:/]",
#                                          "_",
#                                          files),
#                                     "_fsom",
#                                     cluster,
#                                     ".fcs"))
# 
#       ### ### ### ### ###
# 
#       labels_tmp <- labels[file.exists(files_tmp)]
#       files_tmp <- files_tmp[file.exists(files_tmp)] #### different order
#       normMethod.normalize(model = clusterRes[[cluster]],
#                            files = files_tmp,
#                            labels = labels_tmp,
#                            outputDir = file.path(outputDir),
#                            prefix = "Norm_",
#                            transformList = NULL,
#                            transformList.reverse = NULL,
#                            removeOriginal = TRUE,
#                            verbose = verbose)
#     }
# 
#     ###
# 
#     strt <- getwd()
#     setwd("Target_Normalized/")
# 
#     #list.files(getwd(), ".fcs")
#     #all.clust <- read.files(getwd(), file.type = ".fcs", do.embed.file.names = FALSE)
# 
#     all.clust <- list()
#     fl.nms <- list.files(path=getwd(), pattern = ".fcs")
#     fl.nms
# 
#     fl.nrws <- list()
#     fl.ncls <- list()
# 
#     for(file in fl.nms) { # Loop to read files into the list
#       tempdata <- exprs(flowCore::read.FCS(file, transformation = FALSE))
# 
#       if(nrow(tempdata) == 1){
#         tempdata <- tempdata[1:nrow(tempdata),1:ncol(tempdata)]
#         tempdata <- as.data.table(t(tempdata))
#         tempdata
#       }
# 
#       if(nrow(tempdata) > 1){
#         tempdata <- tempdata[1:nrow(tempdata),1:ncol(tempdata)]
#         tempdata <- as.data.table(tempdata)
#       }
# 
#       file <- gsub(".fcs", "", file)
#       all.clust[[file]] <- tempdata
#       fl.ncls[[file]] <- ncol(tempdata)
#       fl.nrws[[file]] <- nrow(tempdata)
#       rm(tempdata)
#     }
# 
#     if(length(unique(unlist(fl.ncls))) > 1){
#       stop("Error -- when combining the 'file-cluster' FCS files back into 'file' datasets, differenting number of columns were found in one 'file-cluster' FCS file.")
#     }
# 
#     # all.clust[[18]]
# 
#     names(all.clust)
#     names(all.clust[[1]])
# 
#     for(i in c(1:length(all.clust))){
#       a <- names(all.clust)[[i]]
#       all.clust[[i]]$NORM_TARGET_FILENAME <- a
#       all.clust[[i]]$NORM_TARGET_FILENUM <- i
#     }
# 
#     names(all.clust)
#     names(all.clust[[1]])
# 
#     setwd(strt)
#     getwd()
# 
# 
#     ### Combine clusters into one final fcs file
#     for(file in files){
#       if(verbose) message("Rebuilding ",file)
# 
#       ff <- flowCore::read.FCS(file)
# 
#       for(cluster in unique(fsom$metaclustering)){
#         file_name <- file.path(outputDir,
#                                paste0("Norm_",gsub("[:/]","_",file),
#                                       "_fsom",cluster,".fcs"))
#         if (file.exists(file_name)) {
#           ff_subset <- flowCore::read.FCS(file_name)
#           flowCore::exprs(ff)[cellClusterIDs[[file]] == cluster,] <- flowCore::exprs(ff_subset)
#         }
#       }
# 
#       if(!is.null(transformList.reverse)){
#         ff <- flowCore::transform(ff, transformList.reverse)
#         # ff@description <- meta[[file]][["description_original"]]
#         # ff@parameters <- meta[[file]][["parameters_original"]]
#       }
# 
# 
#       # Adapt to real min and max because this gets strange values otherwise
#       ff@parameters@data[,"minRange"] <- apply(ff@exprs, 2, min)
#       ff@parameters@data[,"maxRange"] <- apply(ff@exprs, 2, max)
#       ff@parameters@data[,"range"] <- ff@parameters@data[,"maxRange"] -
#         ff@parameters@data[,"minRange"]
# 
#       if(clean){
#         file.remove(file.path(outputDir,
#                               paste0("Norm_",gsub("[:/]","_",file),
#                                      "_fsom",unique(fsom$metaclustering),".fcs")))
#       }
# 
#       suppressWarnings(flowCore::write.FCS(ff,
#                                            file=file.path(outputDir,
#                                                           paste0(prefix,gsub(".*/","",file)))))
#     }
# 
#     if(mem.ctrl == TRUE){
#       gc()
#     }
# 
#     ###########################################################################################
#     ### Read FCS files back in and create merged DT
#     ###########################################################################################
# 
#     message("do.align -- Step 4/5. Re-import aligned files and creating 'modified' aligned data.table")
# 
#     ###
# 
#     setwd(starting.dir)
#     setwd("tmp-target/")
#     setwd("Target_Normalized")
# 
#     # all.clust.strt <- all.clust
#     # all.clust <- all.clust.strt
#     names(all.clust[[1]])
# 
#     sve <- all.clust
# 
#     for(a in names(all.clust)){
#       # a <- names(all.clust)[[1]]
#       b <- a
#       colnames(all.clust[[a]])[colnames(all.clust[[a]]) %in% c("NORM_TARGET_FILENAME")] <- c("RAWFILENAME")
#       a <- gsub('\\_fsom.*', '', a, perl=TRUE)
#       all.clust[[b]][["NORM_TARGET_FILENAME"]] <- a
#     }
# 
#     #all.clust
#     #names(all.clust)
#     #names(all.clust[[1]])
#     #names(all.clust[[70]])
# 
#     if(mem.ctrl == TRUE){
#       gc()
#     }
# 
#     nme.list <- list()
# 
#     for(i in c(1:length(all.clust))){
#       temp <- all.clust[[i]]
#       nme.list[[i]] <- ncol(temp)
#       rm(temp)
#     }
# 
#     if(length(unique(unlist(nme.list))) < 1){
#       stop("Error -- differing number of column names in each file")
#     }
# 
#     aligned.dat <- data.table::rbindlist(all.clust, fill = TRUE)
# 
#     if(mem.ctrl == TRUE){
#       gc()
#     }
# 
#     aligned.dat[, channels, with=FALSE]
# 
#     ### Replace target.dat with aligned.dat -- just the channels aligned
#     mod.dat <- target.dat
# 
#     aligned.dat <- aligned.dat[order(aligned.dat[['temp-pre-alignment-barcode']]),]
# 
#     if(is.unsorted(aligned.dat[['temp-pre-alignment-barcode']])){
#       stop("Error -- pre-alignment cellular barcodes not in correct order")
#     }
# 
# 
#     if(mem.ctrl == TRUE){
#       gc()
#     }
# 
#     for(a in channels){
#       # mod.dat[[a]] <- aligned.dat[[a]]
#       # names(mod.dat)[names(mod.dat) == a] <- paste0(a, "_align")
# 
#       ## Test function for 'adding' data
#       i <- paste0(a, "_aligned")
#       mod.dat[[i]] <- aligned.dat[[a]]
#     }
# 
#     if(nrow(mod.dat) != nrow(aligned.dat)){
#       stop("Error -- differing rows in 'mod.dat' and 'aligned.dat'")
#     }
# 
#     mod.dat[["temp-post-alignment-barcode"]] <- aligned.dat[['temp-pre-alignment-barcode']]
# 
#     if(mem.ctrl == TRUE){
#       gc()
#     }
# 
#     names(mod.dat)
# 
#     ###########################################################################################
#     ### Return mod.dat and cleanup
#     ###########################################################################################
# 
#     message("do.align -- Step 5/5. Cleanup, writing FCS files, and returning aligned data")
# 
#     ## WD and Directory stuff
# 
#     setwd(starting.dir)
#     list.dirs(getwd(), full.names = FALSE, recursive = FALSE)
# 
#     unlink("tmp-train", recursive = TRUE)
#     unlink("tmp-target", recursive = TRUE)
# 
#     ### Replacing batch names with originals
# 
#     mod.target.list <- target.list
#     target.list <- starting.target.list
# 
# 
#     mod.sample.col <- sample.col
#     sample.col <- start.sample.col
# 
# 
#     target.dat$ORDEREDFILE <- NULL
#     mod.dat$ORDEREDFILE <- NULL
# 
#     # final.target.list <- gsub(".*___", "", starting.target.list)
#     #
#     # if(all(starting.target.list == final.target.list) == FALSE){
#     #   stop("Error -- somehow the list of target files (batches) has become mixed up")
#     # }
#     #
#     # if(all(unique(target.dat$ORDEREDFILE) == unique(mod.dat$ORDEREDFILE)) == FALSE){
#     #   stop("Error")
#     # }
#     #
#     # #
#     #
#     # for(i in c(1:length(target.list))){
#     #   a <- target.list[[i]]
#     #   b <- starting.target.list[[i]]
#     #
#     #   target.dat$ORDEREDFILE <- NULL
#     #   mod.dat$ORDEREDFILE <- NULL
#     # }
#     #
#     # target.list <- final.target.list
# 
# 
#     ### Write reference FCS files
# 
#     if(write.ref.fcs == TRUE){
#       setwd(starting.dir)
#       dir.create("CytoNorm_ref_output")
#       setwd("CytoNorm_ref_output")
# 
#       ref.fsom$fsom$filenums
# 
#       for(i in ref.fsom$fsom$filenums){
#         a <- ref.fsom$fsom$files[[i]]
#         write.files(dat = ref.dt[ref.dt[["File"]] == i,],
#                     file.prefix = paste0("Raw_RefFile_", i, "_", a),
#                     write.csv = FALSE,
#                     write.fcs = TRUE)
#       }
#       setwd(starting.dir)
#     }
# 
#     if(mem.ctrl == TRUE){
#       gc()
#     }
# 
#     ### Write target FCS files
# 
#     if(write.target.fcs == TRUE){
#       setwd(starting.dir)
#       dir.create("CytoNorm_target_output")
#       setwd("CytoNorm_target_output")
# 
#       for(i in target.list){
#         write.files(dat = target.dat[target.dat[[sample.col]] == i,],
#                     file.prefix = paste0("Raw_", i),
#                     write.csv = FALSE,
#                     write.fcs = TRUE)
#       }
# 
#       for(i in unique(mod.dat[[sample.col]])){
#         write.files(dat = mod.dat[mod.dat[[sample.col]] == i,],
#                     file.prefix = paste0("Aligned_", i),
#                     write.csv = FALSE,
#                     write.fcs = TRUE)
#       }
#       setwd(starting.dir)
#     }
# 
#     if(mem.ctrl == TRUE){
#       gc()
#     }
# 
#     ### Return
# 
#     return(mod.dat)
#     message("do.align complete")
# 
#   }
}
