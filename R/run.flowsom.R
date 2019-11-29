#' run.flowsom - ...
#'
#' @usage run.flowsom(x, ...)
#'
#' @param x data.frame. Input sample. No default.
#' @param meta.k Numeric. Number of clusters to create. No default.
#' @param clustering.cols Vector of column names to use for clustering. No default.
#' @param clust.seed Numeric. Clustering seed for reproducibility. No default.
#' @param meta.seed Numeric. Metaclustering seed for reproducibility. No default.
#' @param clust.name Character. Name of the resulting 'cluster' parameter. Defaults to "FlowSOM_cluster".
#' @param meta.clust.name Character. Name of the resulting 'metacluster' parameter. Defaults to "FlowSOM_metacluster".
#'
#' This function runs FlowSOM on a dataframe with cells (rows) vs markers (columns), and returns 'res' with result columns
#'
#' @export

run.flowsom <- function(data,
                        xdim = 10,
                        ydim = 10,
                        meta.clust = 20,
                        
                        clust.seed = 42,
                        meta.seed = 42,
                        
                        clustering.cols,
                        clust.name = "FlowSOM_cluster",
                        meta.clust.name = "FlowSOM_metacluster"){
  
  ## Test data
  
      as.matrix(names(demo.start))
  
      data <- demo.start
      xdim <- 10
      ydim <- 10
      meta.clust <- 20
      clust.seed <- 42
      meta.seed <- 42
      clustering.cols <- c(5,6,8,13,17,18,29,21,23,26)
      clust.name <- "FlowSOM_cluster",
      meta.clust.name <- "FlowSOM_metacluster")
  
  ## Data conversion from dataframe or datatable into flowFrame
      metadata <- data.frame(name=dimnames(x)[[2]], desc=paste('column',dimnames(x)[[2]],'from dataset'))
      x.ff <- new("flowFrame",
                  exprs=as.matrix(x), # in order to create a flow frame, data needs to be read as matrix
                  parameters=AnnotatedDataFrame(metadata))
      
      head(exprs(x.ff))
      x_FlowSOM <- x.ff
  
  ## Run FlowSOM
      fSOM <- FlowSOM(fileName,
                      ## Input options:
                      #compensate = TRUE, 
                      #transform = TRUE, 
                      #toTransform=c(8:18),
                      #scale = TRUE,
                      ## SOM options:
                      colsToUse = c(9,12,14:18),# can be column numbers or characters
                      xdim = 7, # 
                      ydim = 7,
                      # Metaclustering options:
                      nClus = 10)
  
  
  
  ## Save output
  
  
}                        
                        
                        
                        
                        
                        
                        meta.k,
                        clustering.cols, # names of columns to cluster
                        clust.seed,
                        meta.seed,
                        clust.name = "FlowSOM_cluster",
                        meta.clust.name = "FlowSOM_metacluster"){
  
  






run.flowsom <- function(x,
                        meta.k,
                        clustering.cols, # names of columns to cluster
                        clust.seed,
                        meta.seed,
                        clust.name = "FlowSOM_cluster",
                        meta.clust.name = "FlowSOM_metacluster"){

  #### TEST VALUES
      #x <- cell.dat
      #meta.k <- 40
      #clustering.cols <- ClusteringCols
      #clust.seed <- 42
      #meta.seed <- 42
      #clust.name <- "FlowSOM_cluster"
      #meta.clust.name <- "FlowSOM_metacluster"

  ##
      #head(x)
      #dimnames(x)[[2]]

  ## Create FCS file metadata - column names with descriptions
  metadata <- data.frame(name=dimnames(x)[[2]], desc=paste('column',dimnames(x)[[2]],'from dataset'))

  ## Create FCS file metadata - ranges, min, and max settings -- by default, they are commented out (adjust ranges manually in FlowJo)
      #metadata$range <- apply(apply(data,2,range),2,diff) # throws an error because of word entry -- hopefully is ok
      #metadata$minRange <- apply(data,2,min)
      #metadata$maxRange <- apply(data,2,max)

  ## Create flowframe with data
  x.ff <- new("flowFrame",
                 exprs=as.matrix(x), # in order to create a flow frame, data needs to be read as matrix
                 parameters=AnnotatedDataFrame(metadata))

  head(exprs(x.ff))

  x_FlowSOM <- x.ff

  # choose markers for FlowSOM analysis
  FlowSOM_cols <- clustering.cols

  ### 4.3. - Run FlowSOM

  ## set seed for reproducibility
  set.seed(clust.seed)

  ## run FlowSOM (initial steps prior to meta-clustering)
  FlowSOM_out <- FlowSOM::ReadInput(x_FlowSOM, transform = FALSE, scale = FALSE)
  FlowSOM_out <- FlowSOM::BuildSOM(FlowSOM_out, colsToUse = FlowSOM_cols)
  FlowSOM_out <- FlowSOM::BuildMST(FlowSOM_out)

  ### some warnings will be returned because of the 'SampleName' and 'GroupName' entries

  ## Optional visualization

  #FlowSOM::PlotStars(FlowSOM_out) # won't plot names if arial font problems still exists, # if 'sample name' is in there, can't plot

  #set.seed(42)
  #FlowSOM::PlotStars(FlowSOM_out,view="tSNE")

  #print(colnames(FlowSOM_out$map$medianValues))
  #FlowSOM::PlotMarker(FlowSOM_out,"BUV395.CD11b")
  #FlowSOM::PlotNumbers(UpdateNodeSize(FlowSOM_out,reset=TRUE))

  ## extract cluster labels (pre meta-clustering) from output object
  labels_pre <- FlowSOM_out$map$mapping[, 1]
  labels_pre
  length(labels_pre)
  nrow(cell.dat)

  flowsom.res.original <- labels_pre

  ## run meta-clustering
  FlowSOM_out_meta <- FlowSOM::metaClustering_consensus(FlowSOM_out$map$codes, k = meta.k, seed = meta.seed)

  # note: In the PREVIOUS version of FlowSOM, the meta-clustering function
  # FlowSOM::metaClustering_consensus() does not pass along the seed argument
  # correctly, so results are not reproducible. We use the internal function
  # ConsensusClusterPlus::ConsensusClusterPlus() to get around this.

  # seed <- 1234
  # out <- ConsensusClusterPlus::ConsensusClusterPlus(t(out$map$codes), maxK = FlowSOM_kvalue, seed = seed)
  # out <- out[[FlowSOM_kvalue]]$consensusClass

  # However, this
  # HAS BEEN fixed in the next update of FlowSOM (version 1.5); then the following
  # (simpler) code can be used instead:

  ## Optional visualisation
  # FlowSOM::PlotMarker(FlowSOM_out,"BUV395.CD11b", backgroundValues = as.factor(FlowSOM_out_meta))


  ## extract META (?) cluster labels from output object
  labels <- FlowSOM_out_meta[labels_pre]

  ## summary of cluster sizes and number of clusters
  table(labels)
  length(table(labels))

  ## save ORIGINAL cluster labels
  flowsom.res.original <- data.frame("labels_pre" = labels_pre)
  colnames(flowsom.res.original)[grepl('labels_pre',colnames(flowsom.res.original))] <- clust.name

  dim(x)
  dim(flowsom.res.original)
  head(flowsom.res.original)

  ## save META cluster labels
  flowsom.res.meta <- data.frame("labels" = labels)
  colnames(flowsom.res.meta)[grepl('labels',colnames(flowsom.res.meta))] <- meta.clust.name

  dim(x)
  dim(flowsom.res.meta)
  head(flowsom.res.meta)

  assign("flowsom.res.original", flowsom.res.original, envir = globalenv())
  assign("flowsom.res.meta", flowsom.res.meta, envir = globalenv())

}
