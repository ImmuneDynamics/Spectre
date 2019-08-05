# Spectre::flowsom


run.flowsom <- function(x,
                        meta.k,
                        clustering.cols, # names of columns to cluster
                        clust.seed,
                        meta.seed,
                        clust.name,
                        meta.clust.name){

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




