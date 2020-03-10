#' run.flowsom - ...
#'
#' @usage run.flowsom(dat, clustering.cols, meta.k, xdim, ydim, clust.seed, meta.seed, clust.name, meta.clust.name, ...)
#'
#' @param dat NO DEFAULT. data.frame. Input sample.
#' @param clustering.cols NO DEFAULT. Vector of column names to use for clustering. It is possible to use a vector of column numbers here but this is not recommended.
#' @param meta.k DEFAULT = 20. Numeric. Number of clusters to create. If set to zero (0), no metaclusters will be created.
#' @param xdim DEFAULT = 10. Numeric. Number of first level clusters across the x-axis. xdim x ydim = total number of first level FlowSOM clusters.
#' @param ydim DEFAULT = 10. Numeric. Number of first level clusters across the y-axis. xdim x ydim = total number of first level FlowSOM clusters.
#' @param clust.seed DEFAULT = 42 Numeric. Clustering seed for reproducibility.
#' @param meta.seed DEFAULT = 42 Numeric. Metaclustering seed for reproducibility.
#' @param clust.name DEFAULT = "FlowSOM_cluster". Character. Name of the resulting 'cluster' parameter.
#' @param meta.clust.name DEFAULT = "FlowSOM_metacluster". Character. Name of the resulting 'metacluster' parameter.
#'
#' This function runs FlowSOM on a dataframe with cells (rows) vs markers (columns), and returns 'res' with result columns
#'
#' @export

run.flowsom <- function(dat,
                        clustering.cols, # names of columns to cluster
                        meta.k = 20,
                        xdim = 10,
                        ydim = 10,
                        clust.seed = 42,
                        meta.seed = 42,
                        clust.name = "FlowSOM_cluster",
                        meta.clust.name = "FlowSOM_metacluster"){

  #### TEST VALUES
      # dat <- demo.start
      #
      # ##
      # ColumnNames <- as.matrix(unname(colnames(dat))) # assign reporter and marker names (column names) to 'ColumnNames'
      # ColumnNames
      # ClusteringColNos <- c(5,6,8,9,11,13,17:19,21:29,32)
      # ClusteringCols <- ColumnNames[ClusteringColNos]
      #
      # clustering.cols <- ClusteringCols
      #
      # xdim <- 10
      # ydim <- 10
      # meta.k <- 40
      #
      # clust.seed <- 42
      # meta.seed <- 42
      # clust.name <- "FlowSOM_cluster"
      # meta.clust.name <- "FlowSOM_metacluster"

  ##
      #head(dat)
      #dimnames(dat)[[2]]
  
  ## Check that necessary packages are installed
  if(!is.element('Spectre', installed.packages()[,1])) stop('Spectre is required but not installed')
  if(!is.element('flowCore', installed.packages()[,1])) stop('flowCore is required but not installed')
  if(!is.element('Biobase', installed.packages()[,1])) stop('Biobase is required but not installed')
  if(!is.element('FlowSOM', installed.packages()[,1])) stop('FlowSOM is required but not installed')
  if(!is.element('data.table', installed.packages()[,1])) stop('data.table is required but not installed')

  ## Require packages
  require(Spectre)
  require(flowCore)
  require(Biobase)
  require(FlowSOM)
  require(data.table)
  
  ## Remove non-numeric
      head
      dat.start <- dat

      nums <- unlist(lapply(dat, is.numeric))
      dat <- as.data.frame(dat)[ , nums]
      dat[clustering.cols]

  ## Create FCS file metadata - column names with descriptions
  metadata <- data.frame(name=dimnames(dat)[[2]], desc=paste('column',dimnames(dat)[[2]],'from dataset'))

  ## Create flowframe with data
  dat.ff <- new("flowFrame",
                 exprs=as.matrix(dat), # in order to create a flow frame, data needs to be read as matrix
                 parameters=Biobase::AnnotatedDataFrame(metadata))

  head(flowCore::exprs(dat.ff))

  dat_FlowSOM <- dat.ff

  # choose markers for FlowSOM analysis
  FlowSOM_cols <- clustering.cols

  ### 4.3. - Run FlowSOM

  ## set seed for reproducibility
  set.seed(clust.seed)

  ## run FlowSOM (initial steps prior to meta-clustering)
  FlowSOM_out <- FlowSOM::ReadInput(dat_FlowSOM, transform = FALSE, scale = FALSE)

  FlowSOM_out <- FlowSOM::BuildSOM(FlowSOM_out,
                                   colsToUse = FlowSOM_cols,
                                   xdim = xdim,
                                   ydim = ydim)

  FlowSOM_out <- FlowSOM::BuildMST(FlowSOM_out)

  ## extract cluster labels (pre meta-clustering) from output object
  labels_pre <- FlowSOM_out$map$mapping[, 1]
  labels_pre
  length(labels_pre)
  nrow(dat)

  flowsom.res.original <- labels_pre

  if (meta.k != 0) {
    ## run meta-clustering
    FlowSOM_out_meta <- FlowSOM::metaClustering_consensus(FlowSOM_out$map$codes, k = meta.k, seed = meta.seed)

    ## extract META (?) cluster labels from output object
    labels <- FlowSOM_out_meta[labels_pre]

    ## summary of cluster sizes and number of clusters
    table(labels)
    length(table(labels))

    ## save META cluster labels
    flowsom.res.meta <- data.frame("labels" = labels)
    colnames(flowsom.res.meta)[grepl('labels',colnames(flowsom.res.meta))] <- meta.clust.name

    dim(dat)
    dim(flowsom.res.meta)
    head(flowsom.res.meta)

    assign("flowsom.res.meta", flowsom.res.meta, envir = globalenv())

    dat.start <- cbind(dat.start, flowsom.res.meta)       # Add results to dat
  }

  ## save ORIGINAL cluster labels
  flowsom.res.original <- data.frame("labels_pre" = labels_pre)
  colnames(flowsom.res.original)[grepl('labels_pre',colnames(flowsom.res.original))] <- clust.name

  dim(dat)
  dim(flowsom.res.original)
  head(flowsom.res.original)

  assign("flowsom.res.original", flowsom.res.original, envir = globalenv())

  dat.start <- cbind(dat.start, flowsom.res.original)   # Add results to dat

  dat.start <- data.table::as.data.table(dat.start) # Make dat a data.table for future manipulation

}
