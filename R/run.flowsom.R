#' Run the FlowSOM algorithm
#'
#' Method to run the FlowSOM clustering algorithm.
#' This function runs FlowSOM on a data.table with cells (rows) vs markers (columns) with new columns for FlowSOM clusters and metaclusters.
#' Output data will be "flowsom.res.original" (for clusters) and "flowsom.res.meta" (for metaclusters).
#' Uses the R packages "FlowSOM" for clustering, "flowCore" for handling .fcs files, "Biobase" for creating a flow frame, "data.table" for handling data.table format.
#'
#' @param dat NO DEFAULT. data.frame. Input sample.
#' @param use.cols NO DEFAULT. Vector of column names to use for clustering.
#' @param xdim DEFAULT = 14. Numeric. Number of first level clusters across the x-axis. xdim x ydim = total number of first level FlowSOM clusters.
#' @param ydim DEFAULT = 14. Numeric. Number of first level clusters across the y-axis. xdim x ydim = total number of first level FlowSOM clusters.
#' @param meta.k DEFAULT = 'auto'. If set to 'auto', then number of metaclusters will be determined automatically. Alternatively, can specify the desired number of metaclusters to create. If set to zero (0), no metaclusters will be created.
#' @param clust.seed DEFAULT = 42 Numeric. Clustering seed for reproducibility.
#' @param meta.seed DEFAULT = 42 Numeric. Metaclustering seed for reproducibility.
#' @param clust.name DEFAULT = "FlowSOM_cluster". Character. Name of the resulting 'cluster' parameter.
#' @param meta.clust.name DEFAULT = "FlowSOM_metacluster". Character. Name of the resulting 'metacluster' parameter.
#' @param mem.ctrl DEFAULT = TRUE. Runs gc() (garbage collection) after a number of steps to free up memory that hasn't been released quickly enough.
#'
#' @usage run.flowsom(dat, use.cols, meta.k, xdim, ydim, clust.seed, meta.seed, clust.name, meta.clust.name)
#'
#' @examples
#' # Run FlowSOM on demonstration dataset
#' res <- Spectre::run.flowsom(Spectre::demo.asih,
#'                             use.cols = names(demo.asinh)[c(11:19)])
#'
#' @author
#' Thomas Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#' Felix Marsh-Wakefield, \email{felix.marsh-wakefield@@sydney.edu.au}
#'
#' @import data.table
#'
#' @export

run.flowsom <- function(dat,
                        use.cols, # names of columns to cluster
                        xdim = 14,
                        ydim = 14,
                        meta.k = 'auto',
                        clust.seed = 42,
                        meta.seed = 42,
                        clust.name = "FlowSOM_cluster",
                        meta.clust.name = "FlowSOM_metacluster",
                        mem.ctrl = TRUE){

  ### Check that necessary packages are installed
      if(!is.element('Spectre', installed.packages()[,1])) stop('Spectre is required but not installed')
      if(!is.element('flowCore', installed.packages()[,1])) stop('flowCore is required but not installed')
      if(!is.element('Biobase', installed.packages()[,1])) stop('Biobase is required but not installed')
      if(!is.element('FlowSOM', installed.packages()[,1])) stop('FlowSOM is required but not installed')
      if(!is.element('data.table', installed.packages()[,1])) stop('data.table is required but not installed')

  ### Require packages
      require(Spectre)
      require(flowCore)
      require(Biobase)
      require(FlowSOM)
      require(data.table)

  ### Test
      # dat <- Spectre::demo.asinh
      # use.cols <- names(demo.asinh)[c(11:19)]
      # xdim = 14
      # ydim = 14
      # meta.k = 'auto'
      # clust.seed = 42
      # meta.seed = 42
      # clust.name = "FlowSOM_cluster"
      # meta.clust.name = "FlowSOM_metacluster"
      # mem.ctrl = TRUE

  ### Prepare starting and using data
      message("Preparing data")

      dat.start <- dat
      clustering.cols <- use.cols

  ### Check selected columns are numeric
      if(any(unlist(lapply(dat[, use.cols, with = FALSE], is.numeric)) == FALSE)) {
        stop('Non-numeric column selected for analysis. Check use.cols parameter.')
      }

      dat <- dat[,use.cols, with = FALSE]

  ### Create flowFrame metadata (column names with descriptions) plus flowFrame
      metadata <- data.frame(name=dimnames(dat)[[2]], desc=paste('column',dimnames(dat)[[2]],'from dataset'))
      dat.ff <- new("flowFrame",
                     exprs=as.matrix(dat), # in order to create a flow frame, data needs to be read as matrix
                     parameters=Biobase::AnnotatedDataFrame(metadata))

      head(flowCore::exprs(dat.ff))
      dat_FlowSOM <- dat.ff

      rm(dat)
      rm(dat.ff)

      if(mem.ctrl == TRUE){
        gc()
      }

  ### Run FlowSOM clustering
      message("Starting FlowSOM")

      FlowSOM_cols <- clustering.cols # clustering cols
      set.seed(clust.seed) # set seed for reproducibility

      ## run FlowSOM (initial steps prior to meta-clustering)
      FlowSOM_out <- FlowSOM::ReadInput(dat_FlowSOM, transform = FALSE, scale = FALSE)

      FlowSOM_out <- FlowSOM::BuildSOM(FlowSOM_out,
                                       colsToUse = FlowSOM_cols,
                                       xdim = xdim,
                                       ydim = ydim)

      FlowSOM_out <- FlowSOM::BuildMST(FlowSOM_out)

      ## extract cluster labels (pre meta-clustering) from output object
      labels_pre <- FlowSOM_out$map$mapping[, 1]

      if(mem.ctrl == TRUE){
        gc()
      }

      flowsom.res.original <- labels_pre

      ## save ORIGINAL cluster labels
      flowsom.res.original <- data.frame("labels_pre" = labels_pre)
      colnames(flowsom.res.original)[grepl('labels_pre',colnames(flowsom.res.original))] <- clust.name

      dat.start <- cbind(dat.start, flowsom.res.original)   # Add results to dat

  ### Metaclustering

      if(meta.k != 0) {

        ## Auto number of MCs
        if(meta.k == 'auto'){
          FlowSOM_out_meta <- FlowSOM::MetaClustering(FlowSOM_out$map$codes,
                                                      method="metaClustering_consensus", 
                                                      seed=meta.seed)
        }

        ## Define number of MCs
        if(meta.k != 'auto'){
          FlowSOM_out_meta <- FlowSOM::metaClustering_consensus(FlowSOM_out$map$codes, k = meta.k, seed = meta.seed)
        }

        labels <- FlowSOM_out_meta[labels_pre]
        flowsom.res.meta <- data.frame("labels" = labels)
        colnames(flowsom.res.meta)[grepl('labels',colnames(flowsom.res.meta))] <- meta.clust.name

        message("Binding metacluster labels to starting dataset")
        dat.start <- cbind(dat.start, flowsom.res.meta)       # Add results to dat
        rm(flowsom.res.meta)
      }

  ### Return
      if(mem.ctrl == TRUE){gc()}
      message("Binding cluster labels to starting dataset")
      dat.start <- data.table::as.data.table(dat.start) # Make dat a data.table for future manipulation

      return(dat.start)
}
