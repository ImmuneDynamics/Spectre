#' Run FlowSOM
#'
#' Method to run the FlowSOM clustering algorithm.
#' This function runs FlowSOM on a data.table with cells (rows) vs markers
#' (columns) with new columns for FlowSOM clusters and metaclusters.
#' Output data will be "flowsom.res.original" (for clusters) and 
#' "flowsom.res.meta" (for metaclusters).
#' Uses the R packages "FlowSOM" for clustering, "flowCore" for handling
#' .fcs files, "Biobase" for creating a flow frame, "data.table" for
#' handling data.table format.
#'
#' @param dat NO DEFAULT. data.frame. Input sample.
#' @param use.cols NO DEFAULT. Vector of column names to use for clustering.
#' @param xdim DEFAULT = 14. Numeric. Number of first level clusters across 
#' the x-axis. xdim x ydim = total number of first level FlowSOM clusters.
#' @param ydim DEFAULT = 14. Numeric. Number of first level clusters across 
#' the y-axis. xdim x ydim = total number of first level FlowSOM clusters.
#' @param meta.k DEFAULT = 'auto'. If set to 'auto', then number of 
#' metaclusters will be determined automatically. Alternatively, can specify 
#' the desired number of metaclusters to create. If set to zero (0), 
#' no metaclusters will be created.
#' @param max.meta DEFAULT = 20. Only used if meta.k is set to 'auto'. 
#' This parameter indicates the maximum number of metaclusters FlowSOM will 
#' try out when determining the optimal number of metaclusters for the dataset.
#' @param clust.seed DEFAULT = 42 Numeric. Clustering seed for reproducibility.
#' @param meta.seed DEFAULT = 42 Numeric. Metaclustering seed for 
#' reproducibility.
#' @param clust.name DEFAULT = "FlowSOM_cluster". Character. Name of the 
#' resulting 'cluster' parameter.
#' @param meta.clust.name DEFAULT = "FlowSOM_metacluster". Character. 
#' Name of the resulting 'metacluster' parameter.
#' @param mem.ctrl DEFAULT = TRUE. Runs gc() (garbage collection) after a 
#' number of steps to free up memory that hasn't been released quickly enough.
#' @param verbose DEFAULT = TRUE. Logical. Whether to print progress messages.
#' 
#' @usage run.flowsom(dat, use.cols, xdim=14, ydim=14, meta.k='auto', 
#' max.meta=20, clust.seed=42, meta.seed=42, clust.name="FlowSOM_cluster", 
#' meta.clust.name="FlowSOM_metacluster", mem.ctrl=TRUE, verbose=TRUE)
#'
#' @examples
#' # Run FlowSOM on demonstration dataset
#' res <- Spectre::run.flowsom(Spectre::demo.clustered,
#' use.cols = c("NK11_asinh", "CD3_asinh", 
#' "CD45_asinh", "Ly6G_asinh", "CD11b_asinh", 
#' "B220_asinh", "CD8a_asinh", "Ly6C_asinh", 
#' "CD4_asinh"))
#'
#' @author
#' Thomas Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#' Felix Marsh-Wakefield, \email{felix.marsh-wakefield@@sydney.edu.au}
#'
#' @import data.table
#' @import FlowSOM
#' @import flowCore
#'
#' @export
#' 
run.flowsom <- function(
    dat, use.cols,
    xdim = 14,
    ydim = 14,
    meta.k = "auto",
    max.meta = 20,
    clust.seed = 42,
    meta.seed = 42,
    clust.name = "FlowSOM_cluster",
    meta.clust.name = "FlowSOM_metacluster",
    mem.ctrl = TRUE,
    verbose = TRUE
) {
    ### Test
    # dat <- Spectre::demo.clustered
    # use.cols <- names(Spectre::demo.clustered)[c(11:19)]
    # xdim = 14
    # ydim = 14
    # meta.k = 'auto'
    # clust.seed = 42
    # meta.seed = 42
    # clust.name = "FlowSOM_cluster"
    # meta.clust.name = "FlowSOM_metacluster"
    # mem.ctrl = TRUE

    ### Prepare starting and using data
    if (verbose) {
        message("Preparing data")
    }

    ### Check selected columns are numeric
    .check_numeric_columns(dat, use.cols)

    # Create flowFrame metadata (column names with descriptions) 
    # plus flowFrame object

    metadata <- data.frame(
        name = use.cols, 
        desc = paste("column", use.cols, "from dataset")
    )
    # don't use flowcore::flowFrame() as it set min and max and range
    # values which screws up clustering.
    dat.ff <- new(
        "flowFrame",
        exprs = as.matrix(dat[, use.cols, with = FALSE]),
        parameters = Biobase::AnnotatedDataFrame(metadata)
    )

    ### Run FlowSOM clustering
    
    if (verbose) {
        message("Starting FlowSOM")
    }

    set.seed(clust.seed) # set seed for reproducibility

    ## run FlowSOM (initial steps prior to meta-clustering)
    FlowSOM_out <- FlowSOM::ReadInput(
        dat.ff, 
        transform = FALSE, 
        scale = FALSE,
        silent = !verbose
    )

    if (mem.ctrl) {
        if (verbose) {
            message("Freeing up memory")
        }
        rm(dat.ff)
        rm(metadata)
        gc()
    }

    FlowSOM_out <- FlowSOM::BuildSOM(
        FlowSOM_out,
        colsToUse = use.cols,
        xdim = xdim,
        ydim = ydim,
        silent = !verbose
    )

    FlowSOM_out <- FlowSOM::BuildMST(
        FlowSOM_out,
        silent = !verbose
    )

    ## extract cluster labels (pre meta-clustering) from output object
    if (verbose) {
        message("Adding cluster labels to input dataset")
    }

    cluster_labels <- FlowSOM_out$map$mapping[, 1]
    dat_out <- .add_or_replace_column(
        dt = data.table::copy(dat),
        col_name = clust.name,
        values = as.factor(cluster_labels)
    )

    ### Metaclustering

    if (meta.k != 0) {
        if (verbose) {
            message("Starting metaclustering")
        }

        if (meta.k == "auto") {
            FlowSOM_out_meta <- FlowSOM::MetaClustering(
                FlowSOM_out$map$codes,
                method = "metaClustering_consensus",
                max = max.meta,
                seed = meta.seed
            )
        } else {
            FlowSOM_out_meta <- FlowSOM::metaClustering_consensus(
                FlowSOM_out$map$codes, 
                k = meta.k, 
                seed = meta.seed
            )
        }
        
        # FlowSOM_out_meta will label each "cluster" with a metacluster id
        # not label each cell with a metacluster id
        # Hence, need to map the metacluster id to the cells via the 
        # "cluster" labels.
        meta_labels <- FlowSOM_out_meta[cluster_labels]

        if (verbose) {
            message("Adding metacluster labels to input dataset")
        }

        dat_out <- .add_or_replace_column(
            dt = dat_out,
            col_name = meta.clust.name,
            values = as.factor(meta_labels)
        )
    }

    return(dat_out)
}
