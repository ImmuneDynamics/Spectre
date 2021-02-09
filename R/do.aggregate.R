#' Aggregate data using data.table functions -- summarises the mean, median, or sum expression of each marker on each cluster/population
#'
#' @usage do.aggregate(dat, use.cols, by, func, dt)
#'
#' @param dat NO DEFAULT. A data.table (or data.frame) containing the data to be aggregated.
#' @param use.cols NO DEFAULT. The columns to be used for aggregation.
#' @param by NO DEFAULT. The column that contains the cluster/population identities.
#' @param func DEFAULT = 'median'. Character, method of aggregation. Can be 'median', 'mean', or 'sum'.
#' @param dt DEFAULT = TRUE. If TRUE returns data.table with observations and the use.cols, if false returns data.frame with use.cols, and rownames as the observations
#'
#' @return Returns a new data.table with values of use.cols aggregated per cluster/population.
#'
#' @examples
#' exp <- do.aggregate(dat = Spectre::demo.clustered,
#'                     use.cols = names(demo.clustered)[c(11:19)],
#'                     by = "Population")
#'
#' # Typically followed up by making an expression heatmap
#' Spectre::make.pheatmap(dat = exp,
#'                        file.name = "Pheatmap.png",
#'                        plot.title = "Pheatmap",
#'                        sample.col = "Population",
#'                        plot.cols = names(exp)[c(2:10)])
#'
#' @import data.table
#'
#' @export

do.aggregate <- function(dat,
                             use.cols,
                             by,
                             func = 'median',
                             dt = TRUE
) {

  ### Testing
      # dat <- as.data.table(demo.umap)
      # use.cols <- c("UMAP_42_X", "UMAP_42_Y")
      # by = "FlowSOM_metacluster"
      # func = 'median'
      # dt = FALSE

  ### Setup

      if(func == 'median'){
        res <- dat[, lapply(.SD, 'median', na.rm=TRUE), by=by, .SDcols=use.cols  ]
      }

      if(func == 'mean'){
        res <- dat[, lapply(.SD, 'mean', na.rm=TRUE), by=by, .SDcols=use.cols  ]
      }

      if(func == 'sum'){
        res <- dat[, lapply(.SD, 'sum', na.rm=TRUE), by=by, .SDcols=use.cols  ]
      }

  ### Setup result

      if(dt == TRUE){
        res
      }

      if(dt == FALSE){
        rnms <- res[[by]]
        res <- as.data.frame(res)
        rownames(res) <- rnms
        res[[by]] <- NULL
        res
      }

  ### Return
      return(res)

}
