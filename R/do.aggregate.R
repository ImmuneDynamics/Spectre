#' Aggregate data using data.table functions
#'
#' @export

do.aggregate <- function(dat,
                         use.cols,
                         by,
                         func = 'median', # 'median', 'mean', 'sum'
                         dt = TRUE # if TRUE returns data.table with observations and the use.cols, if false returns data.frame with use.cols, and rownames as the observations
) {
  
  ### Demo
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