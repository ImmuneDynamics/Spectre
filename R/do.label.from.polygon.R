#' do.label.from.polygon
#' 
#' @import data.table
#' 
#' @export

do.label.from.polygon <- function(spatial.dat,
                                  cell.dat,
                                  mask,
                                  labels,
                                  
                                  name = "Label",
                                  id.col = 'ID',
                                  roi.col = 'ROI'#,
                                  #x.col = 'x',
                                  #y.col = 'y'
){
  
  ###
  
  # spatial.dat <- spatial.dat
  # cell.dat <- cell.dat
  # mask <- 'obj_mask'
  #
  # name = "Label"
  # labels <- c("HTC", "CTL", "Background", "Hepatocyte", "cDC", "CD11bpos", "Other immune", "Macrophages")
  #
  # id.col = 'ID'
  # roi.col = 'ROI'
  # x.col = 'x'
  # y.col = 'y'
  
  TEMP_LABEL_PLACEHOLDER1 <- labels
  TEMP_LABEL_PLACEHOLDER1 <- as.data.table(TEMP_LABEL_PLACEHOLDER1)
  
  ###
  rois <- names(spatial.dat)
  
  start.dat.list <- list()
  dat.list <- list()
  mask.list <- list()
  
  unq.labels <- list()
  res.list <- list()
  
  ###
  
  for(i in rois){
    unq.labels[[i]] <- unique(as.data.frame(spatial.dat[[i]]$MASKS[[mask]]$polygons))
    rm(i)
  }
  
  all.labs <- unique(rbindlist(unq.labels, fill = TRUE))
  all.labs <- all.labs[order(all.labs)]
  
  all.labs <- cbind(all.labs, TEMP_LABEL_PLACEHOLDER1)
  all.labs
  
  for(i in rois){
    # i <- rois[[8]]
    dat.list[[i]] <- cell.dat[cell.dat[[roi.col]] == i,]
    coordinates(dat.list[[i]]) <- ~ x + y
    
    mask.list[[i]] <- spatial.dat[[i]]$MASKS[[mask]]$polygons
    proj4string(dat.list[[i]]) <- proj4string(mask.list[[i]])
    
    res <- over(dat.list[[i]], mask.list[[i]])
    res.list[[i]] <- res[[1]]
  }
  
  ###
  
  res.dt <- unlist(res.list)
  res.dt <- as.data.table(res.dt)
  res.dt
  
  if(nrow(cell.dat) != nrow(res.dt)){
    stop("Result dt rows are inconsistent with starting dt rows")
  }
  
  res.dt <- do.embed.columns(res.dt, "res.dt", all.labs, 'obj_mask')
  res.dt
  
  res.dt$res.dt <- NULL
  names(res.dt)[length(names(res.dt))] <- name
  res.dt
  
  return.dat <- cbind(cell.dat, res.dt)
  
  ###
  return(return.dat)
}


