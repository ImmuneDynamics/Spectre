do.logicle <- function(dat, 
                       use.cols,
                       w = 0.5,
                       t = 262144,
                       m = 4.5,
                       a = 0) {
  
  require(flowCore)
  require(Biobase)
  require(data.table)
  
  values <- dat[,use.cols, with=FALSE]
  
  val.mat <- as.matrix(values)
  
  ff <- new("flowFrame", exprs = val.mat)
  
  if (is.null(w) || is.null(t) || is.null(m) || is.null(a)) {
    # Automatic estimated logicle function
    lgcl <- estimateLogicle(ff, channels = use.cols) 
  } else {
    # User defined logicle function
    trans_lgcl <- logicleTransform(w = w, t = t, m = m, a = a)
    lgcl <- transformList(use.cols, trans_lgcl)
  }
  
  # Do transformation
  trans.val <- transform(ff, lgcl)
  
  trans.val.dt <- data.table(exprs(trans.val))
  names(trans.val.dt) <- paste(names(trans.val.dt), "logicle", sep = "_")
  
  dat_copy <- cbind(dat, trans.val.dt)
  return(dat_copy)
}