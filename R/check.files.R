#' check.files
#' 
#' @param dat List of data.tables, typically one per 'sample'
#' 
#' @export

check.files <- function(dat){
  
  ###
  
  files <- names(dat)
  
  file.table <- data.table('FileName' = files)
  setorderv(file.table, "FileName")
  file.table
  
  colcheck.list <- list()
  colname.list <- list()
  
  ###
  
  for(i in files){
    colcheck.list[[i]] <- dat[[i]][1,]
    colname.list[[i]] <- names(dat[[i]])
  }
  
  colcheck.list
  colname.list
  
  ### Condensing column check data
  
  colcheck.dat <- rbindlist(colcheck.list, fill = TRUE)
  colcheck.dat <- as.data.table(!is.na(colcheck.dat))
  colcheck.dat
  
  ### Condensing column names summary
  
  all.colnames <- unique(unlist(colname.list))
  all.colnames
  
  colnames.summary <- list()
  
  for(a in all.colnames){
    # a <- all.colnames[1]
    colnames.summary[[a]] <- length(which(colcheck.dat[[a]]))
  }
  
  colnames.summary <- rbindlist(list(colnames.summary))
  colnames.summary
  
  ###
  
  colcheck.dat$FileName <- files
  colcheck.dat
  
  ### Making colname table
  
  colname.table <- data.table('Column name' = all.colnames, 'New name' = all.colnames)
  colname.table
  
  ### Condensing test data
  
  tempdat <- rbindlist(dat, fill = TRUE)
  tempdat
  
  ### Output
  
  res <- list('Column check' = colcheck.dat,
              'Column summary' = colnames.summary,
              'All column names' = all.colnames,
              'Column name table' = colname.table,
              'Test data' = tempdat)
  
  ### Return
  
  return(res)
}