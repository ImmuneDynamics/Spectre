### merge.files

merge.files <- function(x,
                        remove.duplicates = TRUE
){

  cell.dat <- plyr::rbind.fill(x)

  if(remove.duplicates == TRUE){
    cell.dat <- cell.dat[!duplicated(cell.dat), ]  # remove rows containing duplicate values within rounding
  }

  assign("cell.dat", cell.dat, envir = globalenv())
  ifelse(exists("cell.dat"), "All files merged into 'cell.dat'", "ERROR in merging files into 'cell.dat'")
}
