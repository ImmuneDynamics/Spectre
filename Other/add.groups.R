### add.groups

add.groups <- function(x, # list of files
                       grp.names,
                       grp.nums){

  ## Add 'GroupName' and 'GroupNo' keywords to dataset
  num.of.groups <- c(1:length(grp.names))
  for(a in num.of.groups){
    for(i in c(grp.nums[[a]])){
      x[[i]]$GroupName <- grp.names[[a]]
      x[[i]]$GroupNo <- a
    }
  }
  assign(x = "data.list", value = x, envir = globalenv())
}
