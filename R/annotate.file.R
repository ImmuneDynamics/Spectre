### annotate.file
annotate.files <- function(x, # list of files
                          col.name, # what do you want to name the column
                          d, # annotations to add
                          assignments, # specifies where the annotations are added
                          add.num){ # add a 'number' for each annotation?

  ## For testing
      #x = data.list
      #col.name = "Sample"
      #d = sample.names
      #assignments = sample.assign
      #add.num = TRUE

  ## Add keywords to dataset
  num.of.groups <- c(1:length(d))
  for(a in num.of.groups){
    for(i in c(assignments[[a]])){
      x[[i]][[col.name]] <- d[[a]]

      if(add.num == TRUE){
        n <- paste0(col.name, "_Num")
        x[[i]][[n]] <- a
      }
    }
  }

  x[[1]]
  x[[12]]

  assign("data.list", x, envir = globalenv())
}

