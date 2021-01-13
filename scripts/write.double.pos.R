
write.double.pos <- function(dat,
                             sample.col,
                             pop.col,
                             markers,
                             cutoffs,
                             and.or = 'and'
                             ){
  ### Test data
      # 
      # dat <- Spectre::demo.clustered
      # sample.col <- "Sample"
      # pop.col <- "Population"
      # markers <- c("CD45_asinh", "CD11b_asinh")
      # cutoffs <- c(3, 3)
      # 
      # and.or <- 'and'
  
  ### Loop per sample
  
      res.list <- list()
  
      samps <- unique(dat[[sample.col]])
      
      for(i in samps){
        # i <- samps[[1]]
        temp <- dat[dat[[sample.col]] == i,]
        
        ## Or (either A+, B+, or A+B+)
            if(and.or == 'or'){
              res <- temp[temp[[markers[1]]] > cutoffs[1] |
                          temp[[markers[2]]] > cutoffs[2],]
            }
        
        ## And (i.e. double positive only)
            if(and.or == 'and'){
              x <- temp[temp[[markers[1]]] > cutoffs[1],]
              res <- x[x[[markers[2]]] > cutoffs[2],]
            }
        
        ## Results
            perc.pos <- nrow(res) / nrow(temp) * 100
            res.list[[i]] <- perc.pos
      }
      
      t <- unlist(res.list)
      
      rtn <- cbind(as.data.table(names(t)), as.data.table(t))
      names(rtn) <- c(sample.col, paste0("Percent ", markers[1], " ", and.or ," ", markers[2], " positive"))    
      rtn
      
  ### Return
      return(rtn)
}

