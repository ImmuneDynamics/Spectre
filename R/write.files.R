# write.files

write.files <- function(x,      # data to save
                        file.prefix,

                        divide.by = NULL,

                        write.csv = TRUE,
                        write.fcs = FALSE){

  ####### TESTING

    # x <- cell.dat
    # file.prefix = "Demo"
    #
    # divide.by = "Sample"
    #
    # write.csv = TRUE
    # write.fcs = TRUE


  ### If writing FCS files, setup numeric only data frames

      if(write.fcs == TRUE){
        ## Convert character entries into numeric -- so that it can be written as an FCS file
            #nums <- unlist(lapply(x, is.numeric))
            #char <- unlist(lapply(x, is.character))

        w <- x[sapply(x, function(x) !is.numeric(x))] # select all character and factor columns
        head(w)

        for(a in colnames(w)){w[[a]] <- as.numeric(factor(w[[a]]))}
        head(w)

        x.for.fcs <- x
        x.for.fcs[sapply(x.for.fcs, function(x.for.fcs) !is.numeric(x.for.fcs))] <- w

        head(x.for.fcs)
        head(x)
      }


  ### If divide.by = NULL, then write all data

    if(is.null(divide.by) == TRUE){
      ## Write CSV (using fwrite)
      if(write.csv == TRUE){fwrite(x = x, file = paste0(file.prefix, "_all_data.csv"))}

      ## Write FCS
      if(write.fcs == TRUE){

        ## Check data and data column names
        head(x.for.fcs)
        dimnames(x.for.fcs)[[2]]

        ## Create FCS file metadata - column names with descriptions
        metadata <- data.frame(name=dimnames(x.for.fcs)[[2]], desc=paste('column',dimnames(x.for.fcs)[[2]],'from dataset'))

        ## Create FCS file metadata - ranges, min, and max settings
        metadata$range <- apply(apply(x.for.fcs,2,range),2,diff)
        metadata$minRange <- apply(x.for.fcs,2,min)
        metadata$maxRange <- apply(x.for.fcs,2,max)

        ## Create flowframe with tSNE data
        x.ff <- new("flowFrame", exprs=as.matrix(x.for.fcs), parameters=AnnotatedDataFrame(metadata))

        ## Save flowframe as .fcs file -- save data (with new tSNE parameters) as FCS
        new_file_name_fcs <- paste0(file.prefix, "_all_data.fcs")
        write.FCS(x.ff, new_file_name_fcs)
      }
    }



  ### If dividing by a column

    if(is.null(divide.by) == FALSE){

      ## Prep keyword data to divide by
      divide.list <- unique(x[[divide.by]])

        ## Write CSV (using fwrite)
        if(write.csv == TRUE){
          for(a in divide.list){
            data_subset <- subset(x, x[[divide.by]] == a)
            dim(data_subset)
            fwrite(x = data_subset, file = paste0(file.prefix, "_", divide.by, "_", a, ".csv"))
          }
        }

        ## Write FCS
        if(write.fcs == TRUE){

          divide.list.fcs <- unique(x.for.fcs[[divide.by]])

          for(a in divide.list.fcs){
            data_subset <- subset(x.for.fcs, x.for.fcs[[divide.by]] == a)
            dim(data_subset)

            ## Check data and data column names
            head(data_subset)
            dimnames(data_subset)[[2]]

            ## Create FCS file metadata - column names with descriptions
            metadata <- data.frame(name=dimnames(data_subset)[[2]], desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))

            ## Create FCS file metadata - ranges, min, and max settings
            metadata$range <- apply(apply(data_subset,2,range),2,diff)
            metadata$minRange <- apply(data_subset,2,min)
            metadata$maxRange <- apply(data_subset,2,max)

            ## Create flowframe with tSNE data
            data_subset.ff <- new("flowFrame", exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata))

            ## Create flowframe with tSNE data
            data_subset.ff <- new("flowFrame", exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata))

            ## Save flowframe as .fcs file -- save data (with new tSNE parameters) as FCS
            write.FCS(data_subset.ff,  paste0(file.prefix, "_", divide.by, "_", divide.list[[a]], ".fcs"))
            }
        }
    }
}





