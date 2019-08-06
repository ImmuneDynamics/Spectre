

# write.files

write.files <- function(x,      # data to save
                        file.name,
                        by.all,

                        by.sample,
                        sample.col,

                        by.group,
                        group.col,

                        write.csv,
                        write.fcs){


  ####### TESTING
    #x <- cell.dat.sub
    #file.name = "TAXXX"
    #by.all = TRUE

    #by.samples = TRUE
    #sample.col = "FileName"

    #by.groups = TRUE
    #group.col = "GroupName"

    #write.csv = TRUE
    #write.fcs = TRUE


  ### All samples ####
    if(by.all == TRUE){

      ## Write CSV (using fwrite)
      if(write.csv == TRUE){fwrite(x = x, file = paste0(file.name, ".csv"))}

      ## Write FCS
      if(write.fcs == TRUE){
        ## Check data and data column names
        head(x)
        dimnames(x)[[2]]

        ## Create FCS file metadata - column names with descriptions
        metadata <- data.frame(name=dimnames(x)[[2]], desc=paste('column',dimnames(x)[[2]],'from dataset'))

        ## Create FCS file metadata - ranges, min, and max settings
        #metadata$range <- apply(apply(x,2,range),2,diff)
        metadata$minRange <- apply(x,2,min)
        metadata$maxRange <- apply(x,2,max)

        ## Create flowframe with tSNE data
        x.ff <- new("flowFrame", exprs=as.matrix(x), parameters=AnnotatedDataFrame(metadata))

        ## Save flowframe as .fcs file -- save data (with new tSNE parameters) as FCS
        new_file_name_fcs <- paste0(file.name, ".fcs")
        write.FCS(x.ff, new_file_name_fcs)
      }
    }


  ### By sample
    if(by.sample == TRUE){
      sample.list <- unique(x[[sample.col]])

      ## Subset data
      for(a in sample.list){
        data_subset <- subset(x, x[[sample.col]] == a)
        dim(data_subset)

          ## Write CSV (using fwrite)
          if(write.csv == TRUE){fwrite(x = data_subset, file = paste0(file.name, "_", "Sample_", a, ".csv"))}

          ## Write FCS
          if(write.fcs == TRUE){

            ## Check data and data column names
            head(data_subset)
            dimnames(data_subset)[[2]]

            ## Create FCS file metadata - column names with descriptions
            metadata <- data.frame(name=dimnames(data_subset)[[2]], desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))

            ## Create FCS file metadata - ranges, min, and max settings
            #metadata$range <- apply(apply(data_subset,2,range),2,diff)
            metadata$minRange <- apply(data_subset,2,min)
            metadata$maxRange <- apply(data_subset,2,max)

            ## Create flowframe with tSNE data
            data_subset.ff <- new("flowFrame", exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata))

            ## Save flowframe as .fcs file -- save data (with new tSNE parameters) as FCS
            write.FCS(data_subset.ff, paste0(file.name, "_", "Sample_", a, ".fcs"))
          }
        }
      }

  ### By group

    if(by.group == TRUE){
      group.list <- unique(x[[group.col]])

      ## Subset data
      for(a in group.list){
        data_subset <- subset(x, x[[group.col]] == a)
        dim(data_subset)

        ## Write CSV (using fwrite)
        if(write.csv == TRUE){fwrite(x = data_subset, file = paste0(file.name, "_", "Group_", a, ".csv"))}

        ## Write FCS
        if(write.fcs == TRUE){

          ## Check data and data column names
          head(data_subset)
          dimnames(data_subset)[[2]]

          ## Create FCS file metadata - column names with descriptions
          metadata <- data.frame(name=dimnames(data_subset)[[2]], desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))

          ## Create FCS file metadata - ranges, min, and max settings
          #metadata$range <- apply(apply(data_subset,2,range),2,diff)
          metadata$minRange <- apply(data_subset,2,min)
          metadata$maxRange <- apply(data_subset,2,max)

          ## Create flowframe with tSNE data
          data_subset.ff <- new("flowFrame", exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata))

          ## Save flowframe as .fcs file -- save data (with new tSNE parameters) as FCS
          write.FCS(data_subset.ff, paste0(file.name, "_", "Group_", a, ".fcs"))
        }
      }
    }

}





