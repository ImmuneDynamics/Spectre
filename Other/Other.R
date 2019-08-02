
### CytoTools::write.files

write.files <- function(x,
                        annot,
                        prim.dir,
                        write.csv,
                        write.fcs,

                        write.all,
                        write.samples,
                        write.groups,

                        samp.col,
                        grp.col)



## Create new output directory
tim <- Sys.time()
tim <- gsub(":", "-", tim)
tim <- gsub(" ", "_", tim)

newdir <- paste0(tim, "_", annot)

setwd(prim.dir)
dir.create(paste0(newdir), showWarnings = FALSE)
setwd(newdir)

## Write 'all' data
if(write.tSNE.merged == 1){

  ## Save data (with new tSNE parameters) as CSV
  write.csv(x = data_subsampled, file = paste0(data.name, "_with_tSNE", ".csv"), row.names=FALSE)

  ## Check data and data column names
  head(data_subsampled)
  dimnames(data_subsampled)[[2]]

  ## Create FCS file metadata - column names with descriptions
  metadata <- data.frame(name=dimnames(data_subsampled)[[2]],
                         desc=paste('column',dimnames(data_subsampled)[[2]],'from dataset')
  )

  ## Create FCS file metadata - ranges, min, and max settings
  #metadata$range <- apply(apply(data_subsampled,2,range),2,diff)
  metadata$minRange <- apply(data_subsampled,2,min)
  metadata$maxRange <- apply(data_subsampled,2,max)

  ## Create flowframe with tSNE data
  data_subsampled.ff <- new("flowFrame",
                            exprs=as.matrix(data_subsampled), # in order to create a flow frame, data needs to be read as matrix
                            parameters=AnnotatedDataFrame(metadata)
  )

  ## Check flow frame
  data_subsampled.ff

  ## Save flowframe as .fcs file -- save data (with new tSNE parameters) as FCS
  new_file_name_fcs <- paste0(data.name, "_with_tSNE", ".fcs")
  write.FCS(data_subsampled.ff, new_file_name_fcs)

  ### There is a delay here -- fcs file ends up in primary directory
}

### 3.7 - Write GROUPED data (with tSNE and FlowSOM parameters) to .csv and .fcs
if(write.tSNE.group == 1){
  for(a in AllGroupNames){
    data_subset <- subset(data_subsampled, data_subsampled[[grp.col]] == a)
    dim(data_subsampled)

    ## write .csv
    write.csv(data_subset, file = paste0(data.name, "_", a, "_with_tSNE.csv", sep = ""), row.names=FALSE)

    ## write .fcs
    metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))
    metadata

    ## Create FCS file metadata - ranges, min, and max settings
    #metadata$range <- apply(apply(data_subset,2,range),2,diff)
    metadata$minRange <- apply(data_subset,2,min)
    metadata$maxRange <- apply(data_subset,2,max)

    data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
    head(data_subset.ff)
    write.FCS(data_subset.ff, paste0(data.name, "_", a, "_with_tSNE", ".fcs"))
  }
}

### 3.8 - Write individual files (with tSNE and FlowSOM parameters) to .csv and .fcs

if(write.tSNE.sep == 1){

  for (a in AllSampleNames) {
    data_subset <- subset(data_subsampled, data_subsampled[[samp.col]] == a)

    ## write .csv
    write.csv(data_subset, file = paste0(data.name, "_", a, "_with_tSNE.csv", sep = ""), row.names=FALSE)

    ## write .fcs
    metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))

    ## Create FCS file metadata - ranges, min, and max settings
    #metadata$range <- apply(apply(data_subset,2,range),2,diff)
    metadata$minRange <- apply(data_subset,2,min)
    metadata$maxRange <- apply(data_subset,2,max)

    data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
    head(data_subset.ff)
    write.FCS(data_subset.ff, paste0(data.name, "_", a, "_with_tSNE", ".fcs"))
  }

  ## Move back to PrimaryDirectory
  setwd(prim.dir)
  getwd()
}

# note -- in the final data output, all parameters are included, but only the subsampled and/or transformed cells

}

