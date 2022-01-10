#' write.files - Write .csv or .fcs files
#'
#' @param dat Dataframe. No default.
#' @param file.prefix Character prefix to add before filename.
#' @param divide.by Character. Do you want to write the whole dataset (use NULL) or split it by sample/group/cluster etc (etner the name of the column to divide by, e.g. "SampleName"). Default NULL. .
#' @param write.csv Defaults to TRUE. Logical. TRUE to write CSV files, FALSE to not.
#' @param write.fcs Defaults to FALSE. Logical. TRUE to write FCS files, FALSE to not.
#'
#' This function writes CSV files from your data using fwrite. Can also write FCS files.
#'
#' @usage write.files(dat, file.prefix, divide.by, write.csv, write.fcs)
#'
#' @export

write.files <- function(dat,      # data to save
                        file.prefix,

                        divide.by = NULL,
                        write.csv = TRUE,
                        write.fcs = FALSE){

  ####### TESTING

      # setwd("/Users/thomasa/Desktop/")
      # dir.create("FCS files")
      # setwd("FCS files")
      #
      # dat <- Spectre::demo.clustered
      # file.prefix = "Demo"
      #
      # divide.by = "Sample"
      #
      # write.csv = TRUE
      # write.fcs = TRUE

  ### Setup

      dat <- as.data.table(dat)

  ### Prep data for FCS files (if required)

      if(write.fcs == TRUE){
        ## Convert character entries into numeric -- so that it can be written as an FCS file
        non.nums <- sapply(dat, function(x) !is.numeric(x))
        nums <- sapply(dat, function(x) is.numeric(x))

        non.num.dat <- dat[,non.nums, with = FALSE] # select all character and factor columns
        for(a in colnames(non.num.dat)){non.num.dat[[a]] <- as.numeric(factor(non.num.dat[[a]]))}
        non.num.dat

        rep.cols <- names(non.num.dat)
        dat.for.fcs <- dat

        for(i in rep.cols){
          dat.for.fcs[[i]] <- non.num.dat[[i]]
        }
      }

  ### Write all data

      if(is.null(divide.by) == TRUE){

        if(write.csv == TRUE){
          fwrite(x = dat, file = paste0(file.prefix, ".csv"))
        }

        if(write.fcs == TRUE){
          head(dat.for.fcs)
          dimnames(dat.for.fcs)[[2]]

          ## Create FCS file metadata - column names with descriptions
          metadata <- data.frame(name=dimnames(dat.for.fcs)[[2]], desc=paste('column',dimnames(dat.for.fcs)[[2]],'from dataset'))

          ## Create FCS file metadata - ranges, min, and max settings
          metadata$range <- apply(apply(dat.for.fcs,2,range),2,diff)
          metadata$minRange <- apply(dat.for.fcs,2,min)
          metadata$maxRange <- apply(dat.for.fcs,2,max)

      ## Create flowframe with tSNE data
      dat.ff <- new("flowFrame", exprs=as.matrix(dat.for.fcs), parameters=Biobase::AnnotatedDataFrame(metadata))

          ## Save flowframe as .fcs file -- save data (with new tSNE parameters) as FCS
          write.FCS(dat.ff, paste0(file.prefix, ".fcs"))
        }
      }

  ### Write divided data

      if(is.null(divide.by) == FALSE){

        divide.list <- unique(dat[[divide.by]])

        if(write.csv == TRUE){
          for(a in divide.list){

            data_subset <- dat[dat[[divide.by]] == a,]
            if(isTRUE(exists("file.prefix"))){fwrite(x = data_subset, file = paste0(file.prefix, "_", divide.by, "_", a, ".csv"))}
            if(isFALSE(exists("file.prefix"))){fwrite(x = data_subset, file = paste0(divide.by, "_", a, ".csv"))}
          }
        }


        if(write.fcs == TRUE){
            if(!is.element('flowCore', installed.packages()[, 1]))
                stop('flowCore is required but not installed')
            
            require(flowCore)
            
            divide.list.fcs <- unique(dat.for.fcs[[divide.by]])
            
            for (a in c(1:length(divide.list.fcs))) {
                nme <- divide.list[[a]]
                num <- divide.list.fcs[[a]]
                
                data_subset_fcs <-
                    dat.for.fcs[dat.for.fcs[[divide.by]] == num, ]
                
                dim(data_subset_fcs)
                
                ## Check data and data column names
                head(data_subset_fcs)
                dimnames(data_subset_fcs)[[2]]
                
                ## Create FCS file metadata - column names with descriptions
                metadata <-
                    data.frame(
                        name = dimnames(data_subset_fcs)[[2]],
                        desc = paste('column', dimnames(data_subset_fcs)[[2]], 'from dataset')
                    )
                
                ## Create FCS file metadata - ranges, min, and max settings
                metadata$range <-
                    apply(apply(data_subset_fcs, 2, range), 2, diff)
                metadata$minRange <- apply(data_subset_fcs, 2, min)
                metadata$maxRange <- apply(data_subset_fcs, 2, max)
                
                data_subset.ff <-
                    new(
                        "flowFrame",
                        exprs = as.matrix(data_subset_fcs),
                        parameters = Biobase::AnnotatedDataFrame(metadata)
                    )
                
                ## Save flowframe as .fcs file -- save data (with new tSNE parameters) as FCS
                if (isTRUE(exists("file.prefix"))) {
                    write.FCS(data_subset.ff,
                              paste0(file.prefix, "_", divide.by, "_", nme, ".fcs"))
                }
                if (isFALSE(exists("file.prefix"))) {
                    write.FCS(data_subset.ff,  paste0(divide.by, "_", nme, ".fcs"))
                }
            }
        }
      }
}
