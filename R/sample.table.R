# make.sample.table


make.sample.table <- function(x, # takes a single dataframe or datatable
                             sample.col.name, # Character, name of the column that specifies sample names # REQUIRED

                             include.groups, # Logical, name of the column that specifies sample names
                             group.col.name # Character, name of the column that specifies group names
                             ){


    ### Create list of sample names
        all.sample.names <- as.matrix(unique(x[sample.col.name])[,1])
        all.group.names <- as.matrix(unique(x[group.col.name])[,1])

    ### Cells per sample
        cells.per.sample = list()

        for(i in c(1:length(all.sample.names))){
          a <- unique(all.sample.names)[i]
          cells.per.sample[[i]] <- nrow(subset(x, x[, samp.col] == a))
        }

        cells.per.sample <- as.matrix(cells.per.sample)

    ### IF NO GROUPS, create sample.table
        if(include.groups == FALSE){
          sample.table <- data.frame(all.sample.names, cells.per.sample)
        }

    ### IF USING GROUPS, create list of group allocations and add them to sample.table
        if(include.groups == TRUE){

          ## Create empty lists
          group.allocations = list()
          group.sample.match = list()

          ## Loop
          for(i in c(1:length(all.sample.names))){
            a <- unique(all.sample.names)[i]
            temp <- subset(x, x[, samp.col] == a)

            ## Add a check to see if all are unique
            group.sample.match[[i]] <- length(unique(temp[grp.col]))

            # uses the group from the first row
            group.allocations[[i]] <- temp[[grp.col]][1]
          }

          group.allocations <- as.matrix(group.allocations)

          ## Create group mismatch results
          mismatch.res <- as.matrix(group.sample.match)
          mismatch.res <- (mismatch.res > 1) # Returns TRUE for each sample that has more than one group assignment

          ## Group mismatch check
          if(all(mismatch.res) == TRUE){
            print("Some of your samples have more than one group assigned to the cells from a single sample. Omitting group allocation data from the table.")
            sample.table <- data.frame(all.sample.names, cells.per.sample)
            colnames(sample.table) <- c("Sample Name", "Cells per sample")
            }

          if(all(mismatch.res) == FALSE){
            sample.table <- data.frame(all.sample.names, group.allocations, cells.per.sample)
            colnames(sample.table) <- c("Sample Name", "Group", "Cells per sample")

            }
          }


    ### Assign
        assign("sample.table", sample.table, envir = globalenv())
        print("'sample.table' has been created")

        assign("all.sample.names", all.sample.names, envir = globalenv())
        print("'all.sample.names' has been created")

        if(include.groups == TRUE){
          if(all(mismatch.res) == FALSE){
            assign("all.group.names", all.group.names, envir = globalenv())
            print("'all.group.names' has been created")
          }
        }
}

