# make.sample.table


make.sample.table <- function(x, # takes a single dataframe or datatable
                               sample.col.name, # Character, name of the column that specifies sample names # REQUIRED

                               include.groups = FALSE, # Logical, name of the column that specifies sample names
                               group.col.name = NULL, # Character, name of the column that specifies group names

                               include.batch = FALSE,
                               batch.col.name = NULL
                               ){

    ### For testing

        #x = cell.dat
        #sample.col.name = "Sample"
        #include.groups = TRUE
        #group.col.name = "Group"
        #include.batch = TRUE
        #batch.col.name = "Batch"

    ### Create list of sample names
        all.sample.names <- as.matrix(unique(x[sample.col.name])[,1])

        if(include.groups == TRUE){
          all.group.names <- as.matrix(unique(x[group.col.name])[,1])
        }

        if(include.batch == TRUE){
          all.batch.names <- as.matrix(unique(x[batch.col.name])[,1])
        }

    ### Cells per sample
        cells.per.sample = list()

        for(i in c(1:length(all.sample.names))){
          a <- unique(all.sample.names)[i]
          cells.per.sample[[i]] <- nrow(subset(x, x[, samp.col] == a))
        }

        cells.per.sample <- as.matrix(cells.per.sample)


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
          group.mismatch.res <- as.matrix(group.sample.match)
          group.mismatch.res <- (group.mismatch.res > 1) # Returns TRUE for each sample that has more than one group assignment
        }


    ### IF USING BATCHES, create list of group allocations and add them to sample.table
        if(include.batch == TRUE){

          ## Create empty lists
          batch.allocations = list()
          batch.sample.match = list()

          ## Loop
          for(i in c(1:length(all.sample.names))){
            a <- unique(all.sample.names)[i]
            temp <- subset(x, x[, samp.col] == a)

            ## Add a check to see if all are unique
            batch.sample.match[[i]] <- length(unique(temp[batch.col]))

            # uses the batch from the first row
            batch.allocations[[i]] <- temp[[batch.col]][1]
          }
          batch.allocations <- as.matrix(batch.allocations)

          ## Create group mismatch results
          batch.mismatch.res <- as.matrix(batch.sample.match)
          batch.mismatch.res <- (batch.mismatch.res > 1) # Returns TRUE for each sample that has more than one group assignment
        }


    ### Create table

        ## Create sample.table with cell counts/sample
          sample.table <- data.frame(all.sample.names, cells.per.sample)
          colnames(sample.table) <- c(sample.col.name, "Rows per sample")

        ## Add group data, if used

          ## Group mismatch check
          if(all(group.mismatch.res) == TRUE){
            print("Some of your samples have more than one group assigned to the cells from a single sample. Omitting group allocation data from the table.")
          }

          if(all(group.mismatch.res) == FALSE){
            sample.table[group.col.name] <- as.data.frame(group.allocations)
          }

        ## Add batch data, if used

          if(all(batch.mismatch.res) == TRUE){
            print("Some of your samples have more than one batch assigned to the cells from a single sample. Omitting batch allocation data from the table.")
          }

          if(all(batch.mismatch.res) == FALSE){
            sample.table[batch.col.name] <- as.data.frame(batch.allocations)
          }

    ### Assign
        assign("sample.table", sample.table, envir = globalenv())
        print("'sample.table' has been created")

        assign("all.sample.names", all.sample.names, envir = globalenv())
        print("'all.sample.names' has been created")

        if(include.groups == TRUE){
          if(all(group.mismatch.res) == FALSE){
            assign("all.group.names", all.group.names, envir = globalenv())
            print("'all.group.names' has been created")
          }
        }

        if(include.batch == TRUE){
          if(all(batch.mismatch.res) == FALSE){
            assign("all.batch.names", all.batch.names, envir = globalenv())
            print("'all.batch.names' has been created")
          }
        }
}

