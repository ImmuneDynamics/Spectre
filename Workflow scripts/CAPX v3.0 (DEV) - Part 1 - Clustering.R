##########################################################################################################
#### DRAFT Cytometry Analysis Pipeline for large and compleX data (CAPX) v3.0 - using the Spectre R package
#### Part 1 - Clustering, dimensionality reduction, save files
##########################################################################################################

    # Thomas Myles Ashhurst, Felix Marsh-Wakefield
    # 2019-08-02
    # Workflow: https://sydneycytometry.org.au/capx
    # Spectre R package: https://sydneycytometry.org.au/spectre

##########################################################################################################
#### 1. Install packages, load packages, and set working directory
##########################################################################################################

    ### 1.1. Load 'Spectre' package (using devtools)
        if(!require('devtools')) {install.packages('devtools')}
        library('devtools')

        #if(!require('sydneycytometry/spectre')) {install_github("sydneycytometry/spectre")}
        library("Spectre")

        # We recommend not to update packages that are dependencies of Spectre


    ### 1.2. Install packages

        ## Basic data manipulation
        if(!require('plyr')) {install.packages('plyr')}
        if(!require('data.table')) {install.packages('data.table')}
        if(!require('tidyr')) {install.packages('tidyr')}
        if(!require('rstudioapi')) {install.packages('rstudioapi')}

        ## Required to download packages from Bioconductor
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")

        ## Basic data manipulation
        if(!require('flowCore')) {BiocManager::install('flowCore')}
        if(!require('flowViz')) {BiocManager::install('flowViz')}
        if(!require('Biobase')) {BiocManager::install('Biobase')}

        ## Clustering and dimensionality reduction
        if(!require('FlowSOM')) {BiocManager::install('FlowSOM')}
        if(!require('Rtsne')) {install.packages("Rtsne")} # for running tSNE
        if(!require('umap')) {install.packages('umap')}

        ## Plotting
        if (!require("ggplot2")){install.packages("ggplot2")} # for plotting tSNE graphs
        if (!require("colorRamps")){install.packages("colorRamps")} # for colour scheme management
        if (!require("ggthemes")){install.packages("ggthemes")} # for plot themes
        if (!require("scales")){install.packages("scales")} # for re-scaling if necessary
        if (!require("RColorBrewer")){install.packages("RColorBrewer")} # for re-scaling if necessary
        if (!require("gridExtra")){install.packages("gridExtra")} # for re-scaling if necessary


    ### 1.3. Load packages from library
        library('plyr')
        library('data.table')
        library('tidyr') # for spread

        library('rstudioapi')
        library('flowViz')
        library('flowCore')
        library('Biobase')
        library('FlowSOM')
        library('Rtsne')
        library('umap')
        library('ggplot2')
        library('scales')
        library('colorRamps')
        library('ggthemes')
        library('RColorBrewer')
        library("gridExtra")

    ## 1.4. Set working directory

        ## Set working directory
        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        getwd()
        PrimaryDirectory <- getwd()
        PrimaryDirectory

        ## Can set manually using these lines, if desired
            #PrimaryDirectory <- "/Users/thomasashhurst/Documents/Github/Public Github/Spectre/Other/Demo_dataset/"
            #setwd(PrimaryDirectory)

        ## Create output directory
        dir.create("Output_CAPX", showWarnings = FALSE)
        setwd("Output_CAPX")
        OutputDirectory <- getwd()
        setwd(PrimaryDirectory)

##########################################################################################################
#### 2. Read and prepare data
##########################################################################################################

    ### Define an experiment name/serial number

        exp.name <- "Demo"

    ### Read samples into workspace and review

        ## List of CSV files in PrimaryDirectory
        list.files(PrimaryDirectory, ".csv")

            ## ADD A PRE-PROCESSING SCRIPT BEFORE THIS ONE -- FILE MERGE etc
            ## HERE WE WANT ONE FILE PER SAMPLE

        ## Import samples (read files into R from disk)
        Spectre::read.files(file.loc = PrimaryDirectory, file.type = ".csv", do.embed.file.names = TRUE)

        ncol.check    # Review number of columns (features, markers) in each sample
        nrow.check    # Review number of rows (cells) in each sample

        ## Review column names and their subsequent values
        name.table
        head(data.list)
        head(data.list[[1]])

        ## Save starting data
        data.start <- data.list


        ##################################################################
        #### READ IN SAMPLE TABLE FOR ALLOCATIONS
        ##################################################################


        sample.table <- read.delim(file = "sample.table.txt")

        file.col <- "Filename"
        sample.col <- "Sample"
        group.col <- "Group"
        batch.col <- "Batch"


        all.sample.names <- sample.table[sample.col]
        all.group.names <- sample.table[group.col]
        all.batch.names <- sample.table[batch.col]
        cells.per.sample <- sample.table["Cells.per.sample"]


        all.file.names == unique(names(data.list))

            #### USE SOME KIND OF LOOP BASED ON WHATS IN THE TABLE -- ADD KEYWORDS/NUMS to the samples

                # make lists from what's in the table

                for(i in c(1:length(all.file.names))){
                  #i <- 1
                    data.list[[i]][[sample.col]] <- NA # fills a new colum
                    data.list[[i]][[sample.col]] <- all.sample.names[i,]

                    data.list[[i]][[group.col]] <- NA # fills a new colum
                    data.list[[i]][[group.col]] <- all.group.names[i,]

                    data.list[[i]][[batch.col]] <- NA # fills a new colum
                    data.list[[i]][[batch.col]] <- all.batch.names[i,]
                }


        head(data.list[[1]])
        head(data.list[[6]])
        head(data.list[[12]])

    ### Merge files

        ## Merge files and review
        Spectre::merge.files(x = data.list)

        head(cell.dat)
        dim(cell.dat)
        as.matrix(unique(cell.dat[["Sample"]]))

        ## Are there any NAs present in cell.dat? Yes if 'TRUE', no if 'FALSE'
        any(is.na(cell.dat))

        ## Cleanup (not necessary, but recommended)
        rm(data.list, data.start, name.table, ncol.check, nrow.check, all.file.names, all.file.nums)


    ################
        #dim(cell.dat)
        #cell.dat.large <- rbind(cell.dat, cell.dat)
        #for(i in c(1:8)){
        #  cell.dat.large <- rbind(cell.dat.large, cell.dat) # x8 or so
        #}
        #cell.dat <- cell.dat.large
    ################



##########################################################################################################
#### Define data and sample variables for analysis
##########################################################################################################


    ### Data subsampling?

        #Spectre::subsample(x = cell.dat,
        #                   method = "per.sample", # or "random
        #                   samp.col = "SampleName",
        #                   targets = c(10000),
        #                   seed = 42)

        #cell.dat.sub <- subsample.res
        #rm(subsample.res)

    ### Define cellular and clustering columns

        ## Save column names
        ColumnNames <- unname(colnames(cell.dat)) # assign reporter and marker names (column names) to 'ColumnNames'

        ## Define columns that are 'valid' cellular markers (i.e. not live/dead, blank channels etc)
        as.matrix(ColumnNames) # view the column 'number' for each parameter
        ValidCellularColsNos <- c(5,6,8,9,11:13,16:19,21:30,32)
        ValidCellularCols <- ColumnNames[ValidCellularColsNos]

        ValidCellularCols  # check that the column names that appear are the ones you want to analyse
        ColumnNames[-ValidCellularColsNos] # Check which columns are being EXCLUDED!
        rm(ValidCellularColsNos)


        ## Define columns for clustering
        as.matrix(ColumnNames)
        ClusteringColNos <- c(5,6,8,9,11,13,17:19,21:29,32)
        ClusteringCols <- ColumnNames[ClusteringColNos] # e.g. [c(11, 23, 10)] to include the markers corresponding to the column numbers 11, 23, 10

        ClusteringCols  # check that the column names that appear are the ones you want to analyse
        ColumnNames[-ClusteringColNos] # Check which columns are being EXCLUDED!
        rm(ClusteringColNos)


        ## Review key data
        head(cell.dat)
        sample.table
        ValidCellularCols
        ClusteringCols

        nrow(cell.dat)

        ## Remove unused parameters



##########################################################################################################
#### Perform clustering
##########################################################################################################

    ### Run FlowSOM
        Spectre::run.flowsom(x = cell.dat,
                             meta.k = 40,
                             clustering.cols = ClusteringCols,
                             clust.seed = 42,
                             meta.seed = 42,
                             clust.name = "FlowSOM_cluster",
                             meta.clust.name = "FlowSOM_metacluster")

        cell.dat <- cbind(cell.dat, flowsom.res.original)   # Add results to cell.dat
        cell.dat <- cbind(cell.dat, flowsom.res.meta)       # Add results to cell.dat
        head(cell.dat)                                      # Check cell.dat to ensure FlowSOM data correctly attached

        rm(flowsom.res.original)                            # Remove results from global environment
        rm(flowsom.res.meta)                                # Remove results from global environment


    ### Perform other clustering approaches if desired

##########################################################################################################
#### Perform downsampling and dimensionality reduction
##########################################################################################################

    ### Subsampling

        ## Review cells per sample
        sample.table

        ## RUn subsampling
        Spectre::subsample(x = cell.dat,
                           method = "per.sample", # or "random
                           samp.col = sample.col,
                           targets = c(rep(100,12)),
                           seed = 42)

        cell.dat.sub <- subsample.res

        nrow(cell.dat.sub)
        rm(subsample.res)

        ######################
        #cell.dat.sub <- cell.dat.large
        ######################

    ### Run UMAP

        Spectre::run.umap(x = cell.dat.sub,
                          use.cols = ClusteringCols,
                          umap.seed = 42)

        cell.dat.sub <- cbind(cell.dat.sub, umap.res) # Merge UMAP results with data

        plot(cell.dat.sub$UMAP_42_X, cell.dat.sub$UMAP_42_Y)


    ### Run tSNE

        #Spectre::run.tsne(x = cell.dat.sub,
        #                  use.cols = ClusteringCols,  #c(5,6,8,9,11,12,13,17:19,21:30,32),
        #                  umap.seed = 42)
        #
        #cell.dat.sub <- cbind(cell.dat.sub, tsne.res) # Merge UMAP results with data

        ## Plot tSNE results
        #plot(cell.dat.sub$tSNE_42_X, cell.dat.sub$tSNE_42_Y)


##########################################################################################################
#### Basic plotting and data examination
##########################################################################################################

    ### Plot single UMAP MFI plot
        names(cell.dat.sub)

        p1 <- Spectre::colour.plot(d = cell.dat.sub,
                             x.axis = "UMAP_42_X",
                             y.axis = "UMAP_42_Y",
                             col.axis = "BV605.Ly6C",
                             title = paste0("All samples", " - ", "BV605.Ly6C"),
                             colours = "spectral",
                             dot.size = 1)
        p1

    ### Plot single UMAP factor plot -- clusters
        names(cell.dat.sub)

        p2 <- Spectre::factor.plot(d = cell.dat.sub,
                                   x.axis = "UMAP_42_X",
                                   y.axis = "UMAP_42_Y",
                                   col.axis = "FlowSOM_metacluster",
                                   title = paste0("All samples", " - ", "Clusters"),
                                   dot.size = 1,
                                   align.xy.by = cell.dat.sub,
                                   align.col.by = cell.dat)

        p2 # TRUE is fine, FALSE also returns 'NULL' when running this line IF using option for centroid.labels as FALSE

        p2 <- Spectre::labelled.factor.plot(d = cell.dat.sub,
                                            x.axis = "UMAP_42_X",
                                            y.axis = "UMAP_42_Y",
                                            col.axis = "FlowSOM_metacluster",
                                            title = paste0("All samples", " - ", "Clusters"),
                                            dot.size = 1,
                                            align.xy.by = cell.dat.sub,
                                            align.col.by = cell.dat)
        p2


    ### Plot single UMAP factor plot -- samples
        names(cell.dat.sub)

        p3 <- Spectre::factor.plot(d = cell.dat.sub,
                                   x.axis = "UMAP_42_X",
                                   y.axis = "UMAP_42_Y",
                                   col.axis = "Sample",
                                   title = "Samples",
                                   dot.size = 0.5) # assumes numeric

        p3


    ### Plot single UMAP factor plot -- groups
        names(cell.dat.sub)

        p4 <- Spectre::factor.plot(d = cell.dat.sub,
                                   x.axis = "UMAP_42_X",
                                   y.axis = "UMAP_42_Y",
                                   col.axis = "Group",
                                   title = "Groups",
                                   dot.size = 0.5) # assumes numeric

        p4

        #ggsave(filename = "p1.jpeg", plot = p1, path = OutputDirectory, width = 9, height = 7)
        #ggsave(filename = "p2.jpeg", plot = p2, path = OutputDirectory, width = 9, height = 7)
        #ggsave(filename = "p3.jpeg", plot = p3, path = OutputDirectory, width = 9, height = 7)
        #ggsave(filename = "p4.jpeg", plot = p4, path = OutputDirectory, width = 9, height = 7)


    ### Create tiles and save

        gp <- grid.arrange(grobs = list(p1, p2, p3, p4), ncol=2, nrow=2) #top = "Main Title"
        gp

        ggsave(filename = "Grid.jpeg", plot = gp, path = OutputDirectory, width = 18, height = 14)

        ### --> FUNCTION WITH XMAX/MIN etc
        ### --> seperate function for COL MAX/MIN etc


##########################################################################################################
#### Save summary statistics and data to disk
##########################################################################################################

    ### Create and save sumtables
        setwd(PrimaryDirectory)
        setwd(OutputDirectory)
        dir.create("Output-SumTables", showWarnings = FALSE)
        setwd("Output-SumTables")
        getwd()

        as.matrix(names(cell.dat))

        Spectre::make.sumtable(x = cell.dat,
                               type = "frequencies",
                               sample.name = "Sample",
                               group.name = "Group",
                               clust.col = "FlowSOM_metacluster",
                               annot.col.nums = c(1,33:39),
                               cells.per.tissue = c(rep(2e+07, 6), rep(1.8e+07, 6))
                               #do.foldchange = TRUE # not active yet
                               #ctrl.group = "Mock"
                               )

        Spectre::make.sumtable(x = cell.dat,
                               type = "expression.per.sample",
                               sample.name = "Sample",
                               group.name = "Group",
                               clust.col = "FlowSOM_metacluster",
                               annot.col.nums = c(1,33:39),
                               fun.type = "median"
                               #do.foldchange = TRUE # not active yet
                               #ctrl.group = "Mock"
                               )

        Spectre::make.sumtable(x = cell.dat,
                               type = "expression.per.marker",
                               sample.name = "Sample",
                               group.name = "Group",
                               clust.col = "FlowSOM_metacluster",
                               annot.col.nums = c(1,33:39),
                               fun.type = "median"
                               #do.foldchange = TRUE # not active yet
                               #ctrl.group = "Mock"
                               )

        ## didn't exclude V1? annot col maybe didn't work
        ## Also some 'clusters' keeping the 'CLUSTER' label
        ## Create a text file in each output folder -- explaining how to read the results

    ### Save data (cell.dat) including clustering results

        setwd(PrimaryDirectory)
        setwd(OutputDirectory)
        dir.create("Output-ClusteredData", showWarnings = FALSE)
        setwd("Output-ClusteredData")
        getwd()

        head(cell.dat)

        ## Write 'all' data
        Spectre::write.files(x = cell.dat,
                             file.prefix = exp.name, # required
                             write.csv = TRUE,
                             write.fcs = TRUE)

        ## Write 'by sample' data
        Spectre::write.files(x = cell.dat,
                             file.prefix = exp.name, # required
                             divide.by = "Sample",
                             write.csv = TRUE,
                             write.fcs = TRUE)

        ## Write 'by group' data
        Spectre::write.files(x = cell.dat,
                             file.prefix = exp.name, # required
                             divide.by = "Group",
                             write.csv = TRUE,
                             write.fcs = TRUE)


        setwd(PrimaryDirectory)

    ### Save subsampled data (cell.dat.sub) includng tSNE/UMAP etc

        setwd(PrimaryDirectory)
        setwd(OutputDirectory)
        dir.create("Output-DimRedData", showWarnings = FALSE)
        setwd("Output-DimRedData")
        getwd()

        head(cell.dat.sub)

        ## Write 'all' data
        Spectre::write.files(x = cell.dat.sub,
                             file.prefix = paste0("DimRed_", exp.name), # required
                             write.csv = TRUE,
                             write.fcs = TRUE) ##### FCS NOT WORK

        ## Write 'by sample' data
        Spectre::write.files(x = cell.dat.sub,
                             file.prefix = paste0("DimRed_", exp.name), # required
                             divide.by = "Sample",
                             write.csv = TRUE,
                             write.fcs = TRUE) ##### FCS NOT WORK

        ## Write 'by group' data
        Spectre::write.files(x = cell.dat.sub,
                             file.prefix = paste0("DimRed_", exp.name), # required
                             divide.by = "Group",
                             write.csv = TRUE,
                             write.fcs = TRUE) ##### FCS NOT WORK

        setwd(PrimaryDirectory)
#
#
#
#
#
#         ### Define sample and group names
#
#         ## Define the sample and group column names
#         names(cell.dat)
#
#         samp.col <- "Sample"
#         grp.col <- "Group"
#         batch.col <- "Batch"
#
#         ## Check to ensure the correct name has been specified above
#         cell.dat[[samp.col]]
#         cell.dat[[grp.col]]
#         cell.dat[[batch.col]]
#
#
#
#
#         ## Make a sample table
#         make.sample.table(x = cell.dat,
#                           sample.col.name = samp.col,
#                           include.groups = TRUE,
#                           group.col.name = grp.col,
#                           include.batch = TRUE,
#                           batch.col.name = batch.col
#         )
#
#         # Check results
#         all.sample.names
#         all.group.names
#         all.batch.names
#
#         sample.table
