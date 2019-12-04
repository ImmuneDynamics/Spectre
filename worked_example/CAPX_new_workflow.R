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

        if(!require('Spectre')) {install_github("sydneycytometry/spectre")}
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
        if(!require('Biobase')) {BiocManager::install('Biobase')}
        if(!require('flowViz')) {BiocManager::install('flowViz')}

        ## Clustering and dimensionality reduction
        if(!require('FlowSOM')) {BiocManager::install('FlowSOM')}
        if(!require('Rtsne')) {install.packages("Rtsne")} # for running tSNE
        if(!require('umap')) {install.packages('umap')}

        ## Plotting and graphing
        if(!require("ggplot2")){install.packages("ggplot2")} # for plotting tSNE graphs
        if(!require("colorRamps")){install.packages("colorRamps")} # for colour scheme management
        if(!require("ggthemes")){install.packages("ggthemes")} # for plot themes
        if(!require("scales")){install.packages("scales")} # for re-scaling if necessary
        if(!require("RColorBrewer")){install.packages("RColorBrewer")} # for re-scaling if necessary
        if(!require("gridExtra")){install.packages("gridExtra")} # for re-scaling if necessary

    ### 1.3. Load packages from library
        library('plyr')
        library('data.table')
        library('tidyr') # for spread
        library('rstudioapi')

        library('flowCore')
        library('Biobase')
        library('flowViz')

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

    ### Create some metadata
        meta.dat <- list()

        list.files(PrimaryDirectory, ".txt")

        meta.dat[["expDetails"]] <- read.delim(file = "experiment.details.txt")
        meta.dat[["sampleDetails"]] <- read.delim(file = "sample.details.txt")

        # ## Warning message:
        # In read.table(file = file, header = header, sep = sep, quote = quote,  :
        #                 incomplete final line found by readTableHeader on 'experiment.details.txt'

    ### Read SAMPLES (data) into workspace and review

        ## List of CSV files in PrimaryDirectory ## ADD A PRE-PROCESSING SCRIPT BEFORE THIS ONE -- FILE MERGE etc             ## HERE WE WANT ONE FILE PER SAMPLE
        list.files(PrimaryDirectory, ".csv")

        ## Import samples (read files into R from disk)
        Spectre::read.files(file.loc = PrimaryDirectory,
                            file.type = ".csv",
                            do.embed.file.names = TRUE)

        ncol.check    # Review number of columns (features, markers) in each sample
        nrow.check    # Review number of rows (cells) in each sample
        name.table    # Review column names and their subsequent values

        ## Check data
        head(data.list)
        head(data.list[[1]])

        ## Save starting data
        data.start <- data.list

    ### Read in sample METADATA (NOT REQUIRED FOR REST OF CAPX SCRIPT)

        ## Specify the column that contains filenames, and which columns of 'sample.table' you want to embed
        meta.dat$sampleDetails

        # to.embed <- meta.dat$sampleDetails[c(1:4)] ## INCLUDING FILENAME
        # file.col <- "Filename"

        Spectre::embed.columns(x = data.list,
                               type = "list",
                               match.to = meta.dat$sampleDetails[c(1)],
                               new.cols = meta.dat$sampleDetails[c(2)],
                               col.name = names(meta.dat$sampleDetails[c(2)]))

        Spectre::embed.columns(x = data.list,
                               type = "list",
                               match.to = meta.dat$sampleDetails[c(1)],
                               new.cols = meta.dat$sampleDetails[c(3)],
                               col.name = names(meta.dat$sampleDetails[c(3)]))

        Spectre::embed.columns(x = data.list,
                               type = "list",
                               match.to = meta.dat$sampleDetails[c(1)],
                               new.cols = meta.dat$sampleDetails[c(4)],
                               col.name = names(meta.dat$sampleDetails[c(4)]))

        head(data.list)

    ### Merge files

        ## Merge files and review
        Spectre::file.merge(x = data.list)

            ### --> can use rbindlist from data.table instead of the plyr option

        str(cell.dat)
        head(cell.dat)
        dim(cell.dat)
        as.matrix(unique(cell.dat[["Sample"]]))

        ## Are there any NAs present in cell.dat? Yes if 'TRUE', no if 'FALSE'
        any(is.na(cell.dat))

        ## Cleanup (not necessary, but recommended)
        rm(data.list, data.start, ncol.check, nrow.check, all.file.names, all.file.nums)

##########################################################################################################
#### Define data and sample variables for analysis
##########################################################################################################

    ### Define key columns

        as.matrix(names(cell.dat))

        ## Define key columns that might be used or dividing data (samples, groups, batches, etc)
        exp.name <- meta.dat[["expDetails"]]$Experiment.Number

        file.col <- "Filename"
        sample.col <- "Sample"
        group.col <- "Group"
        batch.col <- "Batch"

        ## Create a list of column names
        ColumnNames <- as.matrix(unname(colnames(cell.dat))) # assign reporter and marker names (column names) to 'ColumnNames'
        ColumnNames

    ### Define cellular and clustering columns

        ## Define columns that are 'valid' cellular markers (i.e. not live/dead, blank channels etc)
        ColumnNames # view the column 'number' for each parameter

        ValidCellularColsNos <- c(5,6,8,9,11:13,16:19,21:30,32)
        ValidCellularCols <- ColumnNames[ValidCellularColsNos]

        ValidCellularCols  # check that the column names that appear are the ones you want to analyse
        ColumnNames[-ValidCellularColsNos] # Check which columns are being EXCLUDED!

        meta.dat[["ValidCellularCols"]] <- ValidCellularCols
        rm(ValidCellularColsNos)

    ### Define columns for clustering

        ## Define columns that are 'valid' cellular markers (i.e. not live/dead, blank channels etc)
        ColumnNames
        ClusteringColNos <- c(5,6,8,9,11,13,17:19,21:29,32)
        ClusteringCols <- ColumnNames[ClusteringColNos] # e.g. [c(11, 23, 10)] to include the markers corresponding to the column numbers 11, 23, 10

        ClusteringCols  # check that the column names that appear are the ones you want to analyse
        ColumnNames[-ClusteringColNos] # Check which columns are being EXCLUDED!

        meta.dat[["ClusteringCols"]] <- ClusteringCols
        rm(ClusteringColNos)

    ### Check

        head(cell.dat)
        meta.dat

        # TEST IF THE KEY METADATA EXISTS

    ### Create a blank log dataframe in meta.data
        meta.dat[["analysis.log"]] <- data.frame("Function" = NA, "Parameter" = NA, "Value" = NA)

##########################################################################################################
#### Perform clustering
##########################################################################################################

    ### Run FlowSOM
        Spectre::run.flowsom(x = cell.dat,
                             meta.k = 40,
                             clustering.cols = meta.dat$ClusteringCols,
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
        meta.dat$sampleDetails

        Spectre::subsample(x = cell.dat,
                           method = "per.sample", # or "random
                           samp.col = sample.col,
                           targets = c(rep(100,12)),
                           seed = 42)

        cell.dat.sub <- subsample.res

        nrow(cell.dat.sub)
        rm(subsample.res)

    ### Run UMAP
        Spectre::run.umap(x = cell.dat.sub,
                          use.cols = meta.dat$ClusteringCols,
                          umap.seed = 42)

        cell.dat.sub <- cbind(cell.dat.sub, umap.res) # Merge UMAP results with data
        plot(cell.dat.sub$UMAP_42_X, cell.dat.sub$UMAP_42_Y)

    ### Run tSNE
        # Spectre::run.tsne(x = cell.dat.sub)

    ### Run Monocle
        # Spectre::run.monocle(x = cell.dat.sub)


    ### Review results

        p1 <- Spectre::colour.plot(d = cell.dat.sub,
                                   x.axis = "UMAP_42_X",
                                   y.axis = "UMAP_42_Y",
                                   col.axis = "BV605.Ly6C",
                                   title = paste0("All samples", " - ", "BV605.Ly6C"),
                                   align.xy.by = cell.dat.sub,
                                   align.col.by = cell.dat.sub
                                   )
        p1


##########################################################################################################
#### Save data to disk
##########################################################################################################

    ### Save data (cell.dat) including clustering results

        setwd(OutputDirectory)
        head(cell.dat)

        ## Write 'all' data
        Spectre::write.files(x = cell.dat,
                             file.prefix= paste0("Clustered_", exp.name), # required
                             write.csv = TRUE,
                             write.fcs = TRUE)

        # Spectre::write.files(x = cell.dat,
        #                      divide.by = "Group", ## Will add the terms as a prefix
        #                      file.prefix= paste0("Clustered", exp.name), # required
        #                      write.csv = FALSE,
        #                      write.fcs = TRUE)
        #
        # Spectre::write.files(x = cell.dat,
        #                      divide.by = "Sample",
        #                      file.prefix= paste0("Clustered", exp.name), # required
        #                      write.csv = FALSE,
        #                      write.fcs = TRUE)

        setwd(PrimaryDirectory)

    ### Save subsampled data (cell.dat.sub) includng tSNE/UMAP etc

        setwd(OutputDirectory)
        getwd()

        head(cell.dat.sub)

        ## Write 'all' data
        Spectre::write.files(x = cell.dat.sub,
                             file.prefix = paste0("DimRed_", exp.name), # required
                             write.csv = TRUE,
                             write.fcs = TRUE)

        # Spectre::write.files(x = cell.dat,
        #                      divide.by = "Group", ## Will add the terms as a prefix
        #                      file.prefix= paste0("DimRed_", exp.name), # required
        #                      write.csv = FALSE,
        #                      write.fcs = TRUE)
        #
        # Spectre::write.files(x = cell.dat,
        #                      divide.by = "Sample",
        #                      file.prefix= paste0("DimRed_", exp.name), # required
        #                      write.csv = FALSE,
        #                      write.fcs = TRUE)

        setwd(PrimaryDirectory)


##########################################################################################################
#### Save summary statistics to disk
##########################################################################################################

    ### Create and save sumtables
        setwd(PrimaryDirectory)
        setwd(OutputDirectory)
        dir.create("Output-SumTables", showWarnings = FALSE)
        setwd("Output-SumTables")
        getwd()

        as.matrix(names(cell.dat))
        meta.dat$sampleDetails

        Spectre::make.sumtable(x = cell.dat,
                               sample.col = "Sample",
                               clust.col = "FlowSOM_metacluster",
                               annot.col.nums = c(1,33:39),

                               do.frequencies = TRUE,
                               do.exp.per.marker = TRUE,
                               do.exp.per.sample = TRUE,

                               cells.per.tissue = meta.dat$sampleDetails$Cells.per.sample, ## MODIFY TO DO MATCHING
                               fun.type = "median",

                               do.foldchange = TRUE,
                               group.name = "Group",
                               control.group = "Mock"
                              )

        setwd(PrimaryDirectory)



##########################################################################################################
#### Print tSNE/UMAP plots to disk
##########################################################################################################










##########################################################################################################
##########################################################################################################
##########################################################################################################







        ##

        # Spectre::make.sumtable(x = cell.dat,
        #                        type = "frequencies",
        #                        sample.name = "Sample",
        #                        group.name = "Group",
        #                        clust.col = "FlowSOM_metacluster",
        #                        annot.col.nums = c(1,33:39),
        #                        cells.per.tissue = meta.dat$sampleDetails$Cells.per.sample
        #                        #do.foldchange = TRUE # not active yet
        #                        #ctrl.group = "Mock"
        # )
        #
        #
        # Spectre::make.sumtable(x = cell.dat,
        #                        type = "expression.per.sample",
        #                        sample.name = "Sample",
        #                        group.name = "Group",
        #                        clust.col = "FlowSOM_metacluster",
        #                        annot.col.nums = c(1,33:39),
        #                        fun.type = "median"
        #                        #do.foldchange = TRUE # not active yet
        #                        #ctrl.group = "Mock"
        #                        )
        #
        # Spectre::make.sumtable(x = cell.dat,
        #                        type = "expression.per.marker",
        #                        sample.name = "Sample",
        #                        group.name = "Group",
        #                        clust.col = "FlowSOM_metacluster",
        #                        annot.col.nums = c(1,33:39),
        #                        fun.type = "median"
        #                        #do.foldchange = TRUE # not active yet
        #                        #ctrl.group = "Mock"
        #                        )

        ## didn't exclude V1? annot col maybe didn't work
        ## Also some 'clusters' keeping the 'CLUSTER' label
        ## Create a text file in each output folder -- explaining how to read the results




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
