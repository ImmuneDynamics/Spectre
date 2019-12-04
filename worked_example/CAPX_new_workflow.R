##########################################################################################################
#### DRAFT Cytometry Analysis Pipeline for large and compleX data (CAPX) v3.0 - using the Spectre R package
#### Part 1 - Clustering, dimensionality reduction, save files
##########################################################################################################

    # Thomas Myles Ashhurst, Felix Marsh-Wakefield, Givanna Putri
    # 2019-12-02
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

    ### 1.2. Install other packages

        ## Required to download packages from Bioconductor
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")

        ## Basic data manipulation
        if(!require('flowCore')) {BiocManager::install('flowCore')}
        if(!require('Biobase')) {BiocManager::install('Biobase')}
        if(!require('flowViz')) {BiocManager::install('flowViz')}

        ## Clustering and dimensionality reduction
        if(!require('FlowSOM')) {BiocManager::install('FlowSOM')}

            ## Plotting and graphing
            if(!require("colorRamps")){install.packages("colorRamps")} # for colour scheme management
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

    ### 1.4. Set working directory

        ## Set working directory
        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located

        # setwd("/Users/Tom/Google Drive (t.ashhurst@centenary.org.au)/_Sydney Cytometry/2019_Synced/GitHub/Public github/Spectre - workflow scripts/")

        getwd()
        PrimaryDirectory <- getwd()
        PrimaryDirectory


        # setwd("../../Spectre - workflow scripts/")
        # getwd()
        # PrimaryDirectory <- getwd()
        # PrimaryDirectory

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

    ### Read sample metadata and embed in sample data
        meta.dat <- read.delim(file = "sample.details.txt")

        Spectre::embed.columns(x = data.list,
                               type = "list",
                               match.to = meta.dat[c(1)],
                               new.cols = meta.dat[c(2)],
                               col.name = names(meta.dat[c(2)]))

        Spectre::embed.columns(x = data.list,
                               type = "list",
                               match.to = meta.dat[c(1)],
                               new.cols = meta.dat[c(3)],
                               col.name = names(meta.dat[c(3)]))

        Spectre::embed.columns(x = data.list,
                               type = "list",
                               match.to = meta.dat[c(1)],
                               new.cols = meta.dat[c(4)],
                               col.name = names(meta.dat[c(4)]))


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
        exp.name <- "TAdemo"

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

        CellularColsNos <- c(5,6,8,9,11:13,16:19,21:30,32)
        CellularCols <- ColumnNames[CellularColsNos]

        CellularCols  # check that the column names that appear are the ones you want to analyse
        ColumnNames[-CellularColsNos] # Check which columns are being EXCLUDED!

    ### Define columns for clustering

        ## Define columns that are 'valid' cellular markers (i.e. not live/dead, blank channels etc)
        ColumnNames
        ClusteringColNos <- c(5,6,8,9,11,13,17:19,21:29,32)
        ClusteringCols <- ColumnNames[ClusteringColNos] # e.g. [c(11, 23, 10)] to include the markers corresponding to the column numbers 11, 23, 10

        ClusteringCols  # check that the column names that appear are the ones you want to analyse
        ColumnNames[-ClusteringColNos] # Check which columns are being EXCLUDED!

    ### Check

        head(cell.dat)
        CellularCols
        ClusteringCols
        meta.dat


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
        # meta.dat
        #
        # as.matrix(unique(cell.dat["Sample"]))
        #
        # Spectre::subsample(x = cell.dat,
        #                    method = "per.sample", # or "random
        #                    samp.col = sample.col,
        #                    targets = c(rep(100,12)),
        #                    seed = 42)
        #
        # cell.dat.sub <- subsample.res
        # nrow(cell.dat.sub)
        #
        # rm(subsample.res)

    ### Run UMAP

        cell.dat.sub <- cell.dat

        Spectre::run.umap(x = cell.dat.sub,
                          use.cols = ClusteringCols,
                          umap.seed = 42)

        cell.dat.sub <- cbind(cell.dat.sub, umap.res) # Merge UMAP results with data
        plot(cell.dat.sub$UMAP_42_X, cell.dat.sub$UMAP_42_Y)

    ### Review results
        Spectre::colour.plot(d = cell.dat.sub,
                             x.axis = "UMAP_42_X",
                             y.axis = "UMAP_42_Y",
                             col.axis = "BV605.Ly6C",
                             title = paste0("All samples", " - ", "BV605.Ly6C"),
                             align.xy.by = cell.dat.sub,
                             align.col.by = cell.dat.sub,
                             colours = "inferno"
                             )


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

        ## Write 'subsample' data
        Spectre::write.files(x = cell.dat.sub,
                             file.prefix = paste0("DimRed_", exp.name), # required
                             write.csv = TRUE,
                             write.fcs = TRUE)

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

        meta.dat
        as.matrix(names(cell.dat))

        Spectre::make.sumtable(x = cell.dat,
                               sample.col = "Sample",
                               clust.col = "FlowSOM_metacluster",
                               annot.col.nums = c(1,33:39),

                               do.frequencies = TRUE,
                               do.exp.per.marker = TRUE,
                               do.exp.per.sample = TRUE,

                               cells.per.tissue = meta.dat$Cells.per.sample, ## MODIFY TO DO MATCHING
                               fun.type = "median",

                               do.foldchange = TRUE,
                               group.col = "Group",
                               control.group = "Mock"
                              )

        setwd(PrimaryDirectory)


##########################################################################################################
#### Print tSNE/UMAP plots to disk
##########################################################################################################

    ### Loop for cellular markers etc

        setwd(PrimaryDirectory)
        setwd(OutputDirectory)
        dir.create("Output-ColourPlots")
        setwd("Output-ColourPlots")

        plots <- CellularCols

        ## Plot for all data
            for(a in plots){
              p <- Spectre::colour.plot(d = cell.dat.sub, x.axis = "UMAP_42_X", y.axis = "UMAP_42_Y",
                                        col.axis = a, title = a, colours = "spectral", dot.size = 1)
              ggsave(p, filename = paste0("All_samples_", a, ".png"), width = 9, height = 7)
            }

        ## Plot divided by groups

        ## Plot divided by samples
            ## Plot by sample
            #for(i in unique(cell.dat.sub[["FileName"]])){
            #  # subset
            #}
            ## Plot by group

    ### Loop for clusters etc

        head(cell.dat.sub)
        factors <- c("FlowSOM_metacluster")
        setwd(OutputDirectory)

        ## Plot for all data
            for(a in factors){
              p <- Spectre::factor.plot(d = cell.dat.sub,
                                        x.axis = "UMAP_42_X",
                                        y.axis = "UMAP_42_Y",
                                        col.axis = "FlowSOM_metacluster",
                                        title = "Cluster",
                                        dot.size = 1,
                                        add.labels = TRUE) # assumes numeric

              ggsave(p, filename = paste0("All_samples_", a, ".png"), width = 9, height = 7)
            }

        setwd(PrimaryDirectory)

    ### Loop for samples/groups/batches etc


