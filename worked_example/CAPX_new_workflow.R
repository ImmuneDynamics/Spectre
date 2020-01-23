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

    ### 1.1. Load 'Spectre' package (using devtools) and the dependencies that Spectre requires
        if(!require('devtools')) {install.packages('devtools')}
        library('devtools')

        #install_github("sydneycytometry/spectre")
        install_github("sydneycytometry/spectre", ref = 'development') # option to install the development verison if required
        library("Spectre")

    ### 1.2. Install BioConductor packages

        ## Install BiocManager to download packages from Bioconductor
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")

        ## Download additional BioConductor packages
        if(!require('flowCore')) {BiocManager::install('flowCore')}
        if(!require('Biobase')) {BiocManager::install('Biobase')}
        if(!require('flowViz')) {BiocManager::install('flowViz')}
        if(!require('FlowSOM')) {BiocManager::install('FlowSOM')}

    ### 1.3. Load packages

        library(Spectre)
        Spectre::check.packages() # --> change so that message at the end is "All required packages have been successfully installed"
        Spectre::load.packages() # --> change so that message at the end is "All required packages have been successfully loaded"

        session_info()

    ### 1.4. Set number of threads for data.table functions

        getDTthreads()

    ### 1.4. Set working directory

        ## Set working directory
        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        getwd()
        PrimaryDirectory <- getwd()
        PrimaryDirectory

        ## Can set manually using these lines, if desired
            #PrimaryDirectory <- "/Users/Tom/Desktop/TAXXX"
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

        ## List of CSV files in PrimaryDirectory # HERE WE WANT ONE FILE PER SAMPLE
        list.files(PrimaryDirectory, ".csv")

        ## Import samples (read files into R from disk)
        Spectre::read.files(file.loc = PrimaryDirectory,
                            file.type = ".csv",
                            do.embed.file.names = TRUE)

        ## Some checks
        ncol.check    # Review number of columns (features, markers) in each sample
        nrow.check    # Review number of rows (cells) in each sample
        name.table    # Review column names and their subsequent values

        head(data.list)
        head(data.list[[1]])

        ## Save starting data
        data.start <- data.list

    ### Read sample metadata and embed in sample data
        meta.dat <- read.delim(file = "metadata/sample.details.txt")

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

        str(cell.dat)
        head(cell.dat)
        dim(cell.dat)
        as.matrix(unique(cell.dat[["Sample"]]))

        ## Are there any NAs present in cell.dat? Yes if 'TRUE', no if 'FALSE'
        any(is.na(cell.dat))

        ## Cleanup (not necessary, but recommended)
          #rm(data.list, data.start, ncol.check, nrow.check, all.file.names, all.file.nums)

##########################################################################################################
#### 3. Define data and sample variables for analysis
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

    ### Checks

        head(cell.dat)
        CellularCols
        ClusteringCols
        meta.dat


##########################################################################################################
#### 4. Perform clustering
##########################################################################################################

    ### Run FlowSOM
        Spectre::run.flowsom(x = cell.dat,
                             xdim = 10,
                             ydim = 10,
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

#########################################################################################################
#### 5. Perform downsampling and dimensionality reduction
##########################################################################################################

    ### Subsampling
        meta.dat
        as.matrix(unique(cell.dat[["Sample"]]))

        Spectre::subsample(x = cell.dat,
                           method = "per.sample", # or "random
                           samp.col = sample.col,
                           targets = c(rep(500,12)),
                           seed = 42)

        cell.dat.sub <- subsample.res
        nrow(cell.dat.sub)

        rm(subsample.res)

    ### Run UMAP

        Spectre::run.umap(x = cell.dat.sub,
                          use.cols = ClusteringCols,
                          umap.seed = 42)

        cell.dat.sub <- cbind(cell.dat.sub, umap.res) # Merge UMAP results with data
        plot(cell.dat.sub$UMAP_X, cell.dat.sub$UMAP_Y)

    ### Review results
        Spectre::colour.plot(d = cell.dat.sub,
                             x.axis = "UMAP_X",
                             y.axis = "UMAP_Y",
                             col.axis = "BV605.Ly6C",
                             title = paste0("All samples", " - ", "BV605.Ly6C"),
                             align.xy.by = cell.dat.sub,
                             align.col.by = cell.dat.sub,
                             colours = "magma"
                             )

    ### Generate some multi.plots
        setwd(OutputDirectory)

        ## Plot by 'Sample', coloured by 'Group'
        multi.plot(d = cell.dat.sub,
                   type = "factor",
                   x.axis = "UMAP_X",
                   y.axis = "UMAP_Y",
                   col.axis = group.col,
                   plot.by = sample.col,
                   align.xy.by = cell.dat.sub,
                   align.col.by = cell.dat.sub,
                   dot.size = 3
                   )

        ## Plot the whole dataset, one plot per marker in one image
        multi.marker.plot(d = cell.dat.sub,
                          x.axis = "UMAP_X",
                          y.axis = "UMAP_Y",
                          plot.by = c(CellularCols),
                          align.xy.by = cell.dat.sub,
                          align.col.by = cell.dat.sub,
                          colours = "magma",
                          figure.title = "Markers",
                          dot.size = 3)

        ## Plot by 'Sample', coloured by 'marker' levels
        for(i in CellularCols){
            multi.plot(d = cell.dat.sub,
                       type = "colour",
                       x.axis = "UMAP_X",
                       y.axis = "UMAP_Y",
                       col.axis = i,
                       colour = "magma",
                       plot.by = sample.col,
                       figure.title = i,
                       align.xy.by = cell.dat.sub,
                       align.col.by = cell.dat.sub,
                       dot.size = 3
                       )
            }


##########################################################################################################
#### 6. Save data to disk
##########################################################################################################

    ### Save data (cell.dat) including clustering results

        setwd(PrimaryDirectory)
        setwd(OutputDirectory)
        head(cell.dat)
        head(cell.dat.sub)

        ## Write 'large' dataset
        Spectre::write.files(x = cell.dat,
                             file.prefix= paste0("Clustered_", exp.name), # required
                             write.csv = TRUE,
                             write.fcs = TRUE)

        Spectre::write.files(x = cell.dat,
                             file.prefix= paste0("Clustered_", exp.name), # required
                             divide.by = "Sample",
                             write.csv = TRUE,
                             write.fcs = TRUE)


        ## Write 'subsample' dataset
        Spectre::write.files(x = cell.dat.sub,
                             file.prefix = paste0("DimRed_", exp.name), # required
                             write.csv = TRUE,
                             write.fcs = TRUE)

        Spectre::write.files(x = cell.dat.sub,
                             file.prefix = paste0("DimRed_", exp.name), # required
                             divide.by = "Sample",
                             write.csv = TRUE,
                             write.fcs = TRUE)

        setwd(PrimaryDirectory)

##########################################################################################################
#### 7. Save summary statistics to disk
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
                               cells.per.tissue = meta.dat$Cells.per.sample, ## CHECK THE ORDER OF OCCURANCE OF THE unique entries in the sample.col COLUMN -- the order or cell counts in the vector MUST be the same

                               do.exp.per.marker = TRUE,
                               do.exp.per.sample = TRUE,
                               fun.type = "median",

                               do.foldchange = TRUE,
                               group.col = "Group",
                               control.group = "Mock"
                               )

        setwd(PrimaryDirectory)
