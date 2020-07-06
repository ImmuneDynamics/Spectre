##########################################################################################################
#### Spectre - General Discovery Workflow
#### Clustering, dimensionality reduction, plotting, and summarise data
##########################################################################################################

    # Spectre R package: https://sydneycytometry.org.au/spectre
    # Thomas Myles Ashhurst, Felix Marsh-Wakefield, Givanna Putri

##########################################################################################################
#### 1. Load packages, and set working directory
##########################################################################################################

    ### 1.1. Install 'Spectre' package (using devtools) and the dependencies that Spectre requires

        # For instructions on installing Spectre, please visit https://wiki.centenary.org.au/display/SPECTRE

    ### 1.2. Load packages

        library(Spectre)
        Spectre::package.check()    # Check that all required packages are installed
        Spectre::package.load()     # Load required packages

    ### 1.3. Set number of threads for data.table functions

        getDTthreads()

    ### 1.4. Set working directory

        ## Set PrimaryDirectory
        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        getwd()
        PrimaryDirectory <- getwd()
        PrimaryDirectory

        ## Set 'input' directory
        setwd(PrimaryDirectory)
        setwd("data/")
        InputDirectory <- getwd()
        InputDirectory
        setwd(PrimaryDirectory)

        ## Set 'metadata' directory
        setwd(PrimaryDirectory)
        setwd("metadata/")
        MetaDirectory <- getwd()
        MetaDirectory
        setwd(PrimaryDirectory)

        ## Create output directory
        dir.create("Output_Spectre", showWarnings = FALSE)
        setwd("Output_Spectre")
        OutputDirectory <- getwd()

    ### 1.5. Create "Output-info" directory and save session info
        dir.create("Output-info", showWarnings = FALSE)
        setwd("Output-info")

        sink(file = "session_info.txt", append=TRUE, split=FALSE, type = c("output", "message"))
        session_info()
        sink()

        setwd(PrimaryDirectory)

##########################################################################################################
#### 2. Read and prepare data
##########################################################################################################

    ### 2.1. Read SAMPLES (data) into workspace and review

        ## List of CSV files in InputDirectory
        setwd(InputDirectory)
        list.files(InputDirectory, ".csv")

        ## Import samples (read files into R from disk)
        data.list <- Spectre::read.files(file.loc = InputDirectory,
                                         file.type = ".csv",
                                         do.embed.file.names = TRUE)

        ## Some checks
        ncol.check    # Review number of columns (features, markers) in each sample
        nrow.check    # Review number of rows (cells) in each sample
        name.table    # Review column names and their subsequent values

        head(data.list)
        head(data.list[[1]])

    ### 2.2. Merge files

        ## Merge files and review
        cell.dat <- Spectre::do.merge.files(dat = data.list)

        ## Some checks
        str(cell.dat) # check the structure of cell.dat
        head(cell.dat) # check the first and last rows of cell dat
        dim(cell.dat) # check the dimensionality of cell.dat
        as.matrix(unique(cell.dat[["Sample"]])) # check the column names of cell.dat

        ## Are there any NAs present in cell.dat? Yes if 'TRUE', no if 'FALSE'
        any(is.na(cell.dat))

    ### 2.3. Read sample metadata and embed in sample data

        ## Read in metadata
        setwd(MetaDirectory)
        meta.dat <- fread("sample.details.csv")
        meta.dat
        setwd(PrimaryDirectory)

        names(meta.dat)
        meta.dat <- meta.dat[,c(1:4)]
        meta.dat

        cell.dat <- do.embed.columns(dat = cell.dat,
                                     base.col = "FileName",
                                     add.dat = meta.dat,
                                     add.by = "Filename")

        cell.dat

        ## Cleanup (not necessary, but recommended)
        rm(data.list, all.file.names, all.file.nums)

##########################################################################################################
#### 3. Define data and sample variables for analysis
##########################################################################################################

    ### Define key columns

        as.matrix(names(cell.dat))

        ## Define key columns that might be used or dividing data (samples, groups, batches, etc)
        exp.name <- "TAdemo"

        file.col <- "FileName"
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

    ### Define downsample targets

        as.matrix(unique(cell.dat[["Sample"]]))
        data.frame(table(cell.dat[["Sample"]])) # Check number of cells per sample.

        down.samp.targets <- c(rep(500,12))

    ### Checks

        head(cell.dat)
        CellularCols
        ClusteringCols
        meta.dat

    ### Save CellularCols and ClusteringCols

        setwd(OutputDirectory)
        write.csv(x = CellularCols, file = "Output-info/CellularCols.csv", row.names = FALSE)
        write.csv(x = ClusteringCols, file = "Output-info/ClusteringCols.csv", row.names = FALSE)
        write.csv(x = meta.dat, file = "Output-info/metadata.csv", row.names = FALSE)
        setwd(PrimaryDirectory)

##########################################################################################################
#### 4. Perform clustering and dimensionality reduction
##########################################################################################################

    ### Run FlowSOM
        cell.dat <- Spectre::run.flowsom(dat = cell.dat,
                                         xdim = 10,
                                         ydim = 10,
                                         meta.k = 10,
                                         use.cols = ClusteringCols)

        head(cell.dat)    # Check cell.dat to ensure FlowSOM data correctly attached

    ### Subsampling
        meta.dat
        as.matrix(unique(cell.dat[[sample.col]]))

        cell.dat.sub <- Spectre::do.subsample(dat = cell.dat,
                                              targets = down.samp.targets,
                                              divide.by = sample.col)

        nrow(cell.dat.sub)

    ### Run UMAP
        cell.dat.sub <- Spectre::run.umap(dat = cell.dat.sub,
                                          use.cols = ClusteringCols,
                                          umap.seed = 42)


    ### Preview results (without saving to disk)
        Spectre::make.colour.plot(dat = cell.dat.sub,
                                  x.axis = "UMAP_X",
                                  y.axis = "UMAP_Y",
                                  col.axis = "FlowSOM_metacluster",
                                  col.type = 'factor',
                                  add.label = TRUE,
                                  save.to.disk = FALSE)

##########################################################################################################
#### 5. Save data to disk
##########################################################################################################

    ### Save data (cell.dat) including clustering results
        setwd(OutputDirectory)
        dir.create("Output-data")
        setwd("Output-data")

        head(cell.dat)
        head(cell.dat.sub)

    ### Write 'large' dataset
        Spectre::write.files(dat = cell.dat,
                             file.prefix= paste0("Clustered_", exp.name), # required
                             write.csv = TRUE,
                             write.fcs = TRUE)

        Spectre::write.files(dat = cell.dat,
                             file.prefix= paste0("Clustered_", exp.name), # required
                             divide.by = "Sample",
                             write.csv = FALSE,
                             write.fcs = TRUE)

    ### Write 'subsample' dataset
        Spectre::write.files(dat = cell.dat.sub,
                             file.prefix = paste0("DimRed_", exp.name), # required
                             write.csv = TRUE,
                             write.fcs = TRUE)

        Spectre::write.files(dat = cell.dat.sub,
                             file.prefix = paste0("DimRed_", exp.name), # required
                             divide.by = "Sample",
                             write.csv = FALSE,
                             write.fcs = TRUE)

        setwd(PrimaryDirectory)

##########################################################################################################
#### Save expression heatmaps
##########################################################################################################

    ### Create and save expression pattern data (cluster x marker, per sample)
        setwd(OutputDirectory)
        dir.create("Output-expression-heatmaps")
        setwd("Output-expression-heatmaps")

        meta.dat
        as.matrix(names(cell.dat))

        Spectre::write.sumtables(dat = cell.dat,
                                 sample.col = sample.col,
                                 pop.col = "FlowSOM_metacluster",
                                 group.col = group.col,
                                 annot.col = "Batch",
                                 measure.col = CellularCols,
                                 do.proportions = FALSE,
                                 do.mfi.per.marker = FALSE,
                                 do.mfi.per.sample = TRUE)

        setwd(OutputDirectory)
        setwd("Output-expression-heatmaps")

    ### Create an expression heatmap

        ## Expression per cluster
        exp <- read.csv("SumTable-MFI-PerSample/SumTable-MFI-AllSamples.csv")

        library(pheatmap)
        make.pheatmap(dat = exp,
                      file.name = "MFI.heatmap.png",
                      plot.title = "Expression per cluster",
                      sample.col = "FlowSOM_metacluster",
                      plot.cols = CellularCols)

#########################################################################################################
#### Create and save some tSNE/UMAP plots for the whole dataset
#########################################################################################################

    ### Plot some 'all data' factor plots
        setwd(OutputDirectory)
        dir.create("Output-plots")
        setwd("Output-plots")

        Spectre::make.colour.plot(dat = cell.dat.sub,
                                  x.axis = "UMAP_X",
                                  y.axis = "UMAP_Y",
                                  col.axis = sample.col,
                                  col.type = 'factor')

        Spectre::make.colour.plot(dat = cell.dat.sub,
                                  x.axis = "UMAP_X",
                                  y.axis = "UMAP_Y",
                                  col.axis = group.col,
                                  col.type = 'factor')

        Spectre::make.colour.plot(dat = cell.dat.sub,
                                  x.axis = "UMAP_X",
                                  y.axis = "UMAP_Y",
                                  col.axis = batch.col,
                                  col.type = 'factor')

        Spectre::make.colour.plot(dat = cell.dat.sub,
                                  x.axis = "UMAP_X",
                                  y.axis = "UMAP_Y",
                                  col.axis = "FlowSOM_metacluster",
                                  col.type = 'factor',
                                  add.label = TRUE)

    ### Make some factor multiplots

        Spectre::make.multi.plot(dat = cell.dat.sub,
                                 col.type = "factor",
                                 x.axis = "UMAP_X",
                                 y.axis = "UMAP_Y",
                                 plot.by = sample.col,
                                 divide.by = group.col,
                                 figure.title = "Sample by Group")

        Spectre::make.multi.plot(dat = cell.dat.sub,
                                 col.type = "factor",
                                 x.axis = "UMAP_X",
                                 y.axis = "UMAP_Y",
                                 plot.by = group.col,
                                 divide.by = sample.col,
                                 figure.title = "Group by Sample")

        Spectre::make.multi.plot(dat = cell.dat.sub,
                                 col.type = "factor",
                                 x.axis = "UMAP_X",
                                 y.axis = "UMAP_Y",
                                 plot.by = group.col,
                                 divide.by = batch.col,
                                 figure.title = "Group by Batch")

        Spectre::make.multi.plot(dat = cell.dat.sub,
                                 col.type = "factor",
                                 x.axis = "UMAP_X",
                                 y.axis = "UMAP_Y",
                                 figure.title = "Multiplot - clusters by sample",
                                 plot.by = "FlowSOM_metacluster",
                                 divide.by = sample.col)

        Spectre::make.multi.plot(dat = cell.dat.sub,
                                 col.type = "factor",
                                 x.axis = "UMAP_X",
                                 y.axis = "UMAP_Y",
                                 figure.title = "Multiplot - clusters by group",
                                 plot.by = "FlowSOM_metacluster",
                                 divide.by = group.col)

    ### Make some expression multiplots

        Spectre::make.multi.plot(dat = cell.dat.sub,
                                 x.axis = "UMAP_X",
                                 y.axis = "UMAP_Y",
                                 plot.by = CellularCols,
                                 add.density = TRUE)

#########################################################################################################
#### Create tSNE/UMAP plots for each sample and group seperately
#########################################################################################################

    ### Plot some marker-oriented plots -- one set per group, X/Y/colour min and max set by the full dataset
        setwd(OutputDirectory)
        setwd("Output-plots")
        dir.create("by-group")
        setwd("by-group")

        for(i in as.matrix(unique(cell.dat.sub[[group.col]]))){
          temp <- cell.dat.sub[cell.dat.sub[[group.col]] == i,]

          Spectre::make.multi.plot(dat = temp,
                                   x.axis = "UMAP_X",
                                   y.axis = "UMAP_Y",
                                   plot.by = CellularCols,
                                   figure.title = paste0("Multi plot - group ", i, " - "),
                                   align.xy.by = cell.dat.sub,
                                   align.col.by = cell.dat.sub,
                                   add.density = TRUE)
        }


    ### Plot some marker-oriented plots -- one set per group, X/Y/colour min and max set by the full dataset
        setwd(OutputDirectory)
        setwd("Output-plots")
        dir.create("by-sample")
        setwd("by-sample")

        for(i in as.matrix(unique(cell.dat.sub[[sample.col]]))){
          temp <- cell.dat.sub[cell.dat.sub[[sample.col]] == i,]

          Spectre::make.multi.plot(dat = temp,
                                   x.axis = "UMAP_X",
                                   y.axis = "UMAP_Y",
                                   plot.by = CellularCols,
                                   figure.title = paste0("Multi plot - sample ", i, " - "),
                                   align.xy.by = cell.dat.sub,
                                   align.col.by = cell.dat.sub,
                                   add.density = TRUE)
        }

