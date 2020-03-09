##########################################################################################################
#### Spectre - General Discovery Workflow
#### Part 1/3 - Clustering, dimensionality reduction, plotting, and summarise data
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

    ### 2.2. Read sample metadata and embed in sample data

        ## Read in metadata
        setwd(MetaDirectory)
        meta.dat <- read.csv(file = "sample.details.csv")
        setwd(PrimaryDirectory)

        ### Embed sample metadata
        data.list <- Spectre::do.embed.columns(x = data.list,
                                           type = "list",
                                           match.to = meta.dat[c(1)],
                                           new.cols = meta.dat[c(2)],
                                           col.name = names(meta.dat[c(2)]))

        data.list <- Spectre::do.embed.columns(x = data.list,
                                           type = "list",
                                           match.to = meta.dat[c(1)],
                                           new.cols = meta.dat[c(3)],
                                           col.name = names(meta.dat[c(3)]))

        data.list <- Spectre::do.embed.columns(x = data.list,
                                           type = "list",
                                           match.to = meta.dat[c(1)],
                                           new.cols = meta.dat[c(4)],
                                           col.name = names(meta.dat[c(4)]))

        head(data.list)

    ### 2.3. Merge files

        ## Merge files and review
        cell.dat <- Spectre::do.merge.files(dat = data.list)

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
#### 4. Perform clustering
##########################################################################################################

    ### Run FlowSOM
        cell.dat <- Spectre::run.flowsom(dat = cell.dat,
                                         xdim = 10,
                                         ydim = 10,
                                         meta.k = 20,
                                         clustering.cols = ClusteringCols)

        head(cell.dat)    # Check cell.dat to ensure FlowSOM data correctly attached

#########################################################################################################
#### 5. Perform downsampling and dimensionality reduction
##########################################################################################################

    ### Subsampling
        meta.dat
        as.matrix(unique(cell.dat[["Sample"]]))

        cell.dat.sub <- Spectre::do.subsample(dat = cell.dat,
                                               method = "per.sample", # or "random
                                               samp.col = sample.col,
                                               targets = c(rep(500,12)),
                                               seed = 42)

        nrow(cell.dat.sub)

    ### Run UMAP
        cell.dat.sub <- Spectre::run.umap(dat = cell.dat.sub,
                                          use.cols = ClusteringCols,
                                          umap.seed = 42)


    ### Preview results (without saving to disk)
        Spectre::make.colour.plot(dat = cell.dat.sub,
                                  x.axis = "UMAP_X",
                                  y.axis = "UMAP_Y",
                                  col.axis = "BV605.Ly6C",
                                  align.xy.by = cell.dat.sub,
                                  align.col.by = cell.dat.sub,
                                  save.to.disk = FALSE)

##########################################################################################################
#### 6. Save data and summary data to disk
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

    ### Create and save sumtables
        setwd(OutputDirectory)
        dir.create("Output-sumtables")
        setwd("Output-sumtables")

        meta.dat
        as.matrix(names(cell.dat))

        Spectre::write.sumtables(x = cell.dat,
                                 sample.col = "Sample",
                                 pop.col = "FlowSOM_metacluster_42",
                                 measure.col = CellularCols,
                                 annot.col = names(cell.dat)[c(33:37)],
                                 group.col = "Group",
                                 do.frequencies = TRUE,
                                 cell.counts = meta.dat$Cells.per.sample,
                                 do.mfi.per.sample = TRUE,
                                 do.mfi.per.marker = TRUE)

#########################################################################################################
#### 7. Create and save some plots
##########################################################################################################

    ### Plot some sample-oriented plots
        setwd(OutputDirectory)
        dir.create("Output-plots")
        setwd("Output-plots")

        Spectre::make.factor.plot(dat = cell.dat.sub,
                                 x.axis = "UMAP_X",
                                 y.axis = "UMAP_Y",
                                 col.axis = "Sample",
                                 align.xy.by = cell.dat.sub,
                                 align.col.by = cell.dat.sub)

        Spectre::make.multi.plot(dat = cell.dat.sub,
                                 type = "factor",
                                 x.axis = "UMAP_X",
                                 y.axis = "UMAP_Y",
                                 col.axis = group.col,
                                 plot.by = sample.col,
                                 align.xy.by = cell.dat.sub,
                                 align.col.by = cell.dat.sub,
                                 dot.size = 1)

    ### Plot some cluster-oriented plots
        setwd(OutputDirectory)
        dir.create("Plots-clusters")
        setwd("Plots-clusters")

        Spectre::make.factor.plot(d = cell.dat.sub,
                                   x.axis = "UMAP_X",
                                   y.axis = "UMAP_Y",
                                   col.axis = "FlowSOM_metacluster_42",
                                   align.xy.by = cell.dat.sub,
                                   align.col.by = cell.dat.sub,
                                  add.label = TRUE)

    ### Plot some marker-oriented plots
        setwd(OutputDirectory)
        dir.create("Plots-markers")
        setwd("Plots-markers")

        dir.create("All samples")
        setwd("All samples")

        Spectre::make.multi.marker.plot(dat = cell.dat.sub,
                                   x.axis = "UMAP_X",
                                   y.axis = "UMAP_Y",
                                   plot.by = c(CellularCols),
                                   align.xy.by = cell.dat.sub,
                                   align.col.by = cell.dat.sub,
                                   figure.title = "Markers",
                                   dot.size = 1,
                                   save.each.plot = TRUE)

    ### Plot some marker-oriented plots -- one set per sample
        setwd(OutputDirectory)
        setwd("Plots-markers")
        dir.create("By sample")
        setwd("By sample")

        for(i in as.matrix(unique(cell.dat.sub[[sample.col]]))){
          for(a in ClusteringCols){
            Spectre::make.colour.plot(dat = cell.dat.sub[cell.dat.sub[[sample.col]] == i,],
                                 x.axis = "UMAP_X",
                                 y.axis = "UMAP_Y",
                                 col.axis = a,
                                 title = paste0(i, "_", a),
                                 align.xy.by = cell.dat.sub,
                                 align.col.by = cell.dat.sub)
          }
        }

    ### Plot some marker-oriented plots -- one set per group
        setwd(OutputDirectory)
        setwd("Plots-markers")
        dir.create("By group")
        setwd("By group")

        for(i in as.matrix(unique(cell.dat.sub[[group.col]]))){
          for(a in ClusteringCols){
            Spectre::make.colour.plot(dat = cell.dat.sub[cell.dat.sub[[group.col]] == i,],
                                 x.axis = "UMAP_X",
                                 y.axis = "UMAP_Y",
                                 col.axis = a,
                                 title = paste0(i, "_", a),
                                 align.xy.by = cell.dat.sub,
                                 align.col.by = cell.dat.sub)
          }
        }


