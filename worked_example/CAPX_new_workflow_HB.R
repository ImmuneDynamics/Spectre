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

        #install_github("sydneycytometry/spectre")
        install_github("sydneycytometry/spectre", ref = 'development') # option to install the development verison if required
        library("Spectre")

    ### 1.2. Install other packages

        ## Install BiocManager to download packages from Bioconductor
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")

        ## Download additional BioConductor packages
        if(!require('flowCore')) {BiocManager::install('flowCore')}
        if(!require('Biobase')) {BiocManager::install('Biobase')}
        if(!require('flowViz')) {BiocManager::install('flowViz')}
        if(!require('FlowSOM')) {BiocManager::install('FlowSOM')}

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
        
        
        for(i in names(data.list)){
          data.list[[i]]$SampleID <- NULL
        }
        
        ## Save starting data
        data.start <- data.list

    ### Read sample metadata and embed in sample data
        meta.dat <- read.csv(file = "Metadata/SampleMetadata-BM.csv")
        meta.dat
        
        for(i in c(2:length(names(meta.dat)))){
          Spectre::embed.columns(x = data.list,
                                 type = "list",
                                 match.to = meta.dat[c(1)],
                                 new.cols = meta.dat[c(i)],
                                 col.name = names(meta.dat[c(i)]))
        }

        head(data.list)

    ### Merge files

        ## Merge files and review
        Spectre::file.merge(x = data.list)

        str(cell.dat)
        head(cell.dat)
        dim(cell.dat)
        
        unique(cell.dat$SampleName)
        
        as.matrix(unique(cell.dat[["SampleName"]]))

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
        exp.name <- "HnsB"

        file.col <- "FileName"
        sample.col <- "SampleName"
        group.col <- "Group"
        batch.col <- "Batch"

        ## Create a list of column names
        ColumnNames <- as.matrix(unname(colnames(cell.dat))) # assign reporter and marker names (column names) to 'ColumnNames'
        ColumnNames

    ### Define cellular and clustering columns

        ## Define columns that are 'valid' cellular markers (i.e. not live/dead, blank channels etc)
        ColumnNames # view the column 'number' for each parameter

        CellularColsNos <- c(5:14,16:17,19:22,25:26,28:35,37:38,46:49,51:58)
        CellularCols <- ColumnNames[CellularColsNos]

        CellularCols  # check that the column names that appear are the ones you want to analyse
        ColumnNames[-CellularColsNos] # Check which columns are being EXCLUDED!

        setwd(PrimaryDirectory)
        write.csv(x = CellularCols, file = "Metadata/BM_cellular_cols.csv")
        
    ### Define columns for clustering

        ## Define columns that are 'valid' cellular markers (i.e. not live/dead, blank channels etc)
        ColumnNames
        ClusteringColNos <- c(5,7,9:13,16:17,19:21,28:34,37:38,46:49,51:53,56:58)
        ClusteringCols <- ColumnNames[ClusteringColNos] # e.g. [c(11, 23, 10)] to include the markers corresponding to the column numbers 11, 23, 10

        ClusteringCols  # check that the column names that appear are the ones you want to analyse
        ColumnNames[-ClusteringColNos] # Check which columns are being EXCLUDED!

        setwd(PrimaryDirectory)
        write.csv(x = ClusteringCols, file = "Metadata/BM_clustering_cols.csv")
        
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
                             meta.k = 42,
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
                           targets = c(rep(2725,17)),
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
        Spectre::factor.plot(d = cell.dat.sub,
                                     x.axis = "UMAP_X",
                                     y.axis = "UMAP_Y",
                                     col.axis = group.col,
                                     title = paste0("Groups"),
                                     align.xy.by = cell.dat.sub,
                                     align.col.by = cell.dat.sub
                                     )

        
##########################################################################################################
#### Save data to disk
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
                             write.csv = FALSE,
                             write.fcs = TRUE)
        
        
        ## Write 'subsample' dataset
        Spectre::write.files(x = cell.dat.sub,
                             file.prefix = paste0("DimRed_", exp.name), # required
                             write.csv = TRUE,
                             write.fcs = TRUE)
        
        Spectre::write.files(x = cell.dat.sub,
                             file.prefix = paste0("DimRed_", exp.name), # required
                             divide.by = "Sample",
                             write.csv = FALSE,
                             write.fcs = TRUE)
        
        setwd(PrimaryDirectory)
  
##########################################################################################################
#### Plotting
##########################################################################################################
        
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
                   dot.size = 1.5
                   )
        
        
        ## Plot by 'Sample', coloured by 'Batch'
        multi.plot(d = cell.dat.sub,
                   type = "factor",
                   x.axis = "UMAP_X",
                   y.axis = "UMAP_Y",
                   col.axis = batch.col,
                   plot.by = sample.col,
                   align.xy.by = cell.dat.sub,
                   align.col.by = cell.dat.sub,
                   dot.size = 1.5
        )

        ## Plot the whole dataset, one plot per marker in one image
        multi.marker.plot(d = cell.dat.sub,
                          x.axis = "UMAP_X",
                          y.axis = "UMAP_Y",
                          plot.by = c(CellularCols),
                          align.xy.by = cell.dat.sub,
                          align.col.by = cell.dat.sub,
                          figure.title = "Markers",
                          dot.size = 1.5)

        ## Plot by 'Sample', coloured by 'marker' levels
        for(i in CellularCols){
            multi.plot(d = cell.dat.sub,
                       type = "colour",
                       x.axis = "UMAP_X",
                       y.axis = "UMAP_Y",
                       col.axis = i,
                       plot.by = sample.col,
                       figure.title = i,
                       align.xy.by = cell.dat.sub,
                       align.col.by = cell.dat.sub,
                       dot.size = 1.5
                       )
            }

        
        multi.plot(d = cell.dat.sub,
                   type = "factor",
                   x.axis = "UMAP_X",
                   y.axis = "UMAP_Y",
                   col.axis = "FlowSOM_metacluster",
                   plot.by = sample.col,
                   align.xy.by = cell.dat.sub,
                   align.col.by = cell.dat.sub,
                   dot.size = 1.5
        )

        p <- labelled.factor.plot(d = cell.dat.sub,
                             x.axis = "UMAP_X",
                             y.axis = "UMAP_Y",
                             col.axis = "FlowSOM_metacluster",
                             title = "FlowSOM_metacluster",
                             align.xy.by = cell.dat.sub,
                             align.col.by = cell.dat.sub,
                             dot.size = 1)
        
        ggsave(filename = paste0("Clusters", ".png"),
               plot = p,
               path = OutputDirectory,
               width = 9,
               height = 7,
               limitsize = FALSE)
        
        p <- colour.plot(d = cell.dat.sub,
                                  x.axis = "UMAP_X",
                                  y.axis = "UMAP_Y",
                                  col.axis = "176Yb_Ly6C",
                                  title = "176Yb_Ly6C",
                                  align.xy.by = cell.dat.sub,
                                  align.col.by = cell.dat.sub,
                                  dot.size = 1)
        
        ggsave(filename = paste0("176Yb_Ly6C", ".png"),
               plot = p,
               path = OutputDirectory,
               width = 9,
               height = 7,
               limitsize = FALSE)
        


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
                           sample.col = "SampleName",
                           clust.col = "FlowSOM_metacluster",
                           annot.col.nums = c(45,60:69),
                           
                           do.frequencies = TRUE,
                           cells.per.tissue = NULL, ## CHECK THE ORDER OF OCCURANCE OF THE unique entries in the sample.col COLUMN -- the order or cell counts in the vector MUST be the same
                           
                           do.exp.per.marker = TRUE,
                           do.exp.per.sample = TRUE,
                           fun.type = "median",
                           
                           do.foldchange = TRUE,
                           group.col = "Group",
                           control.group = "Air"
                           )
    

### Make an expression heatmap
    
    
    
    
    
### Make some cell freq heatmaps
    
    setwd(OutputDirectory)
    setwd("Output-SumTables")
    setwd("Output_CellNums")
    
    list.files(path = getwd(), pattern = ".csv")
    
    dat <- read.csv("SumTable_Proportions_FoldChangeLog2.csv")
    dat
    
    dat[1] <- NULL
    dat
    
    Batches <- c(1,1,1,1,2,2,
                 1,2,2,2,1,2,
                 2,2,2,2,2) # in order of sample (air, Smoke, Smoke + FT)
    
    dat <- cbind(dat, Batches)
    dat
    
    
    ### Frequency FOLD CHANGE heatmap loop
    
    as.matrix(names(dat))
    Spectre::make.pheatmap(x = dat,
                           file.name = "CellFreqFOLD.png",
                           plot.title = "Cell Frequency Fold Change (Log 2)",
                           sample.col = "SampleName",
                           annot.cols = c(2,45),
                           rmv.cols = c(1), 
                           is.fold = TRUE, 
                           dendrograms = "none", row.sep = c(6,12)) 
    
### MFI heatmap loop
    
    setwd(OutputDirectory)
    setwd("Output-SumTables/")
    setwd("Output_MFI_per_marker/")
    
    files <- list.files(path = getwd(), pattern = ".csv")
    files
    
    # Heatmap loop
    
    for(i in files){
      dat <- read.csv(i)
      dat
      
      Batches <- c(1,1,1,1,2,2,
                   1,2,2,2,1,2,
                   2,2,2,2,2) # in order of sample (air, Smoke, Smoke + FT)
      
      Groups <- c(rep("Air", 6),
                  rep("Smoke", 6),
                  rep("Smoke + FT", 5)
      )
      
      dat <- cbind(dat, Batches, Groups)
      dat
      
      
      i <- gsub(x = i, pattern = ".csv", "")
      
      as.matrix(names(dat))
      Spectre::make.pheatmap(x = dat,
                             file.name = paste0(i, ".png"),
                             plot.title = paste0(i, ""),
                             sample.col = c(1),
                             annot.cols = c(44,45),
                             rmv.cols = c(1), 
                             row.sep = c(6,12),
                             
                             is.fold = FALSE, 
                             dendrograms = "none") 
    }
    
        
