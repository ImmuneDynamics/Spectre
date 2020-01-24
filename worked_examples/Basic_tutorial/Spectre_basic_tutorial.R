##################################################################################################
##### Spectre basic tutorial #####################################################################
##################################################################################################

### 1. Install and load packages

    library(Spectre)

    Spectre::check.packages()
    Spectre::load.packages()

### 2. Set working directory

    setwd("/Users/thomasa/Desktop/")

### 3. Load data and set preferences

    ## Use demo dataset (included in Spectre package)
    cell.dat <- as.data.table(demo.start)
    cell.dat

    ## Load your own data
    #...

    ## Define columns for clustering
    #...

### 4. Run FlowSOM

    cell.dat <- Spectre::run.flowsom(x = cell.dat,
                                     clustering.cols = )#####

### 5. Run UMAP

    cell.dat <- Spectre::run.umap(x = cell.dat,
                                  use.cols = )####

    plot(cell.dat$UMAP_X, cell.dat$UMAP_Y)

### 6. Save the data to disk

    Spectre::write.files()

### 7. Save some plots

    Spectre::colour.plot()
    Spectre::multi.marker.plot()
    Spectre::multi.plot()


