##########################################################################################################
#### DRAFT Cytometry Analysis Pipeline for large and compleX data (CAPX) v3.0 - using the Spectre R package
#### Part 2 - Plotting, data exploration, cluster annotation
##########################################################################################################

# Thomas Myles Ashhurst, Felix Marsh-Wakefield
# 2019-08-02
# Workflow: https://sydneycytometry.org.au/capx
# Spectre R package: https://sydneycytometry.org.au/spectre

##########################################################################################################
#### 1. Install packages, load packages, and set working directory
#########################################################################################################

    ### 1.1. Load 'Spectre' package (using devtools)
    if(!require('devtools')) {install.packages('devtools')}
    library('devtools')

    #if(!require('sydneycytometry/spectre')) {install_github("sydneycytometry/spectre")}
    library("Spectre")

    # We recommend not to update packages that are dependencies of Spectre

    ### 1.2. Install packages
    if(!require("colorRamps")){install.packages("colorRamps")}
    if(!require("scales")){install.packages("scales")}
    if(!require("RColorBrewer")){install.packages("RColorBrewer")}
    if(!require("pheatmap")){install.packages("pheatmap")}

    ### 1.3. Load packages from library
    library('colorRamps')
    library('scales')
    library('RColorBrewer')
    library("pheatmap")


    ## 1.4. Set working directory

    ## Set working directory
    dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
    getwd()

    setwd("./Output_CAPX/Output-SumTables/")
    list.files(getwd(), ".csv")

    PrimaryDirectory <- getwd()
    PrimaryDirectory

    # ## Can set manually using these lines, if desired
    # PrimaryDirectory <- "/Users/Tom/Google Drive (t.ashhurst@centenary.org.au)/_Sydney Cytometry/2019_Synced/GitHub/Public github/Spectre/Workflow scripts/Output_CAPX/Output-DimRedData/"
    # setwd(PrimaryDirectory)

    ## Create output directory
    dir.create("Output_CAPX_heatmaps", showWarnings = FALSE)
    setwd("Output_CAPX_heatmaps")
    OutputDirectory <- getwd()
    setwd(PrimaryDirectory)

##########################################################################################################
#### 2. Heatmap creation
#########################################################################################################

    ### Heatmaps (cell proportions/counts)
    list.files(PrimaryDirectory, ".csv")

    ##
    cell.prop <- read.csv("SumTable_CellsPerTissue.csv")
    cell.prop.fold  <- read.csv("SumTable_CellsPerTissue_FoldChangeLog2.csv")

    ## Replace any '-Inf' with NA
    cell.prop[cell.prop == "-Inf"] <- NA
    cell.prop.fold[cell.prop.fold == "-Inf"] <- NA

    ### Plot 'normal' heatmap for cell.prop
    cell.prop
    names(cell.prop)

    setwd(OutputDirectory)
    Spectre::make.pheatmap(x = cell.prop,
                           file.name = "Cell numbers.png",
                           sample.col = "Sample",
                           annot.cols = 2,
                           rmv.cols = c(1,3),
                            plot.title = "Cell counts",
                            cell.size = 15
                            )

    ### Plot 'fold' heatmap for cell.prop.fold
    cell.prop.fold
    names(cell.prop.fold)

    setwd(OutputDirectory)
    Spectre::make.pheatmap(x = cell.prop.fold,
                           file.name = "Cell numbers - fold.png",
                           is.fold = TRUE,
                           fold.max.range = 3,
                           fold.min.range = -3,
                           sample.col = "Sample",
                           annot.cols = 3,
                           rmv.cols = c(1,2),
                           plot.title = "Cell counts (fold-change, log2)",
                           cell.size = 15
    )


##########################################################################################################
#### 3. Heatmap loops
#########################################################################################################

    ## Read in data
    #norm.list <- as.matrix(list.files(PrimaryDirectory, ".csv")[c(5)])
    #fold.list <- as.matrix(list.files(PrimaryDirectory, ".csv")[c(3)])
