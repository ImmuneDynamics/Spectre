##########################################################################################################
#### Spectre - statistical, quantitative, and differential analysis
##########################################################################################################

    # Spectre R package: https://github.com/immunedynamics/spectre
    # Thomas Myles Ashhurst, Felix Marsh-Wakefield, Givanna Putri

##########################################################################################################
#### 1. Load packages, and set working directory
##########################################################################################################

    ### Load libraries

        library(Spectre)
        Spectre::package.check()    # Check that all required packages are installed
        Spectre::package.load()     # Load required packages

    ### Set PrimaryDirectory
        
        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        getwd()
        PrimaryDirectory <- getwd()

    ### Create output directory
        
        dir.create("Output_Spectre_stats", showWarnings = FALSE)
        setwd("Output_Spectre_stats")
        OutputDirectory <- getwd()
        setwd(PrimaryDirectory)
        
##########################################################################################################
#### 2. Import summary data
##########################################################################################################

    ### Import data
        
        setwd(PrimaryDirectory)
        sum.dat <- fread("sum.dat.csv")
        
    ### Columns for analysis
        
        sum.dat
        as.matrix(names(sum.dat))
        
        sample.col <- names(sum.dat)[c(1)]
        group.col <- names(sum.dat)[c(2)]
        annot.cols <- names(sum.dat)[c(2:3)]
        
        plot.cols <- names(sum.dat)[c(4:15)]
        as.matrix(plot.cols)
        
    ### Experimental groups
        
        as.matrix(unique(sum.dat[[group.col]]))
        
        grp.order <- c("Mock", "WNV", "WNV + Rx")
        as.matrix(grp.order)
        
        comparisons <- list(c("Mock", "WNV"),
                            c("WNV", "WNV + Rx"),
                            c("Mock", "WNV + Rx")
                            )
        comparisons
        
    ### Setup
    
        variance.test <- 'kruskal.test'
        pairwise.test <- "wilcox.test"
        
    ### Reorder summary data and SAVE
        
        sum.dat <- do.reorder(sum.dat, group.col, grp.order)
        
        sum.dat
        sum.dat[,c(1:5)]
        
##########################################################################################################
#### 3. Output analysis
##########################################################################################################

    setwd(OutputDirectory)
    fwrite(sum.dat, 'sum.dat.csv')
    
    ### Autographs

        for(i in plot.cols){

            make.autograph(sum.dat,
                           x.axis = group.col,
                           y.axis = i,
                           y.axis.label = i,
                           
                           grp.order = grp.order,
                           my_comparisons = comparisons,
                           
                           Variance_test = variance.test,
                           Pairwise_test = pairwise.test,
                           
                           title = i,
                           filename = paste0(i, '.pdf'))
            
        }
        
    ### Create a fold change heatmap
        
        ## Z-score calculation
        
        sum.dat.z <- do.zscore(sum.dat, plot.cols, replace = TRUE)
        
        ## Group 
        
        t.first <- match(grp.order, sum.dat.z[[group.col]])
        t.first <- t.first -1
        t.first
        
        ## Make heatmap
        make.pheatmap(sum.dat.z, 
                      sample.col = sample.col, 
                      plot.cols = plot.cols, 
                      is.fold = TRUE, 
                      plot.title = 'Z-score',
                      annot.cols = annot.cols,
                      dendrograms = 'column',
                      row.sep = t.first,
                      cutree_cols = 3)

##########################################################################################################
#### 4. Output session info
##########################################################################################################

    ### Session info and metadata
        
        setwd(OutputDirectory)
        dir.create("Output - info", showWarnings = FALSE)
        setwd("Output - info")

        sink(file = "session_info.txt", append=TRUE, split=FALSE, type = c("output", "message"))
        session_info()
        sink()

