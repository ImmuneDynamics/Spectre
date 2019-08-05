##### Performing UMAP on cytometry data using the Spectre package #####
    # Thomas Myles Ashhurst
    # 2019-08-02
    # https://sydneycytometry.org.au/spectre

### 1. Install packages, load packages, and set workin directory

    ## 1.1. Install packages
        if(!require('plyr')) {install.packages('plyr')}
        if(!require('data.table')) {install.packages('data.table')}
        if(!require('rstudioapi')) {install.packages('rstudioapi')}
        if(!require('devtools')) {install.packages('devtools')}

        if(!require('flowViz')) {source("https://bioconductor.org/biocLite.R")
          biocLite('flowViz')}
        if(!require('flowCore')) {source("https://bioconductor.org/biocLite.R")
          biocLite('flowCore')}
        if(!require('Biobase')) {source("https://bioconductor.org/biocLite.R")
          biocLite('Biobase')}

        if(!require('FlowSOM')) {source("https://bioconductor.org/biocLite.R") # for running FlowSOM ### POSSIBLY INSTALL FROM GITHUB
          biocLite('FlowSOM')}
        if(!require('Rtsne')) {install.packages("Rtsne")} # for running tSNE
        if(!require('umap')) {install.packages('umap')}

        if (!require("ggplot2")){install.packages("ggplot2")} # for plotting tSNE graphs
        if (!require("colorRamps")){install.packages("colorRamps")} # for colour scheme management
        if (!require("ggthemes")){install.packages("ggthemes")} # for plot themes
        if (!require("scales")){install.packages("scales")} # for re-scaling if necessary
        if (!require("RColorBrewer")){install.packages("RColorBrewer")} # for re-scaling if necessary

    ## 1.2. Load packages
        library('plyr')
        library('data.table')
        library('rstudioapi')
        library('devtools')

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

    ## 1.3. Load 'Spectre' package

        #if(!require('sydneycytometry/spectre')) {install_github("sydneycytometry/spectre")}
        library("Spectre")

    ## 1.4. Set working directory
        PrimaryDirectory <- "/Users/Tom/Downloads/CAPX-2.5/Demo dataset/"
        setwd(PrimaryDirectory)

        list.files(PrimaryDirectory, ".csv")

### 2. Read and prepare data

    ## Read in samples
        Spectre::read.files(file.loc = PrimaryDirectory, file.type = ".csv", do.embed.file.names = TRUE)

    ## Review sample properties
        ncol.check
        nrow.check

        name.table
        head(data.list)
        head(data.list[[1]])

    ## Add sample names/numbers

        #CytoTools::file.annotate()

    ## Add groups
        as.matrix(names(data.list))

            group.names = list()
            group.nums = list()

            group.names[[1]] <- "Mock"
            group.names[[2]] <- "WNV"

            group.nums[[1]] <- c(1:6)
            group.nums[[2]] <- c(7:12)

        Spectre::add.groups(x = data.list, grp.names = group.names, grp.nums = group.nums)

        head(data.list)
        head(data.list[[7]])

    ## Merge samples
        Spectre::merge.files(x = data.list)
        head(cell.dat)

        # Look for any NAs

        # Cleanup (not necessary, but recommended)
            rm(data.list)
            rm(name.table)
            rm(ncol.check)
            rm(nrow.check)
            #rm(group.names)
            #rm(group.nums)

    dim(cell.dat)

    cell.dat.large <- rbind(cell.dat, cell.dat)
    cell.dat.large <- rbind(cell.dat.large, cell.dat) # x8 or so
    cell.dat <- cell.dat.large


    ## Define clustering columns

        ## Save column names
        ColumnNames <- unname(colnames(cell.dat)) # assign reporter and marker names (column names) to 'ColumnNames'
        as.matrix(ColumnNames) # view the column 'number' for each parameter

        ## Define columns for clustering
        ClusteringColNos <- c(5,6,8,9,11:13,17:19,21:29,32)
        ClusteringColNos

        ClusteringCols <- ColumnNames[ClusteringColNos] # e.g. [c(11, 23, 10)] to include the markers corresponding to the column numbers 11, 23, 10

        ##
        ClusteringCols  # check that the column names that appear are the ones you want to analyse
        ColumnNames[-ClusteringColNos] # Check which columns are being EXCLUDED!

        head(cell.dat)

### 3. Perform FlowSOM clustering

    ## FlowSOM
    Spectre::run.flowsom(x = cell.dat,
                         meta.k = 40,
                         clustering.cols = ClusteringCols,
                         clust.seed = 42,
                         meta.seed = 42,
                         clust.name = "FlowSOM_cluster",
                         meta.clust.name = "FlowSOM_metacluster")

    ## Add results
    cell.dat <- cbind(cell.dat, flowsom.res.original)
    cell.dat <- cbind(cell.dat, flowsom.res.meta)

    ## Check cell.dat
    head(cell.dat)

    ## Remove results
    rm(flowsom.res.original)
    rm(flowsom.res.meta)



    ########################

    # save cell.dat
    # create sumtables
    # save sumtables



### 4. Downsample data in preparation for dimensionality reduction

    ## Run subsample
    #Spectre::subsample(x = cell.dat,
    #                   method = "per.sample",
    #                   samp.col = "FileName",
    #                   targets = c(rep(500, 12)),
    #                   seed = 42)
    #
    #cell.dat.sub <- subsample.res
    #rm(subsample.res)

    cell.dat.sub <- cell.dat.large


### 5. Perform UMAP

    ## Check column names
    as.matrix(names(cell.dat.sub))

    ## Run UMAP
    Spectre::run.umap(x = cell.dat.sub,
                  use.cols = ClusteringCols,  #c(5,6,8,9,11,12,13,17:19,21:30,32),
                  umap.seed = 42)

    ## Merge UMAP results with data
    cell.dat.sub <- cbind(cell.dat.sub, umap.res)


### 6. Plot UMAP

    ## Plot UMAP results
    plot(cell.dat.sub$UMAP_42_X, cell.dat.sub$UMAP_42_Y)

    ## Plot single UMAP MFI plot
    p <- Spectre::colour.plot(d = cell.dat.sub,
                         x.axis = "UMAP_42_X",
                         y.axis = "UMAP_42_Y",
                         col.axis = "BV605.Ly6C",
                         title = "MFI",
                         colours = "spectral",
                         dot.size = 1)
    p

    head(cell.dat)

    ## Plot single UMAP factor plot
    p <- Spectre::factor.plot(d = cell.dat.sub,
                              x.axis = "UMAP_42_X",
                              y.axis = "UMAP_42_Y",
                              col.axis = "FileName",
                              title = "Samples",
                              dot.size = 1,
                              add.labels = FALSE) # assumes numeric

    p




    ## add option for centroid detection and labelling

    ## copy the # codes for 'spectral' colours -- use that directly, avoid needing Rcolourbrewer etc






##########################################################################################





    test <- function(d){
      ggplot(data = d, aes(x = d[[x.axis]], y = d[[y.axis]], colour = d[[col.axis]])) +
        geom_point(size = point.size) +
        scale_colour_gradientn(colours = colour.scheme(50),
                               limits = c(quantile(d[[col.axis]], probs = c(min.threshold)), #0.03-01 seems to work well
                                          quantile(d[[col.axis]], probs = c(max.threshold))), #0.97-995 seems to work well
                               oob=squish) +
        labs(colour = col.axis)+
        xlab(x.axis)+
        ylab(y.axis)+
        ggtitle(plot.title) +
        theme(
          panel.background = element_rect(fill = "white", colour = "black", size = 0.5) # change 'colour' to black for informative axis
          #axis.line=element_blank(),
          #axis.text.x=element_blank(),
          #axis.text.y=element_blank(),
          #axis.ticks=element_blank(),
          #axis.title.x=element_blank(),
          #axis.title.y=element_blank(),
          #panel.grid.major = element_blank(),
          #panel.background=element_blank(),
          #panel.border=element_blank(),
          #panel.grid.minor=element_blank(),
          #plot.background=element_blank(),
          #legend.position = "right",
          #legend.text=element_text(size=15), # large = 30 # small = 8
          #legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
          #legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
          #legend.title=element_blank()
          #plot.title = element_text(color="Black", face="bold", size=15, hjust=0) # size 70 for large, # 18 for small
        )
      }


    test(cell.dat.sub)




    #######

    jet <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    #spectral.list <- list(color = colorRampPalette(brewer.pal(11,"Spectral"))(100))
    spectral.list <- colorRampPalette(brewer.pal(11,"Spectral"))(50)
    spectral.list <- rev(spectral.list)

    spectral <- colorRampPalette(c(spectral.list))

    YlOrRd.list <- rev(colorRampPalette(brewer.pal(9,"YlOrRd"))(100))
    YltoRd.cols <- colorRampPalette(c(YlOrRd.list))

    BuPu.list <- rev(colorRampPalette(brewer.pal(9,"BuPu"))(100))
    BuPu.cols <- colorRampPalette(c(BuPu.list))


    # Custom, from spectral, Red to Blue -->
    custom <- colorRampPalette(rev(c("#7F0000", "#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B",
                                     "#FFFFBF",
                                     "#E6F598",
                                     "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2", "#00007F",
                                     "#170a45")))

    colour.palette.list <- (colorRampPalette(brewer.pal(9, "YlGnBu"))(31)) # 256
    colour.palette <- colorRampPalette(c(colour.palette.list))

    fold.palette <- colorRampPalette(rev(c("#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026","black","#023858","#045a8d","#0570b0","#3690c0","#74a9cf","#a6bddb","#d0d1e6","#ece7f2")))
    fold.palette <- colorRampPalette(rev(c("#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026","black")))
    #fold.palette <- colorRampPalette(c(fold.palette.list))


    colour.scheme <- jet
    #colour.scheme <- spectral
    #colour.scheme <- YltoRd.cols
    #colour.scheme <- BuPu.cols
    #colour.scheme <- custom
    #colour.scheme <- colour.palette
    #colour.scheme <- fold.palette


    head(cell.dat.sub)
    filenames <- unique(cell.dat.sub$FileName)


    ## Plot, coloured by samples
    ggplot(data = cell.dat.sub, aes(x = UMAP_42_X, y = UMAP_42_Y, colour = as.factor(cell.dat.sub$GroupName))) +
      geom_point(size = 1)+ # 2 for large # 0.5 for small
      #scale_colour_gradientn(colours = colour.scheme(50)) +
      #scale_colour_manual(name = "FileName", values = c(colour.scheme(length(filenames)))) +
      ggtitle("Samples") +
      theme(
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5)#, # change 'colour' to black for informative axis
        #axis.line=element_blank(),
        #axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks=element_blank(),
        #axis.title.x=element_blank(),
        #axis.title.y=element_blank(),
        #panel.grid.major = element_blank(),
        #panel.background=element_blank(),
        #panel.border=element_blank(),
        #panel.grid.minor=element_blank(),
        #plot.background=element_blank(),
        #legend.position = "right",
        #legend.text=element_text(size=15), # large = 30 # small = 8
        #legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
        #legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
        #legend.title=element_blank(),
        #plot.title = element_text(color="Black", face="bold", size=15, hjust=0) # size 70 for large, # 18 for small
      )




    head(cell.dat.sub)

    ## Plot, colour by MFI

          d <- cell.dat.sub

          x.axis <- "UMAP_42_X"
          y.axis <- "UMAP_42_Y"
          col.axis <- "BV711.SCA.1"

          min.threshold <- 0.01 # 0.01
          max.threshold <- 0.995 # 0.995 default

          plot.title <- "MFI"
          point.size <- 0.5 # 2 for large # 0.5 for small


    ggplot(data = d, aes(x = d[[x.axis]], y = d[[y.axis]], colour = d[[col.axis]])) +
      geom_point(size = point.size) +
      scale_colour_gradientn(colours = colour.scheme(50),
                             limits = c(quantile(d[[col.axis]], probs = c(min.threshold)), #0.03-01 seems to work well
                                        quantile(d[[col.axis]], probs = c(max.threshold))), #0.97-995 seems to work well
                             oob=squish) +
      labs(colour = col.axis)+
      xlab(x.axis)+
      ylab(y.axis)+
      ggtitle(plot.title) +
      theme(
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5) # change 'colour' to black for informative axis
        #axis.line=element_blank(),
        #axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks=element_blank(),
        #axis.title.x=element_blank(),
        #axis.title.y=element_blank(),
        #panel.grid.major = element_blank(),
        #panel.background=element_blank(),
        #panel.border=element_blank(),
        #panel.grid.minor=element_blank(),
        #plot.background=element_blank(),
        #legend.position = "right",
        #legend.text=element_text(size=15), # large = 30 # small = 8
        #legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
        #legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
        #legend.title=element_blank()
        #plot.title = element_text(color="Black", face="bold", size=15, hjust=0) # size 70 for large, # 18 for small
      )




