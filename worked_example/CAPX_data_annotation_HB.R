##########################################################################################################
#### ANNOTATION, SUMMARY DATA, HEATMAPS, AUTOGRAPHS
##########################################################################################################

##########################################################################################################
#### SETUP
##########################################################################################################

### Load packages from library
    library('Spectre')
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
    
    if(!require('pheatmap')) {install.packages('pheatmap')}
    if(!require('ggpubr')) {install.packages('ggpubr')}
    library(pheatmap)
    library(ggpubr)

### Set working directory
    dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
    getwd()
    PrimaryDirectory <- getwd()
    PrimaryDirectory
    
### Create output directory
    dir.create("Output_Annotated", showWarnings = FALSE)
    setwd("Output_Annotated")
    OutputDirectory <- getwd()
    setwd(PrimaryDirectory)

### Read data into R
    list.files(path = PrimaryDirectory, pattern = ".csv")
    cell.dat <- fread("Clustered_HnsB.csv")
    cell.dat
    
    cell.dat.sub <-fread("DimRed_HnsB.csv")
    cell.dat.sub

    
##########################################################################################################
#### ADD ANNOTATIONS
##########################################################################################################
    
### Cluster numbers and population names
    
    clusters <- c(1:42)
    popnames <- c("Unknown",
                  "CD8a+ T cells",
                  "B cells",
                  "B cells",
                  "B cells",
                  "B cells",
                  "7",
                  "8",
                  "NK cells",
                  "10",
                  "CD11b+ cDC",
                  "12",
                  "Ly6G(hi) neutrophils",
                  "Basophils",
                  "Stem and progenitors",
                  "Stem and progenitors",
                  "Stem and progenitors",
                  "PDC",
                  "19",
                  "20",
                  "21",
                  "PDC",
                  "CD8a+ T cells",
                  "Neutrophil progenitors",
                  "25",
                  "26",
                  "27",
                  "PDC",
                  "Ly6G(hi) neutrophils",
                  "Ly6C(hi) monocytes",
                  "Neutrophil progenitors",
                  "Ly6G(hi) neutrophils",
                  "Ly6C(hi) monocytes",
                  "Ly6C(hi) monocytes",
                  "Ly6G(hi) neutrophils",
                  "Ly6G(hi) neutrophils",
                  "37",
                  "Ly6G(hi) neutrophils",
                  "39",
                  "40",
                  "CD11b+ cDC",
                  "Ly6C(hi) monocytes")

### Embed population name columns
    
    ## cell.dat
    Spectre::embed.columns(x = cell.dat, 
                           type = "data.frame",
                           base.name = "FlowSOM_metacluster",
                           col.name = "PopName",
                           match.to = clusters,
                           new.cols = popnames)
    
            cell.dat <- cbind(cell.dat, embed.res)
            cell.dat
    
    
    ## cell.dat.sub
    Spectre::embed.columns(x = cell.dat.sub, 
                           type = "data.frame",
                           base.name = "FlowSOM_metacluster",
                           col.name = "PopName",
                           match.to = clusters,
                           new.cols = popnames)
    
            cell.dat.sub <- cbind(cell.dat.sub, embed.res)
            cell.dat.sub
    
    
##########################################################################################################
#### GENERATE SOME NEW SUMMARY DATA (based on population names)
##########################################################################################################
    
### Set positive cut offs for selected markers
    
    setwd(OutputDirectory)
    dir.create("SumTable-PercentPos")
    setwd("SumTable-PercentPos")
    
    as.matrix(names(cell.dat))
    
    sumtable.percent.pos(x = cell.dat,
                         sample.name = "SampleName",
                         cluster.name = "PopName",
                         Markers = c("127I_IdU","164Dy_SCA_1", "174Yb_MHCII"),
                         Cutoff = c(450, 200, 150)
                         )
    
### Generate summary data
    setwd(OutputDirectory)
    
    as.matrix(names(cell.dat))
    Spectre::make.sumtable(x = cell.dat,
                           sample.col = "SampleName",
                           clust.col = "PopName",
                           annot.col.nums = c(45,60:70),
                           
                           do.frequencies = TRUE, 
                           cells.per.tissue = NULL, ## CHECK THE ORDER OF OCCURANCE OF THE unique entries in the sample.col COLUMN -- the order or cell counts in the vector MUST be the same
                           
                           do.exp.per.marker = TRUE,
                           do.exp.per.sample = TRUE,
                           fun.type = "median",
                           
                           do.foldchange = TRUE,
                           group.col = "Group",
                           control.group = "Air")   
 

##########################################################################################################
#### NEW PLOTS IF REQUIRED
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("New_plots")
    setwd("New_plots")
    
### Make a cluster plot
    p <- labelled.factor.plot(d = cell.dat.sub,
                             x.axis = "UMAP_X",
                             y.axis = "UMAP_Y",
                             col.axis = "PopName",
                             title = "Populations",
                             align.xy.by = cell.dat.sub,
                             align.col.by = cell.dat.sub,
                             dot.size = 1)
    
    ggsave(filename = paste0("Populations", ".png"),
           plot = p,
           path = getwd(),
           width = 9,
           height = 7,
           limitsize = FALSE)  
    
    
    multi.plot(d = cell.dat.sub,
               type = "labelled.factor",
               x.axis = "UMAP_X",
               y.axis = "UMAP_Y",
               col.axis = "PopName",
               plot.by = "Group",
               align.xy.by = cell.dat.sub,
               align.col.by = cell.dat.sub,
               dot.size = 1.5)
    
    
    

##########################################################################################################
#### FREQUENCY DATA
##########################################################################################################

    ### Read in and setup up FREQUENCY data
    
    setwd(OutputDirectory)
    setwd("Output_CellNums")
    
    files <- list.files(path = getwd(), pattern = ".csv")
    files

    cell.freq <- read.csv(files[1], check.names = FALSE)
    cell.freq[[1]] <- NULL
    
    cell.freq.fold <- read.csv(files[3], check.names = FALSE)
    cell.freq.fold[[1]] <- NULL
    
    Batches <- c(1,1,1,1,2,2,
                 1,2,2,2,1,2,
                 2,2,2,2,2) # in order of sample (air, Smoke, Smoke + FT)
    
    Groups <- c(rep("Air", 6),
                rep("Smoke", 6),
                rep("Smoke + FT", 5))
    
    cell.freq <- cbind(cell.freq, Batches, Groups)
    cell.freq.fold <- cbind(cell.freq.fold, Batches, Groups)
    
    as.matrix(names(cell.freq.fold))
    to.plot <- c(3:26)
    
    ## Generic preferneces
    sample.col <- "Sample"
    group.col <- "Groups"
    ctrl.grp <- "Air"
    
    
### FREQUENCY HEAMTPAS 
    setwd(OutputDirectory)
    setwd("Output_CellNums")
    
      dat <- cell.freq.fold
      as.matrix(names(dat))
    
      i <- gsub(x = i, pattern = ".csv", "")
      
      as.matrix(names(dat))
      Spectre::make.pheatmap(x = dat,
                             file.name = paste0("Frequency fold change.png"),
                             plot.title = paste0("Frequency fold change (log 2)"),
                             sample.col = "SampleName",
                             annot.cols = c(27,28),
                             rmv.cols = c(1,2), 
                             row.sep = c(6,12),
                             is.fold = TRUE, 
                             dendrograms = "none") 
    
    
    
### FREQUENCY AUTOGRAPHS
    setwd(OutputDirectory)
    setwd("Output_CellNums")
    
    dat <- cell.freq
    as.matrix(names(dat))
    
    dat[to.plot] <- dat[to.plot]*100
    
    for(a in to.plot){
      # a <- 3
      a <- names(dat)[[a]]
      i <- gsub(x = i, pattern = ".csv", "")
      
      make.autograph(x = dat, 
                     x.axis = group.col, 
                     y.axis = a, 
                     colour.by = "Batches",
                     colours = c("Black", "Red"), 
                     y.axis.label = "Percent positive", 
                     my_comparisons = my.comparisons,
                     title = paste0("Frequency of ", a), 
                     filename = paste0("Frequency of ", a, ".pdf"))
    }
    
    
    
    
    
    
##########################################################################################################
#### PERCENT POSITIVE DATA
##########################################################################################################

### Read in and setup up FREQUENCY data
    
    setwd(OutputDirectory)
    setwd("SumTable-PercentPos")
    
    files <- list.files(path = getwd(), pattern = ".csv")
    files
    
    file.list <- list()
    
    for(i in files){
      file.list[[i]] <- read.csv(i, check.names = FALSE)
    }
    
    file.list[[1]]
    as.matrix(names(file.list[[1]]))
    
    Batches <- c(1,1,1,1,2,2,
                 1,2,2,2,1,2,
                 2,2,2,2,2) # in order of sample (air, Smoke, Smoke + FT)
    
    Groups <- c(rep("Air", 6),
                rep("Smoke", 6),
                rep("Smoke + FT", 5))
    
    for(i in names(file.list)){
      file.list[[i]][1] <- NULL
      file.list[[i]] <- cbind(file.list[[i]], Batches, Groups)
    }
    
    file.list[[1]]
    
    as.matrix(names(file.list[[1]]))
    to.plot <- c(2:25)
    
    ## Generic preferneces
    sample.col <- "Sample"
    group.col <- "Groups"
    ctrl.grp <- "Air"
    
    ## AutoGraph preferences
    file.list[[1]]
    as.matrix(names(file.list[[1]]))
    unique(file.list[[1]][[group.col]])
    
    to.plot <- c(2:25)
    my.comparisons <- list(c("Air", "Smoke"),
                           c("Smoke", "Smoke + FT"),
                           c("Air", "Smoke + FT"))
        
    
### Percent positive FOLD CHANGE heatmap  
    
    setwd(OutputDirectory)
    setwd("SumTable-PercentPos")
    
    for(i in names(file.list)){
      # i <- names(file.list)[[1]]
      
      dat <- file.list[[i]]
      as.matrix(names(dat))
      
      dat <- convert.to.fold(x = dat, 
                              sample.col = sample.col, 
                              group.col = group.col, 
                              ctrl.grp = ctrl.grp, 
                              convert.cols = to.plot, 
                              log2 = TRUE)

      i <- gsub(x = i, pattern = ".csv", "")
      
      as.matrix(names(dat))
      Spectre::make.pheatmap(x = dat,
                             file.name = paste0(i, ".png"),
                             plot.title = paste0(i, "Percent positive"),
                             sample.col = "Sample",
                             annot.cols = c(26,27),
                             rmv.cols = c(1), 
                             row.sep = c(6,12),
                             is.fold = TRUE, 
                             dendrograms = "none") 
    }
    
    
    
### Percent positive AUTOGRAPHS
    
    setwd(OutputDirectory)
    setwd("SumTable-PercentPos")

    for(i in names(file.list)){
      # i <- "127I_IdU.csv"
      
      dat <- file.list[[i]]
      dat
      
      for(a in to.plot){
        # a <- 2
        a <- names(dat)[[a]]
        i <- gsub(x = i, pattern = ".csv", "")
        
        make.autograph(x = dat, 
                       x.axis = group.col, 
                       y.axis = a, 
                       colour.by = "Batches",
                       colours = c("Black", "Red"), 
                       y.axis.label = "Percent positive", 
                       my_comparisons = my.comparisons,
                       title = paste0(i, " - % pos of ", a), 
                       filename = paste0(i, " - PC pos of ", a, ".pdf")
                       )
      }
    }
    

    
    
##########################################################################################################
#### MFI HEATMAPS - per marker, samples vs clusters
##########################################################################################################

    setwd(OutputDirectory)
    setwd("Output_MFI_per_marker/")
    
    files <- list.files(path = getwd(), pattern = ".csv")
    files
    
    file.list <- list()
    
    for(i in files){
      file.list[[i]] <- read.csv(i, check.names = FALSE)
      colnames(file.list[[i]])[1] <- "Sample"
    }
    
    file.list[[1]]
    as.matrix(names(file.list[[1]]))
    
    Batches <- c(1,1,1,1,2,2,
                 1,2,2,2,1,2,
                 2,2,2,2,2) # in order of sample (air, Smoke, Smoke + FT)
    
    Groups <- c(rep("Air", 6),
                rep("Smoke", 6),
                rep("Smoke + FT", 5))
    
    for(i in names(file.list)){
      file.list[[i]] <- cbind(file.list[[i]], Batches, Groups)
    }
    
    file.list[[1]]
    
    as.matrix(names(file.list[[1]]))
    to.plot <- c(2:25)
    
    ## Generic preferneces
    sample.col <- "Sample"
    group.col <- "Groups"
    ctrl.grp <- "Air"
    
    ## AutoGraph preferences
    file.list[[1]]
    as.matrix(names(file.list[[1]]))
    unique(file.list[[1]][[group.col]])
    
    to.plot <- c(2:25)

### MFI heatmap loop
    
    setwd(OutputDirectory)
    setwd("Output_MFI_per_marker/")
 
    # Heatmap loop

    for(i in names(file.list)){
      
      dat <- file.list[[i]]
      
      dat <- convert.to.fold(x = dat, 
                             sample.col = sample.col, 
                             group.col = group.col, 
                             ctrl.grp = ctrl.grp, 
                             convert.cols = to.plot, 
                             log2 = TRUE)
      
      
      i <- gsub(x = i, pattern = ".csv", "")
      
      as.matrix(names(dat))
      Spectre::make.pheatmap(x = dat,
                             file.name = paste0(i, ".png"),
                             plot.title = paste0(i, ""),
                             sample.col = c(1),
                             annot.cols = c(26,27),
                             rmv.cols = c(1), 
                             row.sep = c(6,12),
                             
                             is.fold = TRUE, 
                             dendrograms = "none") 
    }
    

    
    
    
    
    
   
    
    
    