###################################################################################
### Spectre: TIFF to FCS
###################################################################################
        
    ### Load libraries
        
        library('Spectre')
        
        Spectre::package.check(type = 'spatial')
        Spectre::package.load(type = 'spatial')
        
    ### Set PrimaryDirectory
        
        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        getwd()
        PrimaryDirectory <- getwd()
        PrimaryDirectory
        
    ### Set InputDirectory (ROI TIFFs)
        
        setwd(PrimaryDirectory)
        setwd("data/ROIs/")
        InputDirectory <- getwd()
        InputDirectory
        
    ### Set MaskDirectory (ROI mask TIFFs)
        
        setwd(PrimaryDirectory)
        setwd("data/masks/")
        MaskDirectory <- getwd()
        MaskDirectory
        
    ### Create output directory
        
        setwd(PrimaryDirectory)
        dir.create("Output 1 - add masks")
        setwd("Output 1 - add masks")
        OutputDirectory <- getwd()
        OutputDirectory

###################################################################################
### Check ROIs and TIFFs
###################################################################################
    
    ### Initialise the spatial data object with channel TIFF files
        
        setwd(InputDirectory)
        
        rois <- list.dirs(full.names = FALSE, recursive = FALSE)
        as.matrix(rois)
        
        tiff.list <- list()
    
    ### Check channel names
        
        for(i in rois){
            setwd(InputDirectory)
            setwd(i)
            tiff.list[[i]] <- list.files(getwd())
        }
        
        t(as.data.frame(tiff.list))

###################################################################################
### Read in TIFF files and create spatial objects
###################################################################################        
    
    ### Read in ROI channel TIFFs
        
        setwd(InputDirectory)
        spatial.dat <- read.spatial.files(rois = rois, roi.loc = InputDirectory, ext = '.tiff')
        
    ### Check results
        
        str(spatial.dat, 3)
        spatial.dat[[1]]$RASTERS

###################################################################################
### Read in masks files
###################################################################################
    
    ### Define cell mask extension for different mask types
        
        setwd(MaskDirectory)
        
        all.masks <- list.files(pattern = '.tiff')
        as.matrix(all.masks)
    
    ### Import CELL masks
        
        as.matrix(all.masks)
        cell.mask.ext <- '_Cell_mask.tiff'
        
        cell.masks <- list.files(pattern = cell.mask.ext)
        as.matrix(cell.masks)
        
        spatial.dat <- do.add.masks(spatial.dat = spatial.dat,
                                    mask.loc = getwd(),
                                    masks = cell.masks,
                                    mask.label = 'cell.mask',
                                    mask.ext = cell.mask.ext)
    
    ### Review data
        
        str(spatial.dat, 3)
        
        spatial.dat[[1]]$MASKS

###################################################################################
### Generate polygons and outlines
###################################################################################
    
    ### Generate polygons and outlines
    
        spatial.dat <- do.create.outlines(spatial.dat, 'cell.mask')
    
    ### Checks
        
        str(spatial.dat, 3)
        
        spatial.dat[[1]]$MASKS

###################################################################################
### Mask QC plots
###################################################################################       
    
    ### Mask plots
        
        setwd(OutputDirectory)
        dir.create('Plots - cell masks')
        setwd('Plots - cell masks')
        
        as.matrix(names(spatial.dat[[1]]$RASTERS))
        
        base <- 'DNA1_Ir191'
        mask <- 'cell.mask'
        
        for(i in names(spatial.dat)){
            make.spatial.plot(spatial.dat = spatial.dat, 
                              image.roi = i, 
                              image.channel = base, 
                              mask.outlines = 'cell.mask')
        }

###################################################################################
### Extract 'per cell' data
###################################################################################       

    ### Extract cellular data for each cell mask
    
        spatial.dat <- do.extract(spatial.dat, 'cell.mask')
        str(spatial.dat, 3)

###################################################################################
### Save spatial data object as qs file
###################################################################################       
    
    ### Output
        
        setwd(OutputDirectory)
        dir.create('QS file')
        setwd('QS file')
        
        qsave(spatial.dat, "spatial.dat.qs")

###################################################################################
### Write FCS files
###################################################################################       
    
    ### Pull cellular data
    
        cell.dat <- do.pull.data(spatial.dat, 'CellData')
        cell.dat
        
        as.matrix(names(cell.dat))
        cellular.cols <- names(cell.dat)[c(6:20)]
        as.matrix(cellular.cols)
    
    ### Arcsinh transformation
    
        cell.dat <- do.asinh(cell.dat, cellular.cols, cofactor = 1)
    
    ### Invert y axis
        
        all.neg <- function(test) -1*abs(test)
        
        y_invert <- cell.dat[['y']]
        y_invert <- all.neg(y_invert)
        cell.dat[['y_invert']] <- y_invert
        
        cell.dat
    
    ### Write FCS files  
        
        setwd(OutputDirectory)
        dir.create('FCS files')
        setwd('FCS files')
        
        write.files(cell.dat, 'spatial_all', write.csv = FALSE, write.fcs = TRUE)
        write.files(cell.dat, 'spatial', divide.by = 'ROI', write.csv = FALSE, write.fcs = TRUE)

###################################################################################
### Further plots
###################################################################################       
    
    ### Plot settings
        
        as.matrix(names(spatial.dat))
        
        plot.rois <- c("ROI002", "ROI004")
        plot.rois
        
        as.matrix(names(spatial.dat[[1]]$RASTERS))
        exp.plot <- names(spatial.dat[[1]]$RASTERS)[1:15]
        exp.plot
        
        as.matrix(names(spatial.dat[[1]]$RASTERS))
        factor.plot <- names(spatial.dat[[1]]$RASTERS)[0]
        factor.plot
    
    ### Expression data point plots
        
        setwd(OutputDirectory)
        dir.create('Plots - expression')
        setwd('Plots - expression')
        
        for(i in plot.rois){
            for(a in exp.plot){
                make.spatial.plot(spatial.dat = spatial.dat, 
                                  image.roi = i, 
                                  image.channel = a, 
                                  mask.outlines = 'cell.mask')
            }
        }
        
    ### Factor data point plots
        
        setwd(OutputDirectory)
        dir.create('Plots - data point expression')
        setwd('Plots - data point expression')
        
        for(i in plot.rois){
            for(a in exp.plot){
                make.spatial.plot(spatial.dat = spatial.dat, 
                                  image.roi = i, 
                                  image.channel = a, 
                                  mask.outlines = 'cell.mask', 
                                  cell.dat = 'CellData', 
                                  cell.col = a)
            }
            
            for(a in factor.plot){
                make.spatial.plot(spatial.dat = spatial.dat, 
                                  image.roi = i, 
                                  image.channel = a, 
                                  mask.outlines = 'cell.mask', 
                                  cell.dat = 'CellData', 
                                  cell.col = a, 
                                  cell.col.type = 'factor')
            }
        }

