###################################################################################
###  Spectre: spatial 1 - Setup spatial data object
###################################################################################

    ### Load packages

        library(Spectre)

        package.check(type = 'spatial')
        package.load(type = 'spatial')

    ### Set PrimaryDirectory

        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        getwd()
        PrimaryDirectory <- getwd()
        PrimaryDirectory

    ### Set InputDirectory

        setwd(PrimaryDirectory)
        setwd("data/")
        InputDirectory <- getwd()
        InputDirectory
        
    ### Create output directory

        setwd(PrimaryDirectory)
        dir.create("Output 1 - Setup")
        setwd("Output 1 - Setup")
        OutputDirectory <- getwd()
        OutputDirectory
        
    ### Create Mask directory
        
        setwd(PrimaryDirectory)
        dir.create("masks")
        setwd("masks")
        MaskDirectory <- getwd()
        MaskDirectory
        
###################################################################################
### Check ROIs and TIFFs
###################################################################################

    ### Initialise the spatial data object with channel TIFF files

        setwd(InputDirectory)

        rois <- list.dirs(full.names = FALSE, recursive = FALSE)
        as.matrix(rois)

        tiff.list <- list()
        
        for(i in rois){
            
            setwd(InputDirectory)
            
            setwd(i)
            
            tiff.list[[i]] <- list.files(getwd())
            
        }
        
        tiff.list
        
###################################################################################
### Read in TIFF files and create spatial objects
###################################################################################        

        setwd(InputDirectory)
        
        spatial.dat <- read.spatial.files(rois = rois, roi.loc = getwd())
        str(spatial.dat, 3)

###################################################################################
### Create HDF5 file for Ilastik
###################################################################################          
        
        setwd(MaskDirectory)
        
        nms <- names(spatial.dat[[1]]$RASTERS)
        as.matrix(nms)
        
        for.ilastik <- nms[c(13:25)]
        as.matrix(for.ilastik)
        
        write.hdf5(spatial.dat, channels = for.ilastik, random.crop.x = 300, random.crop.y = 300)
        
###################################################################################
### Save spatial.dat object as RDS file
###################################################################################

    ### Save as quick serial (qs) file

        setwd(OutputDirectory)

        qsave(spatial.dat, "spatial.dat.qs")
        
###################################################################################
### Some quick QC spatial plots
###################################################################################

        # ### Choose an ROI to use
        # as.matrix(names(spatial.dat[[1]]$RASTERS))
        # 
        # roi.plot <- names(spatial.dat)[c(1)]
        # roi.plot
        
        