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
        setwd("ROIs/")
        InputDirectory <- getwd()
        InputDirectory
        
    ### Create output directory

        setwd(PrimaryDirectory)
        dir.create("QS")
        setwd("QS")
        OutputDirectory <- getwd()
        OutputDirectory
        
    ### Create Mask directory
        
        setwd(PrimaryDirectory)
        dir.create("masks")
        setwd("masks")
        MaskDirectory <- getwd()
        MaskDirectory
        
    ### Cropped HDF5 for Ilastik
        
        setwd(PrimaryDirectory)
        dir.create("cropped")
        setwd("cropped")
        CroppedDirectory <- getwd()
        CroppedDirectory
        
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
        
        t(as.data.frame(tiff.list))
        
###################################################################################
### Read in TIFF files and create spatial objects
###################################################################################        

        setwd(InputDirectory)
        
        spatial.dat <- read.spatial.files(rois = rois, roi.loc = getwd())
        
        str(spatial.dat, 3)
        spatial.dat[[1]]$RASTERS

###################################################################################
### Create HDF5 files for mask creation
###################################################################################          
        
        nms <- names(spatial.dat[[1]]$RASTERS)
        as.matrix(nms)
        
        for.ilastik <- nms[c(3:15)]
        as.matrix(for.ilastik)
        
        ## Whole ROIs for Ilastik
        setwd(MaskDirectory)
        write.hdf5(spatial.dat, channels = for.ilastik)
        
        ## Cropped ROIs to train Ilastik
        setwd(CroppedDirectory)
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
        
        