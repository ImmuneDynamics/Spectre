###################################################################################
### Spectre: spatial 3 - spatial analysis
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

    ### Set InputDirectory

        setwd(PrimaryDirectory)
        setwd("Output 2 - cellular analysis/")
        InputDirectory <- getwd()
        InputDirectory

    ### Create output directory

        setwd(PrimaryDirectory)
        dir.create("Output 3 - spatial analysis")
        setwd("Output 3 - spatial analysis")
        OutputDirectory <- getwd()

###################################################################################
### Read in data
###################################################################################

    ### Read in data.table
        
        setwd(InputDirectory)
        list.files(getwd(), '.csv')
        
        cell.dat <- fread("cell.dat.csv")
        cell.dat
        
    ### Read in spatial data
        
        setwd(PrimaryDirectory)
        setwd("Output 1 - add masks/QS file/")
        
        list.files(getwd(), '.qs')
    
        spatial.dat <- qread('spatial.dat.qs')
        
    ### Read in data.table

        area.table <- do.calculate.area(spatial.dat)
        area.table
        
        cell.dat$Region <- 'Total'
        cell.dat
        
###################################################################################
### Spatial analysis
###################################################################################

    setwd(OutputDirectory)

    ### Define columns
        
        as.matrix(names(cell.dat))
        
        roi.col <- 'ROI'
        sample.col <- 'Sample'
        group.col <- 'Group'
        
        pop.col <- 'Annotated metacluster'
        region.col <- 'Region'
        
    ### Some checks
        
        as.matrix(unique(cell.dat[[roi.col]]))
        as.matrix(unique(cell.dat[[group.col]]))
        as.matrix(unique(cell.dat[[region.col]]))
        area.table
        
        as.matrix(unique(cell.dat[[pop.col]]))
    
    ### Run area analysis on the DATA.TABLE

        reg.dat <- run.spatial.analysis(dat = cell.dat, 
                                        sample.col = roi.col, 
                                        pop.col = pop.col, 
                                        annot.cols = group.col, 
                                        region.col = region.col, 
                                        area.table = area.table) ## Also calculate on total by default
        
        reg.dat[,c(1:10)]
        
    ### Save regional analysis data
        
        setwd(OutputDirectory)
        fwrite(reg.dat, 'reg.dat.csv')
        
###################################################################################
### Plotting etc 
###################################################################################
        
    setwd(OutputDirectory)
        
    ### Plotting preferences
        
        as.matrix(names(reg.dat))
        
        as.matrix(names(reg.dat)[c(1:10)])
        to.plot <- names(reg.dat)[c(3:ncol(reg.dat))]
        to.plot

    ### Pheatmap

        reg.dat.z <- do.zscore(reg.dat, use.cols = to.plot, replace = TRUE)
        reg.dat.z
        
        reg.dat.z <- reg.dat.z[,colSums(is.na(reg.dat.z))<nrow(reg.dat.z), with = FALSE]
        reg.dat.z
         
        as.matrix(names(reg.dat.z))
        to.plot <- names(reg.dat.z)[c(3:ncol(reg.dat.z))]
        
        # names(reg.dat.z) <- gsub("Cells per region", "", names(reg.dat.z))
        # names(reg.dat.z) <- gsub("Cells per 100 um^2 of region", "", names(reg.dat.z))
        # names(reg.dat.z) <- gsub("Percent of cell type in sample", "", names(reg.dat.z))
        # names(reg.dat.z) <- gsub("Percent of cells in region", "", names(reg.dat.z))

        as.matrix(names(reg.dat.z))

        make.pheatmap(reg.dat.z, 
                      sample.col = 'ROI', 
                      plot.cols = to.plot, 
                      annot.cols = group.col,
                      is.fold = TRUE, 
                      dendrograms = 'both',
                      #row.sep = 2,
                      cutree_rows = 2, 
                      cutree_cols = 3)

    ### AutoGraphs

        setwd(OutputDirectory)
        dir.create("Autographs")
        setwd("Autographs")
        
        # meas.type <- unique(sub(" -- .*", "", names(sum.dat[,..plot.cols])))
        # meas.type
        
        unique(cell.dat$Group)
        as.matrix(names(reg.dat))
        
        for(i in names(reg.dat)[c(3:ncol(reg.dat))]){
            
            meas <- sub(".* -- ", "", i) # population
            pop <- sub(" -- .*", "", i) # measurement
            
            make.autograph(reg.dat,
                           x.axis = 'Group',
                           y.axis = i,
                           y.axis.label = meas,
                           
                           grp.order = c('Control', 'Test'),
                           my_comparisons = list(c('Control', 'Test')), 
                           
                           Pairwise_test = 't.test',
                           
                           title = pop,
                           #subtitle = meas
            )
        }
        
        
