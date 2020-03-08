##########################################################################################################
#### Spectre -- General Discovery Workflow
#### Part 2/3 - Annotation and re-summarise data
##########################################################################################################

    # Spectre R package: https://sydneycytometry.org.au/spectre
    # Thomas Myles Ashhurst, Felix Marsh-Wakefield, Givanna Putri

##########################################################################################################
#### SETUP
##########################################################################################################

    ### Load packages from library

        library(Spectre)
        Spectre::package.check() # --> change so that message at the end is "All required packages have been successfully installed"
        Spectre::package.load() # --> change so that message at the end is "All required packages have been successfully loaded"

        session_info()

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

    ### Determine input directory
        setwd("Output_Spectre/Output-data/")
        InputDirectory <- getwd()
        setwd(PrimaryDirectory)

    ### Create output directory
        setwd(PrimaryDirectory)
        dir.create("Output_Spectre", showWarnings = FALSE)
        setwd("Output_Spectre")

        dir.create("Output-annotated", showWarnings = FALSE)
        setwd("Output-annotated")

        OutputDirectory <- getwd()
        setwd(PrimaryDirectory)

##########################################################################################################
#### LOAD DATA
##########################################################################################################

    ### Read data into R
        setwd(InputDirectory)
        list.files(path = getwd(), pattern = ".csv")

        cell.dat <- fread("Clustered_TAdemo.csv")
        cell.dat

        cell.dat.sub <-fread("DimRed_TAdemo.csv")
        cell.dat.sub

    ### Read in annotation list

        setwd(PrimaryDirectory)
        setwd("metadata")

        list.files(path = getwd(), pattern = ".csv")
        annotations <- read.csv("annotation.csv")
        annotations

    ### Read in other stuff

        setwd(PrimaryDirectory)
        setwd("Output_Spectre/Output-info")
        list.files(path = getwd(), pattern = ".csv")

        meta.dat <- read.csv("metadata.csv")
        meta.dat

        CellularCols <- read.csv("CellularCols.csv")
        CellularCols <- as.vector(CellularCols[[1]])
        CellularCols

        ClusteringCols <- read.csv("ClusteringCols.csv")
        ClusteringCols <- as.vector(ClusteringCols[[1]])
        ClusteringCols

##########################################################################################################
#### ADD ANNOTATIONS
##########################################################################################################

    ### Embed population name columns in cell.dat
        cell.dat <- Spectre::do.embed.columns(x = cell.dat,
                                              type = "data.table",
                                              base.name = "FlowSOM_metacluster_42",
                                              col.name = "Population",
                                              match.to = annotations[,1],
                                              new.cols = annotations[,2])

        cell.dat

    ### Embed population name columns in cell.dat.sub

        cell.dat.sub <- Spectre::do.embed.columns(x = cell.dat.sub,
                                                 type = "data.table",
                                                 base.name = "FlowSOM_metacluster_42",
                                                 col.name = "Population",
                                                 match.to = annotations[,1],
                                                 new.cols = annotations[,2])

        cell.dat.sub

##########################################################################################################
#### GENERATE SOME NEW DATA AND SUMMARY DATA (based on population names)
##########################################################################################################

    ### Write 'large' dataset
        setwd(OutputDirectory)
        dir.create("Annotated-data", showWarnings = FALSE)
        setwd("Annotated-data")

        Spectre::write.files(dat = cell.dat,
                             file.prefix= paste0("Clustered"), # required
                             write.csv = TRUE,
                             write.fcs = TRUE)

        Spectre::write.files(dat = cell.dat,
                             file.prefix= paste0("Clustered"), # required
                             divide.by = "Sample",
                             write.csv = FALSE,
                             write.fcs = TRUE)

    ### Write 'subsample' dataset
        setwd(OutputDirectory)
        dir.create("Annotated-data", showWarnings = FALSE)
        setwd("Annotated-data")

        Spectre::write.files(dat = cell.dat.sub,
                             file.prefix = paste0("DimRed"), # required
                             write.csv = TRUE,
                             write.fcs = TRUE)

        Spectre::write.files(dat = cell.dat.sub,
                             file.prefix = paste0("DimRed"), # required
                             divide.by = "Sample",
                             write.csv = FALSE,
                             write.fcs = TRUE)

        setwd(PrimaryDirectory)


    ### Set positive cut offs for selected markers

        setwd(OutputDirectory)
        dir.create("Annotated-sumtables", showWarnings = FALSE)
        setwd("Annotated-sumtables")
        dir.create("SumTable-PercentPositive")
        setwd("SumTable-PercentPositive")

        as.matrix(names(cell.dat))

        make.density.plot(dat = cell.dat.sub,
                          x.axis = "BV711.SCA.1",
                          y.axis = "BV605.Ly6C",
                          title = paste0("Density-", "BV711.SCA.1"),
                          save.to.disk = TRUE)

        make.density.plot(dat = cell.dat.sub,
                          x.axis = "APC.BrdU",
                          y.axis = "BV605.Ly6C",
                          title = paste0("Density-", "APC.BrdU"),
                          save.to.disk = TRUE)

        # write.sumtable.percent.pos(x = cell.dat,
        #                            sample.name = "Sample",
        #                            cluster.name = "Population",
        #                            Markers = c("BV711.SCA.1","APC.BrdU"),
        #                            Cutoff = c(580, 450)
        #                            )

    ### Generate summary data
        setwd(OutputDirectory)
        dir.create("Annotated-sumtables", showWarnings = FALSE)
        setwd("Annotated-sumtables")

        as.matrix(names(cell.dat))
        cell.counts <- c(meta.dat[[5]])

        Spectre::write.sumtables(x = cell.dat,
                                 sample.col = "Sample",
                                 pop.col = "Population",
                                 measure.col = CellularCols,
                                 annot.col = names(cell.dat)[c(33:37)],
                                 group.col = "Group",
                                 do.frequencies = TRUE,
                                 cell.counts = cell.counts,
                                 do.mfi.per.sample = TRUE,
                                 do.mfi.per.marker = TRUE)

