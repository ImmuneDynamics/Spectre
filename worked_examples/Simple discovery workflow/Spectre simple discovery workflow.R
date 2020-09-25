##########################################################################################################
#### Spectre - Simple Discovery Workflow
#### Clustering, dimensionality reduction, plotting, and summarise data
##########################################################################################################

    # Spectre R package: https://sydneycytometry.org.au/spectre
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

    ### Set 'input' directory
        setwd(PrimaryDirectory)
        setwd("data/")
        InputDirectory <- getwd()
        setwd(PrimaryDirectory)

    ### Set 'metadata' directory
        setwd(PrimaryDirectory)
        setwd("metadata/")
        MetaDirectory <- getwd()
        setwd(PrimaryDirectory)

    ### Create output directory
        dir.create("Output_Spectre", showWarnings = FALSE)
        setwd("Output_Spectre")
        OutputDirectory <- getwd()
        setwd(PrimaryDirectory)

##########################################################################################################
#### 2. Import and prep data
##########################################################################################################

    ### Import data

        setwd(InputDirectory)
        list.files(InputDirectory, ".csv")

        data.list <- Spectre::read.files(file.loc = InputDirectory,
                                         file.type = ".csv",
                                         do.embed.file.names = TRUE)

    ### Check the data

        check <- do.list.summary(data.list)

        check$name.table # Review column names and their subsequent values
        check$ncol.check # Review number of columns (features, markers) in each sample
        check$nrow.check # Review number of rows (cells) in each sample

        data.list[[1]]

    ### Merge data

        cell.dat <- Spectre::do.merge.files(dat = data.list)
        cell.dat

    ### Arcsinh transformation

        as.matrix(names(cell.dat))

        to.asinh <- names(cell.dat)[c(1:9)]
        to.asinh

        cofactor <- 500
        plot.against <- "Ly6C_asinh"

        cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)
        transformed.cols <- paste0(to.asinh, "_asinh")

        setwd(OutputDirectory)
        dir.create("Output - transformed plots")
        setwd("Output - transformed plots")

        for(i in transformed.cols){
          make.colour.plot(do.subsample(cell.dat, 20000), i, plot.against)
        }

##########################################################################################################
#### 3. Add metadata and set some preferences
##########################################################################################################

    setwd(MetaDirectory)

    ### Metadata

        meta.dat <- fread("sample.details.csv")
        meta.dat

        meta.dat <- meta.dat[,c(1:4)]
        meta.dat

        cell.dat <- do.add.cols(cell.dat, "FileName", meta.dat, "Filename", rmv.ext = TRUE)
        cell.dat

    ### Columns

        as.matrix(names(cell.dat))

        cellular.cols <- names(cell.dat)[c(12:20)]
        as.matrix(cellular.cols)

        cluster.cols <- names(cell.dat)[c(12:20)]
        as.matrix(cluster.cols)

        exp.name <- "CNS experiment"
        sample.col <- "Sample"
        group.col <- "Group"
        batch.col <- "Batch"

    ### Subsample targets per group

        data.frame(table(cell.dat[[group.col]])) # Check number of cells per sample.

        sub.targets <- c(2000, 20000) # target subsample numbers from each group
        sub.targets

##########################################################################################################
#### 4. Clustering and dimensionality reduction
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output - clustering")
    setwd("Output - clustering")

    ### Clustering

        cell.dat <- run.flowsom(cell.dat, cluster.cols)
        fwrite(cell.dat, "clustered.data.csv")

    ### Dimensionality reduction

        cell.sub <- do.subsample(cell.dat, sub.targets, group.col)
        cell.sub <- run.umap(cell.sub, cluster.cols)

        fwrite(cell.dat, "clustered.data.DR.csv")

    ### DR plots

        make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", col.type = 'factor', add.label = TRUE)
        make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", cellular.cols)
        make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", group.col, col.type = 'factor')

    ### Expression heatmap

        exp <- do.aggregate(cell.dat, cellular.cols, by = "FlowSOM_metacluster")
        make.pheatmap(exp, "FlowSOM_metacluster", cellular.cols)

##########################################################################################################
#### 5. Annotate clusters
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output - annotation")
    setwd("Output - annotation")

    ### Annotate

        annots <- list("Microglia" = c(2,4),
                       "Macrophages" = c(3),
                       "Neutrophils" = c(1),
                       "NK cells" = c(6),
                       "CD8 T cells" = c(5),
                       "CD4 T cells" = c(7)
                       )

        annots <- do.list.switch(annots)
        names(annots) <- c("Values", "Population")
        setorderv(annots, 'Values')
        annots

    ### Add annotations

        cell.dat <- do.add.cols(cell.dat, "FlowSOM_metacluster", annots, "Values")
        cell.dat
        fwrite(cell.dat, "Annotated.data.csv")

        cell.sub <- do.add.cols(cell.sub, "FlowSOM_metacluster", annots, "Values")
        cell.sub
        fwrite(cell.dat, "Annotated.data.DR.csv")

        make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", "Population", col.type = 'factor', add.label = TRUE)
        make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", "Population", group.col, col.type = 'factor')

##########################################################################################################
#### 6. Summary data
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output - summary data")
    setwd("Output - summary data")

    ### Create summary tables

        write.sumtables(dat = cell.dat,
                        sample.col = sample.col,
                        pop.col = "Population",
                        measure.col = cellular.cols,
                        annot.col = c(group.col, batch.col),
                        group.col = group.col, do.proportions = TRUE,
                        do.mfi.per.sample = FALSE,
                        do.mfi.per.marker = TRUE)

    ### Autographs

        files <- list.files(getwd(), ".csv")
        as.matrix(files)

        files <- files[c(7,10)]
        files

        init <- fread(files[[1]])
        init
        as.matrix(names(init))

        plot.cols <- names(init)[c(6:11)]
        plot.cols

        for(i in c(1:length(files))){
          nm <- files[[i]]
          temp <- fread(files[[i]])
          nm <- gsub(".csv", "", nm)

          for(a in plot.cols){
           make.autograph(temp,
                          x.axis = group.col,
                          y.axis = a,
                          y.axis.label = a,
                          title = a,
                          subtitle = gsub("SumTable-", "", nm),
                          filename = paste0(gsub("SumTable-", "", nm), " ", a, ".png"))
          }
        }


##########################################################################################################
#### 7. Output session info and FCS files
##########################################################################################################

    ### Session info and metadata
        setwd(OutputDirectory)
        dir.create("Output - info", showWarnings = FALSE)
        setwd("Output - info")

        sink(file = "session_info.txt", append=TRUE, split=FALSE, type = c("output", "message"))
        session_info()
        sink()

        write(cellular.cols, "cellular.cols.txt")
        write(cluster.cols, "cluster.cols.txt")

    ### Write FCS files

        setwd(OutputDirectory)
        dir.create("Output - FCS files", showWarnings = FALSE)
        setwd("Output - FCS files")

        write.files(cell.dat,
                    file.prefix = exp.name,
                    divide.by = sample.col,
                    write.csv = FALSE,
                    write.fcs = TRUE)

