##########################################################################################################
#### Spectre: General Discovery Workflow - (2/4) - Batch alignment
##########################################################################################################

    # Spectre R package: https://github.com/ImmuneDynamics/Spectre
    # Thomas Myles Ashhurst, Felix Marsh-Wakefield, Givanna Putri

##########################################################################################################
#### Analysis session setup
##########################################################################################################

    ### Load packages

        library(Spectre)
        Spectre::package.check()    # Check that all required packages are installed
        Spectre::package.load()     # Load required packages

    ### Set DT threads

        getDTthreads()

    ### Set primary directory

        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        getwd()
        PrimaryDirectory <- getwd()
        PrimaryDirectory

    ### Set input directory

        setwd(PrimaryDirectory)
        setwd("Output 1 - data prep/Output 1.2 - transformed data/")
        InputDirectory <- getwd()
        InputDirectory
        setwd(PrimaryDirectory)
        
    ### Set metadata directory
        
        setwd(PrimaryDirectory)
        setwd("metadata/")
        MetaDirectory <- getwd()
        MetaDirectory
        setwd(PrimaryDirectory)
        
    ### Set output directory

        setwd(PrimaryDirectory)
        dir.create("Output 2 - batch alignment", showWarnings = FALSE)
        setwd("Output 2 - batch alignment")
        OutputDirectory <- getwd()
        
        setwd(PrimaryDirectory)

##########################################################################################################
#### Import data
##########################################################################################################

    ### Read in cellular data
        
        setwd(InputDirectory)
        
        list.files(getwd(), ".csv")
        
        cell.dat <- fread("cell.dat.csv")
        cell.dat
        
    ### Read in metadata
        
        setwd(MetaDirectory)
        
        list.files(getwd(), ".csv")
        
        meta.dat <- fread("sample.details.csv")
        meta.dat
        
##########################################################################################################
#### Define columns
##########################################################################################################
        
    setwd(InputDirectory)
        
    ### Define other cols
        
        as.matrix(names(cell.dat))
        
        sample.col <- 'Sample'
        group.col <- 'Group'
        batch.col <- 'Batch'   

    ### Define cellular cols
        
        as.matrix(names(cell.dat))
        
        cellular.cols <- names(cell.dat)[c(11:18)]
        as.matrix(cellular.cols)
        
    ### Define clustering cols
        
        as.matrix(names(cell.dat))
        
        cluster.cols <- names(cell.dat)[c(11:18)]
        as.matrix(cluster.cols)
        
    ### Define columns to align
        
        as.matrix(names(cell.dat))
        
        to.align <- names(cell.dat)[c(11:18)]
        as.matrix(to.align)  

    ### Double check selections
        
        sample.col
        group.col
        batch.col
        
        as.matrix(cellular.cols)
        as.matrix(cluster.cols)
        as.matrix(to.align)
    
##########################################################################################################
#### Initial (pre-alignment) plots
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 2.1 - initial plots")
    setwd("Output 2.1 - initial plots")
        
    ### Pre-alignment UMAP
        
        rm(sub)
        sub <- do.subsample(cell.dat, 10000)
        sub <- run.umap(sub, cluster.cols)
        
    ### Create plots
        
        make.colour.plot(sub, "UMAP_X", "UMAP_Y", batch.col, col.type = 'factor', filename = paste0('Batches.png'))
        make.colour.plot(sub, "UMAP_X", "UMAP_Y", group.col, col.type = 'factor', filename = paste0('Groups.png'))
        
        make.multi.plot(sub, "UMAP_X", "UMAP_Y", cellular.cols, figure.title = 'Celluar markers')
        make.multi.plot(sub, "UMAP_X", "UMAP_Y", cluster.cols, figure.title = 'Clustering markers')
        make.multi.plot(sub, "UMAP_X", "UMAP_Y", to.align, figure.title = 'Markers to align')
        
    ### Cleanup
        
        rm(sub)
        
##########################################################################################################
#### Define and assess reference samples
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 2.2 - reference samples")
    setwd("Output 2.2 - reference samples")
        
    ### Define 'reference' samples
        
        meta.dat
    
        as.matrix(unique(cell.dat[[sample.col]]))
    
        ref.ctrls <- unique(cell.dat[[sample.col]])[c(1,5)]
        ref.ctrls
        
        ref.dat <- do.filter(cell.dat, use.col = sample.col, values = ref.ctrls)
        ref.dat
        
    ### Check reference in ref.dat and cell.dat
        
        unique(ref.dat[[batch.col]])
        unique(cell.dat[[batch.col]])
        
        unique(ref.dat[[sample.col]])
        unique(cell.dat[[sample.col]])
        
    ### UMAP of reference samples
        
        rm(sub)
        sub <- do.subsample(ref.dat, 10000)
        sub <- run.umap(sub, cluster.cols)
        
    ### Create plots
        
        make.colour.plot(sub, "UMAP_X", "UMAP_Y", batch.col, col.type = 'factor', filename = paste0('Batches.png'))
        make.colour.plot(sub, "UMAP_X", "UMAP_Y", group.col, col.type = 'factor', filename = paste0('Groups.png'))
        
        make.multi.plot(sub, "UMAP_X", "UMAP_Y", cellular.cols, figure.title = 'Celluar markers')
        make.multi.plot(sub, "UMAP_X", "UMAP_Y", cluster.cols, figure.title = 'Clustering markers')
        make.multi.plot(sub, "UMAP_X", "UMAP_Y", to.align, figure.title = 'Markers to align')
        
    ### Cleanup
        
        rm(sub)
        
##########################################################################################################
#### Perform and validate COARSE alignment (to improve FlowSOM clustering during fine alignment)
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 2.3 - coarse alignment")
    setwd("Output 2.3 - coarse alignment")
    
    crs.dir <- getwd()
        
    ### Setup
    
        crs.dat <- cell.dat
    
        cellular.cols
        cluster.cols
        to.align
        
        method <- '95p'
        crs.append <- '_coarseAlign'
    
    ### Coarse alignment
        
        setwd(crs.dir)
        
        crs.dat <- run.align(ref.dat = ref.dat, 
                                target.dat = crs.dat, 
                                batch.col = batch.col, 
                                align.cols = to.align, 
                                method = method, 
                                append.name = crs.append, 
                                dir = crs.dir)

        setwd(crs.dir)
        gc()
        
    ### Check results
        
        as.matrix(names(crs.dat))
        crs.dat
        
    ### Plot alignment results
        
        setwd(crs.dir)
        dir.create("A - Alignment plots")
        setwd("A - Alignment plots")
        
        rm(sub)
        sub  <- do.subsample(crs.dat, 100000)
        
        for(i in to.align){
            make.colour.plot(sub, paste0(i, crs.append), i, batch.col)
        }
        
        rm(sub)

    ### Examine clustering
        
        setwd(crs.dir)
        dir.create("B - Clustering check")
        setwd("B - Clustering check")
        
        rm(sub)
        sub <- do.subsample(crs.dat, 10000)
        sub <- run.umap(sub, paste0(cluster.cols, crs.append))
        
        make.colour.plot(sub, "UMAP_X", "UMAP_Y", batch.col, col.type = 'factor', filename = paste0('Batches ', method, ".png"))
        make.multi.plot(sub, "UMAP_X", "UMAP_Y", cellular.cols, figure.title = paste0("Markers - raw - ", method))
        make.multi.plot(sub, "UMAP_X", "UMAP_Y", paste0(cellular.cols, crs.append), figure.title = paste0("Markers - aligned - ", method))
 
        rm(sub)
 
    ### Finalise data and assign as cell.dat
        
        cell.dat <- crs.dat
        rm(crs.dat)
               
    ### Write coarse aligned data
        
        setwd(crs.dir)
        dir.create("C - Coarse aligned data")
        setwd("C - Coarse aligned data")
        
        fwrite(cell.dat, 'cell.dat.csv')
        
        write.files(cell.dat, 
                    file.prefix = "Coarse_aligned", 
                    divide.by = sample.col, 
                    write.csv = FALSE, 
                    write.fcs = TRUE)
        
##########################################################################################################
#### Fine alignment with CytoNorm
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output 2.4 - fine alignment")
    setwd("Output 2.4 - fine alignment")
    
    fine.dir <- getwd()

    ### Settings
    
        fine.dat <- cell.dat
    
        cytonorm.goal <- 'mean'
        cytonorm.nQ <- 101
        
        fine.append <- '_fineAlign'
    
    ### Re-extract reference data (from newly coarse aligned cell.dat)
    
        rm(ref.dat)
        
        ref.dat <- do.filter(fine.dat, use.col = sample.col, values = ref.ctrls)
        ref.dat
        
        unique(ref.dat[[sample.col]])
        unique(ref.dat[[batch.col]])
        
    ### Prep FlowSOM and perform clustering using ref.dat

        setwd(fine.dir)
        align.model <- prep.cytonorm(dat = ref.dat,
                                     cellular.cols = c(cellular.cols, paste0(cellular.cols, crs.append)),
                                     cluster.cols = paste0(cluster.cols, crs.append),
                                     batch.col = batch.col,
                                     dir = fine.dir,
                                     xdim = 14,
                                     ydim = 14,
                                     meta.k = 8)

        str(align.model, 1)

    ### Examine clustering and generate plots

        setwd(fine.dir)
        dir.create("A - Examine referece data clustering")
        setwd("A - Examine referece data clustering")

        sub <- do.subsample(align.model$dt, 10000)
        sub <- run.umap(sub, paste0(cluster.cols, crs.append))

        make.colour.plot(sub, "UMAP_X", "UMAP_Y", "File", col.type = 'factor', filename = "Reference data - batches.png")
        make.colour.plot(sub, "UMAP_X", "UMAP_Y", "prep.fsom.metacluster", col.type = 'factor', add.label = TRUE, filename = "Reference data - metaclusters.png")
       
        make.multi.plot(sub, "UMAP_X", "UMAP_Y", paste0(cellular.cols, crs.append), figure.title = "Reference data - markers - coarse aligned")
        make.multi.plot(sub, "UMAP_X", "UMAP_Y", cellular.cols, figure.title = "Reference data - markers - raw")

        rm(sub)

    ### Train the alignment conversions in the 'align.model' object

        setwd(fine.dir)
        align.model <- train.cytonorm(model = align.model,
                                      align.cols = to.align,
                                      cytonorm.goal = cytonorm.goal,
                                      cytonorm.nQ = cytonorm.nQ,
                                      dir = fine.dir)

        str(align.model, 1)
        gc()

    ### Run cytonorm

        setwd(fine.dir)
        fine.dat <- run.cytonorm(dat = fine.dat,
                                 model = align.model,
                                 batch.col = batch.col,
                                 append.name = fine.append,
                                 dir = fine.dir)

        as.matrix(names(fine.dat))
        fine.dat

    ### Examine results
        
        setwd(fine.dir)
        dir.create("B - Examine cytonorm results")
        setwd("B - Examine cytonorm results")
        
        res <- do.subsample(fine.dat, 10000)
        res <- run.umap(res, paste0(cluster.cols, "_fineAlign"))
        
        make.colour.plot(res, "UMAP_X", "UMAP_Y", batch.col, col.type = 'factor', filename = "Target data - batches.png")
        make.colour.plot(res, "UMAP_X", "UMAP_Y", group.col, col.type = 'factor', filename = "Target data - groups.png")
        
        make.colour.plot(res, "UMAP_X", "UMAP_Y", paste0("Alignment_MC", fine.append), col.type = 'factor', add.label = TRUE, filename = "Target data - metaclusters.png")
        
        make.multi.plot(res, "UMAP_X", "UMAP_Y", paste0(cellular.cols, crs.append), figure.title = "Target data - markers - coarse aligned")
        make.multi.plot(res, "UMAP_X", "UMAP_Y", paste0(cellular.cols, fine.append), figure.title = "Target data - markers - fine aligned")
        make.multi.plot(res, "UMAP_X", "UMAP_Y", cellular.cols, figure.title = "Target data - markers - raw")
        
        rm(res)

    ### Plot comparisons

        sub <- do.subsample(fine.dat, 10000)
        
        crs.append
        fine.append
        
        ## Raw vs coarse
        
            setwd(fine.dir)
            dir.create("C - Comparison plots - raw vs coarse")
            setwd("C - Comparison plots - raw vs coarse")
            
            for(i in cellular.cols){
                a <- paste0(i, crs.append)
                make.colour.plot(sub, a, i, batch.col)
            }
        
        ## Raw vs fine
            
            setwd(fine.dir)
            dir.create("D - Comparison plots - raw vs fine")
            setwd("D - Comparison plots - raw vs fine")

            for(i in cellular.cols){
              a <- paste0(i, fine.append)
              make.colour.plot(sub, a, i, batch.col)
            }

        ## Coarse vs fine    
            
            setwd(fine.dir)
            dir.create("E - Comparison plots - coarse vs fine")
            setwd("E - Comparison plots - coarse vs fine")
    
            for(i in cellular.cols){
              a <- paste0(i, fine.append)
              o <- paste0(i, crs.append)
              make.colour.plot(sub, a, o, batch.col)
            }

    ### Finalsie data
        
        cell.dat <- fine.dat
        rm(fine.dat)
        
    ### Save initial data

        setwd(fine.dir)
        dir.create("F - Fine aligned data")
        setwd("F - Fine aligned data")

        fwrite(cell.dat, "cell.dat.csv")
        
        write.files(cell.dat, 
                    file.prefix = "Fine_aligned", 
                    divide.by = sample.col, 
                    write.csv = FALSE, 
                    write.fcs = TRUE)
        
##########################################################################################################
#### Save session info
##########################################################################################################

    setwd(OutputDirectory)
    dir.create("Output - info", showWarnings = FALSE)
    setwd("Output - info")

    ### Save session info to disk
        
        setwd(OutputDirectory)
        dir.create("Output - info", showWarnings = FALSE)
        setwd("Output - info")
        
        sink(file = "session_info.txt", append=TRUE, split=FALSE, type = c("output", "message"))
        session_info()
        sink()
        