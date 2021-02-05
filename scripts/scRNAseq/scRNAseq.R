##########################################################################################################
#### Spectre - Split into FCS files
##########################################################################################################
    
    # Spectre R package: https://sydneycytometry.org.au/spectre
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
        setwd("filtered_gene_bc_matrices/hg19/")
        InputDirectory <- getwd()
        InputDirectory
        setwd(PrimaryDirectory)
    
    ### Set output directory
    
        setwd(PrimaryDirectory)
        dir.create("Spectre_scRNAseq_output", showWarnings = FALSE)
        setwd("Spectre_scRNAseq_output")
        OutputDirectory <- getwd()
        setwd(PrimaryDirectory)
        
##########################################################################################################
#### Read in scRNAseq data and prepare with Seurat
##########################################################################################################     
   
    ### Code derived from Seurat's excellent analysis tutorial
        
        library(dplyr)
        library(Seurat)
        library(patchwork)
        
        setwd(InputDirectory)
        
        pbmc.data <- Read10X(data.dir = getwd())
        pbmc.data
        
        pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
        pbmc
        
    ### Exploration
        
        pbmc$orig.ident
        pbmc$nCount_RNA
        pbmc$nFeature_RNA
        
        nGenes <- pbmc@assays$RNA@counts@Dim[[1]]
        nCells <- pbmc@assays$RNA@counts@Dim[[2]]
        
        nObs <- length(pbmc@assays$RNA@counts@x)
        
        nGenes * nCells
        
        nObs / nGenes
        
    ### QC etc
        
        # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
        pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
        
        # Visualize QC metrics as a violin plot
        VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
        
        # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
        # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
        
        plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
        plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
        plot1 + plot2
        
        pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
        
        pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
        pbmc <- NormalizeData(pbmc)
        #pbmc[["RNA"]]@data
        
        pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
        
        # Identify the 10 most highly variable genes
        top10 <- head(VariableFeatures(pbmc), 10)
        
        # plot variable features with and without labels
        plot1 <- VariableFeaturePlot(pbmc)
        plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
        plot1 + plot2
        
        all.genes <- rownames(pbmc)
        pbmc <- ScaleData(pbmc, features = all.genes)
        
    ### Variable genes
        var.genes <- pbmc[["RNA"]]@var.features
        var.genes
 
##########################################################################################################
#### Clustering & DR with seurat
##########################################################################################################     

    ### PCA
        
        pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
        
        # Examine and visualize PCA results a few different ways
        print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
        
        DimPlot(pbmc, reduction = "pca")
        
    ### Clustering
        
        pbmc <- FindNeighbors(pbmc, dims = 1:10)
        pbmc <- FindClusters(pbmc, resolution = 0.5)
        
    ### UMAP
        
        # If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
        # 'umap-learn')
        pbmc <- RunUMAP(pbmc, dims = 1:10)
        
        # note that you can set `label = TRUE` or use the LabelClusters function to help label
        # individual clusters
        DimPlot(pbmc, reduction = "umap")
        
##########################################################################################################
#### Convert to data.table
##########################################################################################################     

    ### Convert from Seurat object to data.table
        
        pbmc.dt <- create.dt(pbmc)
        
    ### Review
        
        str(pbmc.dt, 1)
        str(pbmc.dt$data.table)

##########################################################################################################
#### Plotting etc
##########################################################################################################     

    ### Plotting
        
        pbmc.dt$geneNames
        pbmc.dt$meta.data
        pbmc.dt$assays
        pbmc.dt$slots
        pbmc.dt$dim.reds
        
        make.colour.plot(pbmc.dt$data.table, 'umap_UMAP_1', 'umap_UMAP_2')
        
        

        