#' Unit test for creating data.table from other popular SC based objects.

context("Test create data.table from objects")
library(Spectre)
# Install from https://github.com/satijalab/seurat-data
library(SeuratData)
library(Seurat)

test_that("Seurat RNA only object is converted", {
    InstallData("pbmc3k")
    pbmc <- LoadData("pbmc3k", type = "pbmc3k.final")
    
    dat <- create.dt(pbmc)
    
    # Meant to return a list
    expect_equal(typeof(dat), "list")
    
    # Check all the genes are captured
    gene_names <- paste0(rownames(pbmc), "_RNA_counts")
    expect_true(all(gene_names %in% colnames(dat[['data.table']])))
    # There is an element in the list that stores the gene names
    expect_true(all(rownames(pbmc) %in% dat[['geneNames']]))
    
    # Check all cell barcodes are captured
    expect_true(all(colnames(pbmc) %in% dat[['data.table']][['cellNames']]))
    expect_true(all(colnames(pbmc) %in% dat[['cellNames']]))
    
    # Check metadata
    expect_true(all(colnames(pbmc@meta.data) %in% dat[['meta.data']]))
    for (meta_dat in colnames(pbmc@meta.data)) {
        expected <- pbmc@meta.data[[meta_dat]]
        obj <- dat[['data.table']][[meta_dat]]
        expect_true(all.equal(obj, expected))
    }
    
    # Check HVG
    expect_true(all(VariableFeatures(pbmc) %in% dat[['var.features']]))
    expect_true(all(head(VariableFeatures(pbmc), 10) %in% dat[['var.features.top10']]))
    
    # Check assay
    expect_true("_RNA" == dat[['assays']])
    
    # Check all slots are captured
    expect_true(all(c("_counts", "_data", "_scale.data") %in% dat[['slots']]))
    
    # Check dimreds
    expect_true(all(c("pca_", "umap_") %in% dat[['dim.reds']]))
    # By default the PBMC data has 50 PCs.
    pca_coordinates <- Embeddings(pbmc, reduction = "pca")
    for (i in c(1:50)) {
        expected <- pca_coordinates[,paste0('PC_', i) ]
        obj <- dat[['data.table']][[paste0('pca_PC_', i)]]
        names(obj) <- dat[['data.table']][['cellNames']]
        expect_true(identical(expected, obj))
    }
    
    # By default, the PBMC data has 2 UMAP
    umap_coordinates <- Embeddings(pbmc, reduction = "umap")
    for (i in c(1:2)) {
        expected <- umap_coordinates[,paste0('UMAP_', i) ]
        obj <- dat[['data.table']][[paste0('umap_UMAP_', i)]]
        names(obj) <- dat[['data.table']][['cellNames']]
        expect_true(identical(expected, obj))
    }
    
    
})