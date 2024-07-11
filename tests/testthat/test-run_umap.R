# Use the demo data. Subsampled so it is quicker to run.

test_that("Fast UMAP single thread can run", {
    dat <- do.subsample(
        dat = demo.data, 
        targets = 1000, 
        seed = 42
    )
    markers <- names(dat)[11:19]
    
    dat <- run.umap(
        dat, 
        use.cols = markers, 
        umap.x.name = "UWOT_X", 
        umap.y.name = "UWOT_Y",
        fast = TRUE,
        verbose = FALSE,
        n_threads = 1
    )
    
    expect_true("UWOT_X" %in% names(dat))
    expect_true("UWOT_Y" %in% names(dat))

    expect_equal(
        nrow(na.omit(dat, cols = c("UWOT_X", "UWOT_Y"))),
        nrow(dat)
    )
})

test_that("Fast UMAP infer n_threads can run", {
    dat <- do.subsample(
        dat = demo.data, 
        targets = 1000, 
        seed = 42
    )
    markers <- names(dat)[11:19]
    
    dat <- run.umap(
        dat, 
        use.cols = markers, 
        umap.x.name = "UWOT_X", 
        umap.y.name = "UWOT_Y",
        fast = TRUE,
        verbose = FALSE,
        n_threads = 'auto'
    )
    
    expect_true("UWOT_X" %in% names(dat))
    expect_true("UWOT_Y" %in% names(dat))
    
    expect_equal(
        nrow(na.omit(dat, cols = c("UWOT_X", "UWOT_Y"))),
        nrow(dat)
    )
})

test_that("Default UMAP can run", {
    dat <- do.subsample(
        dat = demo.data, 
        targets = 1000, 
        seed = 42
    )
    markers <- names(dat)[11:19]
    
    dat <- run.umap(
        dat, 
        use.cols = markers, 
        umap.x.name = "UMAP_X_2", 
        umap.y.name = "UMAP_Y_2",
        fast = FALSE,
        verbose = FALSE
    )
    
    expect_true("UMAP_X_2" %in% names(dat))
    expect_true("UMAP_Y_2" %in% names(dat))
    
    expect_equal(
        nrow(na.omit(dat, cols = c("UMAP_X_2", "UMAP_Y_2"))),
        nrow(dat)
    )
})



