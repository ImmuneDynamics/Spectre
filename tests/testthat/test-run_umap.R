# Use the demo data. Subsampled so it is quicker to run.
dat <- do.subsample(dat=demo.clustered, targets=50000, seed=42)
markers <- names(dat)[11:19]
dat <- run.fast.umap(dat, use.cols=markers, umap.x.name="UWOT_X", umap.y.name="UWOT_Y")

test_that("UMAP can be ran", {
    expect_true("UWOT_X" %in% names(dat))
    expect_true("UWOT_Y" %in% names(dat))
    
    expect_equal(
        nrow(na.omit(dat, cols=c("UWOT_X", "UWOT_Y"))), 
        nrow(dat)
    )
})