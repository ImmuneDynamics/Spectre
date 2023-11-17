
test_that("FITSNE runs", {
    # Use the demo data. Subsampled so it is quicker to run.
    dat <- do.subsample(dat=demo.data, targets=1000, seed=42)
    markers <- names(dat)[11:19]
    dat <- run.fitsne(dat, use.cols=markers)
    
    expect_true("FItSNE_X" %in% names(dat))
    expect_true("FItSNE_Y" %in% names(dat))
    
    expect_equal(
        nrow(na.omit(dat, cols=c("FItSNE_X", "FItSNE_Y"))), 
        nrow(dat)
    )
})
