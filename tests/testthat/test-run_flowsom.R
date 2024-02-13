test_that("run flowsom works", {
    dat <- Spectre::demo.asinh
    
    dat_clust <- run.flowsom(
        dat = dat,
        use.cols = names(demo.asinh)[c(11:19)],
        xdim = 14,
        ydim = 14,
        meta.k = 'auto',
        clust.seed = 42,
        meta.seed = 42,
        clust.name = "FlowSOM_cluster",
        meta.clust.name = "FlowSOM_metacluster",
        mem.ctrl = TRUE
    )
    
    expect_true("FlowSOM_cluster" %in% names(dat_clust))
    expect_true("FlowSOM_metacluster" %in% names(dat_clust))
    
})
