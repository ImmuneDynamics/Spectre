# Create small test data
dat <- do.subsample(
    dat = demo.clustered, 
    targets = 50000, 
    seed = 42
)
markers <- names(dat)[11:19]

test_that("make.bubble.plot works if all parameters are in order", {
    expect_warning(
        plt <- make.bubble.plot(
            dt = dat, 
            markers_to_plot=markers, 
            cluster = "FlowSOM_metacluster"
        )
    , regexp=NA)
    
    # To make sure the plot is drawn.
    expect_true("ggplot" %in% class(plt))
})

test_that(
    "make.bubble.plot throws warning when using col.scheme in the viridis package", 
    {
        col.scheme = "inferno"
        expect_warning(
            plt <- make.bubble.plot(
                dt = dat, 
                markers_to_plot=markers, 
                cluster = "FlowSOM_metacluster",
                col.scheme=col.scheme
            )
            , regexp=paste("Finding col.scheme", col.scheme, "in viridis package."))
        # To make sure the plot is drawn.
        expect_true("ggplot" %in% class(plt))
    }
)

test_that(
    "make.bubble.plot still works when using col.scheme not in either the viridis or RColorBrewer package", 
    {
        suppressWarnings(plt <- make.bubble.plot(
            dt = dat, 
            markers_to_plot=markers, 
            cluster = "FlowSOM_metacluster",
            col.scheme = "random"
        ))
        # To make sure the plot is drawn.
        expect_true("ggplot" %in% class(plt))
    }
)
