test_that("make colour plot works", {
    dat <- Spectre::demo.clustered
    
    expect_no_error(
        make.colour.plot(
            dat = dat,
            x.axis = 'UMAP_X',
            y.axis = 'UMAP_Y',
            col.axis = 'Population',
            col.type = "factor",
            add.label = TRUE,
            hex = FALSE,
            hex.bins = 30,
            colours = "spectral",
            save.to.disk = FALSE,
        )
    )
    
    
    
})

test_that("fast density colour plot works", {
    dat <- Spectre::demo.clustered
    
    color_schemes <- c("jet", "spectral", "viridis", "inferno", "BuPu")
    
    for (c in color_schemes) {
        expect_no_error(
            fast.colour.plot(
                dat = dat,
                x.axis = 'UMAP_X',
                y.axis = 'UMAP_Y',
                save.to.disk = FALSE,
                colours = c
            )
        )
    }
    
    
})

test_that("fast make colour plot factor works", {
    dat <- Spectre::demo.clustered
    
    expect_no_error(
        fast.colour.plot(
            dat = dat,
            x.axis = 'UMAP_X',
            y.axis = 'UMAP_Y',
            col.axis = 'Population',
            col.type = "factor",
            add.label = TRUE,
            save.to.disk = FALSE,
        )
    )
})

test_that("fast make colour plot continuous works", {
    dat <- Spectre::demo.clustered
    
    color_schemes <- c("jet", "spectral", "viridis", "inferno", "BuPu")
    
    for (c in color_schemes) {
        expect_no_error(
            fast.colour.plot(
                dat = dat,
                x.axis = 'UMAP_X',
                y.axis = 'UMAP_Y',
                col.axis = 'NK11_asinh',
                col.type = "continuous",
                colours = c,
                save.to.disk = FALSE,
            )
        )
    }
})

test_that("fast make colour plot hex works", {
    dat <- Spectre::demo.clustered
    
    color_schemes <- c("jet", "spectral", "viridis", "inferno", "BuPu")
    
    for (c in color_schemes) {
        expect_no_error(
            fast.colour.plot(
                dat = dat,
                x.axis = 'UMAP_X',
                y.axis = 'UMAP_Y',
                col.axis = 'NK11_asinh',
                add.label = TRUE,
                hex = TRUE,
                hex.bins = 30,
                save.to.disk = FALSE,
                colours = c
            )
        )
    }
    
    
})





