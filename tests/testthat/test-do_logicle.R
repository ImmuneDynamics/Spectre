test_that("do.logicle works", {
    dat <- Spectre::demo.asinh
    expect_no_error(
        do.logicle(
            dat = dat,
            use.cols = c("NK11", "CD3", "CD45")
        )
    )
    
    dat_logicle <- do.logicle(
        dat = dat,
        use.cols = c("NK11", "CD3", "CD45")
    )
    
    expect_true(all(
        paste0(c("NK11", "CD3", "CD45"), "_logicle") %in% names(dat_logicle)
    ))
})
