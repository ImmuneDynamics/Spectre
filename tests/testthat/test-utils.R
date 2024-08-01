test_that("add_to_data_table new columns are appended", {
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )
    
    dat_to_add <- data.table(
        marker1_asinh = rnorm(10, 1),
        marker2_asinh = rnorm(10, 2)
    )
    
    res <- add_to_data_table(dat, dat_to_add)
    
    expect_equal(res$marker1, dat$marker1)
    expect_equal(res$marker2, dat$marker2)
    expect_equal(res$marker1_asinh, dat_to_add$marker1_asinh)
    expect_equal(res$marker2_asinh, dat_to_add$marker2_asinh)
})

test_that("add_to_data_table same columns are replaced", {
    dat <- data.table(
        marker1 = rnorm(10, 1),
        marker2 = rnorm(10, 2)
    )
    
    dat_to_add <- data.table(
        marker1 = rnorm(10, 1),
        marker2_asinh = rnorm(10, 2)
    )
    
    res <- add_to_data_table(dat, dat_to_add)
    
    expect_equal(res$marker1, dat_to_add$marker1)
    expect_equal(res$marker2, dat$marker2)
    expect_equal(res$marker2_asinh, dat_to_add$marker2_asinh)
})
