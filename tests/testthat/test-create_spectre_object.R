test_that("default SpectreObject creation works", {
    dt = data.table(x=seq(1,10), y=seq(1,10))
    
    test_obj = SpectreObject(cytometry_data=dt, citeseq_data=dt)
    
    expect_equal(ncol(test_obj@cytometry_data), 2)
    expect_equal(nrow(test_obj@cytometry_data), 10)
    expect_equal(ncol(test_obj@citeseq_data), 2)
    expect_equal(nrow(test_obj@citeseq_data), 10)
    
})

test_that("no cytometry_data slot for SpectreObject creation errors", {
    expect_error(
        SpectreObject()
    )
})

test_that("no citeseq_data slot for SpectreObject creation works", {
    dt = data.table(x=seq(1,10), y=seq(1,10))
    expect_no_error(
        SpectreObject(cytometry_data=dt)
    )
})
