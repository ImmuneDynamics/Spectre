test_that("default SpectreObject creation works", {
    dt = data.table(x=seq(1,10), y=seq(1,10))
    
    test_obj = new("SpectreObject", cytometry_data=dt, citeseq_data=dt)
    
    expect_equal(ncol(test_obj@cytometry_data), 2)
    expect_equal(nrow(test_obj@cytometry_data), 10)
    expect_equal(ncol(test_obj@citeseq_data), 2)
    expect_equal(nrow(test_obj@citeseq_data), 10)
    
})

test_that("default NULL slots for SpectreObject creation works", {
    expect_no_error(
        new("SpectreObject")
    )
})

test_that("missing either of the 2 slots for SpectreObject creation errors", {
    dt = data.table(x=seq(1,10), y=seq(1,10))
    
    expect_error(
        new("SpectreObject", cytometry_data=dt)
    )
    
    expect_error(
        new("SpectreObject", citeseq_data=dt)
    )
})
