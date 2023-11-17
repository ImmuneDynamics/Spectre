library(stringr)

test_that("Earth mover calculation works", {
    dat <- demo.data
    markers <- names(dat)[11:19]
    
    # The samples are 05_mock_05, 11_WNV_05.
    # Looks like the latter 05 bit can be used to identify the mouse the sample came from.
    # So, we can create a pretend column based on the 05 to identify each mouse
    # and pretend that this is a paired experiment where for each mouse, we get
    # the sample from mock and wnv
    
    dat$mouse_id <- str_split_i(dat$Sample, "_", i =3)
    
    
    expect_no_error(
        do.calculate.earth.mover.dist(
            dat = dat, 
            markers = markers, 
            sample_source_col = "mouse_id", 
            batch_id_col = "Batch"
        )
    )
})
