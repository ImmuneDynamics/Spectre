test_that("MDS plot works", {
  dat <- demo.data
  markers <- names(dat)[11:19]
  expect_no_error(
      make.mds.plot(
          dat = dat, 
          sample_col = "Sample", 
          markers = markers, 
          colour_by = "Sample")
    )
  
  expect_no_error(
      make.mds.plot(
          dat = dat, 
          sample_col = "Sample", 
          markers = markers, 
          colour_by = "Batch")
  )
  
  expect_no_error(
      make.mds.plot(
          dat = dat, 
          sample_col = "Sample", 
          markers = markers, 
          colour_by = "Batch",
          add_point_label = FALSE)
  )
})
