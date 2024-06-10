test_that("run phenograph", {
  dat <- Spectre::demo.asinh
  dat <- Spectre::do.subsample(dat, 5000)
  
  # Test phenograph
  expect_no_error(
    dat <- run.phenograph(dat = dat, 
                          use.cols = names(dat)[c(11:19)])
  )
  
  # Make sure "Phenograph_cluster" column was added to data
  expect_true("Phenograph_cluster" %in% colnames(dat))
  
})
