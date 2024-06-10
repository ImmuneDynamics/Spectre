test_that("do.rescale works", {
  dat <- Spectre::demo.asinh
  
  dat_rescaled <- do.rescale(dat, c("NK11", "CD3"))
  
  expect_true(all(
      paste0(c("NK11", "CD3"), "_rescaled") %in% names(dat_rescaled)
    ))
  
  expect_equal(min(dat_rescaled$NK11_rescaled), 0)
  expect_equal(max(dat_rescaled$NK11_rescaled), 1)
  expect_equal(min(dat_rescaled$CD3_rescaled), 0)
  expect_equal(max(dat_rescaled$CD3_rescaled), 1)
  
})
