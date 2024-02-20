test_that("aggregation works", {
  dat <- Spectre::demo.clustered
  
  expect_no_error(
    do.aggregate(
      dat = dat,
      use.cols = names(dat)[c(11:19)],
      by = "Population"
    )
  )
  
})


