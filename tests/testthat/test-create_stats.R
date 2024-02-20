test_that("run statistics", {
  dat <- Spectre::demo.sum
  
  expect_no_error(
    create.stats(dat = dat,
                 use.cols = names(dat)[c(4:15)],
                 sample.col = "Sample",
                 group.col = "Group",
                 comparisons = list(c('Mock', 'WNV'))
    )
  )
  
})
