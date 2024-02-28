test_that("zscore no replacement", {
  dat <- Spectre::demo.sum
  
  expect_no_error(
    dat <- do.zscore(dat = dat,
                     use.cols = names(dat)[c(4:15)])
  )
  
})

test_that("zscore replace data", {
  dat <- Spectre::demo.sum
  
  expect_no_error(
    dat <- do.zscore(dat = dat,
                     use.cols = names(dat)[c(4:15)], 
                     replace = TRUE)
  )
  
})
