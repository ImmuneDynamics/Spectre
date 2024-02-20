test_that("make volcano plot", {
  # Create temporary directory where files will be saved
  tmp.dir <- withr::with_tempdir(getwd(), 
                                 clean = FALSE #leaves directory alone
                                 )
  setwd(tmp.dir)
  
  # Prepare data
  sum.dat <- Spectre::demo.sum
  to.plot <- names(sum.dat)[c(4:15)]
  
  # Calculate statistics
  expect_no_error(
    stat.dat <- create.stats(dat = sum.dat,
                             use.cols = to.plot,
                             sample.col = 'Sample',
                             group.col = 'Group',
                             comparisons = list(c('Mock', 'WNV'))
    )
  )
  
  # Create volcano plot
  expect_no_error(
    make.volcano.plot(dat.p = stat.dat[2,..to.plot],
                      dat.fc = stat.dat[1,..to.plot],
                      vars = to.plot,
                      title = 'Volcano')
  )
  
  # Test expected file was created (and is larger than 0 in size)
  expect_true(
    file.exists(grep("Volcano plot", list.files(), value = TRUE)) && 
      file.size(grep("Volcano plot", list.files(), value = TRUE))>0
  )
  
})
