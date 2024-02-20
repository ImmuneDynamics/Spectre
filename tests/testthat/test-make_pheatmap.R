# Test marker vs cluster heatmap
test_that("make cluster vs marker heatmap", {
  dat <- Spectre::demo.exp
  
  # Create temporary directory where files will be saved
  tmp.dir <- withr::with_tempdir(getwd(), 
                                 clean = FALSE #leaves directory alone
                                 )
  setwd(tmp.dir)
  
  expect_no_error(
    make.pheatmap(dat = dat,
                  file.name = "Expression pheatmap.png",
                  plot.title = "Expression",
                  sample.col = "Population",
                  plot.cols = names(dat)[c(2:10)])
    )
  
  # Test expected file was created (and is larger than 0 in size)
  expect_true(
    file.exists(grep("Expression pheatmap", list.files(), value = TRUE)) && 
      file.size(grep("Expression pheatmap", list.files(), value = TRUE))>0
  )
  
})

# Test z-score heatmap
test_that("make z-score heatmap", {
  z.dat <- do.zscore(dat = Spectre::demo.sum,
                     use.cols = names(Spectre::demo.sum)[c(4:15)],
                     replace = TRUE)
  
  # Create temporary directory where files will be saved
  tmp.dir <- withr::with_tempdir(getwd(), 
                                 clean = FALSE #leaves directory alone
                                 )
  setwd(tmp.dir)
  
  expect_no_error(
    make.pheatmap(dat = z.dat,
                  file.name = "z-score.png",
                  plot.title = "z-score",
                  sample.col = "Sample",
                  plot.cols = names(z.dat)[c(4:15)],
                  annot.cols = names(z.dat)[c(2:3)],
                  is.fold = TRUE,
                  fold.range = c(3,-3)
                  )
    )
  
  # Test expected file was created (and is larger than 0 in size)
  expect_true(
    file.exists(grep("z-score", list.files(), value = TRUE)) && 
      file.size(grep("z-score", list.files(), value = TRUE))>0
  )
  
})



