# Test pca.lite where principal components are added directly to dataset
test_that("add PC to dataset", {
  dat <- Spectre::demo.clustered
  
  expect_no_error(
    pca.dat <- run.pca(
      dat = dat,
      use.cols = c(11:19),
      add.pca.col = TRUE,
      pca.lite = TRUE
    )
  )
  
  # Check additional columns were added to the data
  expect_true(ncol(pca.dat) == (ncol(dat) + min(length(c(11:19)), nrow(dat))))
  
})


# Run PCA on dataset (including file creation)
test_that("run PCA", {
  dat <- Spectre::demo.clustered
  
  # Create temporary directory where files will be saved
  tmp.dir <- withr::with_tempdir(getwd(), 
                                 clean = FALSE #leaves directory alone
                                 )
  setwd(tmp.dir)
  
  expect_no_error(
    run.pca(
      dat = dat,
      use.cols = c(11:19),
      repel = TRUE
    )
  )
  
  # Test expected files were created (and are larger than 0 in size)
  # Scree plot
  expect_true(
    file.exists(grep("Scree plot", list.files(), value = TRUE)) && 
      file.size(grep("Scree plot", list.files(), value = TRUE))>0
  )
  
  # Combined PCA plot
  expect_true(
    file.exists(grep("PCA plot-combined", list.files(), value = TRUE)) && 
      file.size(grep("PCA plot-combined", list.files(), value = TRUE))>0
  )
  
  # Contribution PCA plot
  expect_true(
    file.exists(grep("PCA plot-contribution", list.files(), value = TRUE)) && 
      file.size(grep("PCA plot-contribution", list.files(), value = TRUE))>0
  )
  
  # Individuals PCA plot
  expect_true(
    file.exists(grep("PCA plot-individuals", list.files(), value = TRUE)) && 
      file.size(grep("PCA plot-individuals", list.files(), value = TRUE))>0
  )
  
  # Variables PCA plot
  expect_true(
    file.exists(grep("PCA plot-variables", list.files(), value = TRUE)) && 
      file.size(grep("PCA plot-variables", list.files(), value = TRUE))>0
  )
  
  # Contributions PCA CSV file
  expect_true(
    file.exists(grep("PCA-contributions", list.files(), value = TRUE)) && 
      file.size(grep("PCA-contributions", list.files(), value = TRUE))>0
  )
  
  # Variables PCA CSV file
  expect_true(
    file.exists(grep("PCA-variables", list.files(), value = TRUE)) && 
      file.size(grep("PCA-variables", list.files(), value = TRUE))>0
  )
  
  # Eigenvalues PCA CSV file
  expect_true(
    file.exists(grep("PCA-eigenvalues", list.files(), value = TRUE)) && 
      file.size(grep("PCA-eigenvalues", list.files(), value = TRUE))>0
  )
  
  # Individuals PCA CSV file
  expect_true(
    file.exists(grep("PCA-individuals", list.files(), value = TRUE)) && 
      file.size(grep("PCA-individuals", list.files(), value = TRUE))>0
  )
  
  # Eigen contributions PCA CSV file
  expect_true(
    file.exists(grep("PCA-eig-contrib_2-dim", list.files(), value = TRUE)) && 
      file.size(grep("PCA-eig-contrib_2-dim", list.files(), value = TRUE))>0
  )
  
})


