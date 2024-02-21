test_that("make multi plot", {
  # Create temporary directory where files will be saved
  tmp.dir <- withr::with_tempdir(getwd(), 
                                 clean = FALSE #leaves directory alone
  )
  setwd(tmp.dir)
  
  # Create multi-plot
  dat <- Spectre::demo.clustered
  colour_schemes <- c("jet", "spectral", "viridis", "inferno", "magma", "BuPu")
  
  for(c in colour_schemes) {
    expect_no_error(
      Spectre::make.multi.plot(dat = Spectre::demo.clustered,
                               x.axis = "UMAP_X",
                               y.axis = "UMAP_Y",
                               plot.by = c("Ly6C_asinh", "B220_asinh", "CD45_asinh"), 
                               colours = c)
    )
    
    # Test expected file was created (and is larger than 0 in size)
    expect_true(
      file.exists(grep("Multi plot", list.files(), value = TRUE)) && 
        file.size(grep("Multi plot", list.files(), value = TRUE))>0
    )
  }
  
})

test_that("fast make multi plot", {
  # Create temporary directory where files will be saved
  tmp.dir <- withr::with_tempdir(getwd(), 
                                 clean = FALSE #leaves directory alone
  )
  setwd(tmp.dir)
  
  # Create multi-plot
  dat <- Spectre::demo.clustered
  colour_schemes <- c("jet", "spectral", "viridis", "inferno", "magma", "BuPu")
  
  for(c in colour_schemes) {
    expect_no_error(
      Spectre::fast.multi.plot(dat = Spectre::demo.clustered,
                               x.axis = "UMAP_X",
                               y.axis = "UMAP_Y",
                               plot.by = c("Ly6C_asinh", "B220_asinh", "CD45_asinh"), 
                               colours = c)
      )
    
    # Test expected file was created (and is larger than 0 in size)
    expect_true(
      file.exists(grep("Multi plot", list.files(), value = TRUE)) && 
        file.size(grep("Multi plot", list.files(), value = TRUE))>0
    )
  }
  
})

