
setwd('/Users/thomasa/Documents/GitHub/Spectre/')

devtools::document()




### Install development

devtools::install_github(repo = 'immunedynamics/spectre', ref='v2')

Spectre2::package.check()
Spectre2::package.load()




# demo.clustered
# 
# 
# demo.dat <- do.subsample(demo.clustered, 20000)
# demo.dat
# 
# table(demo.dat$FileName)
# 
# object.size(demo.dat) / 1000 / 1000 # MB
# 
# # make.colour.plot(demo.dat, 'UMAP_X', 'UMAP_Y', save.to.disk = FALSE)
# # make.multi.plot(demo.dat, 'UMAP_X', 'UMAP_Y', 'FileName', 'FileName')
# 
# 
# demo.data <- demo.dat
# save(demo.data, file = "demo.data.RData")
# 
# 
# 
# 
# 
# demo.batches.1
# object.size(demo.batches.1) / 1000 / 1000 # MB
# 
# 
# demo.batches <- do.subsample(demo.batches.1, 20000)
# demo.batches
# 
# object.size(demo.batches) / 1000 / 1000 # MB
# save(demo.batches, file = "demo.batches.RData")
