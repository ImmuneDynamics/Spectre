library(devtools)
install_github("immunedynamics/Spectre", ref='add_fitsne')


library(Spectre)

dat <- demo.clustered
dat.sub <- do.subsample(dat, 30000)

use.cols <- names(dat)[12:19]

dat.reduced <- run.fitsne(dat = dat.sub,
                          use.cols = use.cols)

make.colour.plot(dat = dat.reduced,
                 x.axis = "FItSNE_X",
                 y.axis = "FItSNE_Y",
                 col.axis = 'Population',
                 save.to.disk = TRUE)

make.colour.plot(dat = dat.reduced,
                 x.axis = "UMAP_X",
                 y.axis = "UMAP_Y",
                 col.axis = 'Population',
                 save.to.disk = TRUE)
