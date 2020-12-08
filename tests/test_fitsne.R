library(devtools)
install_github("immunedynamics/Spectre", ref='add_fitsne')


library(Spectre)

dat <- demo.clustered
dat.sub <- do.subsample(dat, 10000)

use.cols <- names(dat)[12:19]

dat.reduced <- run.fitsne(dat = dat.sub,
                          use.cols = use.cols)
