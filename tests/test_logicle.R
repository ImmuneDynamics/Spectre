library(Spectre)

dat = demo.asinh
use.cols = names(dat)[c(2:10)]
dat_lgcl1 = do.logicle(dat, use.cols = use.cols)
dat_lgcl2 = do.logicle(dat, use.cols = use.cols, w = NULL)

setwd("~/Documents/phd/workdir/")

make.colour.plot(dat_lgcl1,
                 x.axis = 'Ly6C_logicle',
                 y.axis = 'CD4_logicle')

make.colour.plot(dat_lgcl1,
                 x.axis = 'Ly6C',
                 y.axis = 'CD4',
                 title = 'Raw')

make.colour.plot(dat_lgcl1,
                 x.axis = 'Ly6C_asinh',
                 y.axis = 'CD4_asinh',
                 title = 'Asinh')
