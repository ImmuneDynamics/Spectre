library(Spectre)
library(data.table)

setwd("~/dropbox/spectre_paper_data_scripts/logicle_transform/")

dat = fread("mass.csv")
use.cols = names(dat)[c(1:57)]
dat_lgcl = Spectre::do.logicle(dat, use.cols = use.cols, 
                               auto.infer.function = FALSE,
                               linearisation.width = 0.5)
make.colour.plot(dat_lgcl,
                 x.axis = '176Yb_Ly6C_logicle',
                 y.axis = '168Er_CD8a_logicle',
                 save.to.disk = TRUE,
                 title = 'Logicle')

make.colour.plot(dat_lgcl,
                 x.axis = '176Yb_Ly6C_asinh',
                 y.axis = '168Er_CD8a_asinh',
                 save.to.disk = TRUE,
                 title = 'asinh')

dat = demo.asinh
use.cols = names(dat)[c(2:10)]
dat_lgcl_auto = do.logicle(dat, use.cols = use.cols)

dat_lgcl_custom = do.logicle(dat, use.cols = use.cols, 
                             auto.infer.function = FALSE,
                             linearisation.width = 1.2)

make.colour.plot(dat_lgcl_custom,
                 x.axis = 'Ly6C_logicle',
                 y.axis = 'CD8a_logicle',
                 save.to.disk = FALSE)

make.colour.plot(dat_lgcl_custom,
                 x.axis = 'Ly6C_asinh',
                 y.axis = 'CD4_asinh',
                 title = 'Asinh',
                 save.to.disk = FALSE)
