library(Spectre)

dat = demo.asinh
use.cols = names(dat)[c(2:10)]
dat_lgcl_auto = do.logicle(dat, use.cols = use.cols)

dat_lgcl_custom = do.logicle(dat, use.cols = use.cols, 
                             auto.infer.function = FALSE,
                             linearisation.width = 0.1)

dat_lgcl_custom = do.logicle(dat, use.cols = use.cols, 
                             auto.infer.function = FALSE,
                             max.scale.val = 10000)

dat_lgcl_custom = do.logicle(dat, use.cols = use.cols, 
                             auto.infer.function = FALSE,
                             full.transform.width = 5)

dat_lgcl_custom = do.logicle(dat, use.cols = use.cols, 
                             auto.infer.function = FALSE,
                             additional.negative.range = 2)

dat_lgcl_custom = do.logicle(dat, use.cols = use.cols, 
                             auto.infer.function = FALSE,
                             linearisation.width = 0.1,
                             max.scale.val = 10000,
                             full.transform.width = 5,
                             additional.negative.range = 2)

make.colour.plot(dat_lgcl_custom,
                 x.axis = 'Ly6C_logicle',
                 y.axis = 'CD4_logicle',
                 save.to.disk = FALSE)
make.colour.plot(dat_lgcl_auto,
                 x.axis = 'Ly6C_logicle',
                 y.axis = 'CD4_logicle',
                 save.to.disk = FALSE)

# For comparison purposes
make.colour.plot(dat_lgcl1,
                 x.axis = 'Ly6C',
                 y.axis = 'CD4',
                 title = 'Raw')

make.colour.plot(dat_lgcl1,
                 x.axis = 'Ly6C_asinh',
                 y.axis = 'CD4_asinh',
                 title = 'Asinh')
