
library("Spectre")
package.check()
package.load()

cell.dat <- Spectre::demo.start
as.matrix(names(cell.dat))
use.cols <- names(cell.dat)[c(2:10)]

cell.dat.1 <- do.asinh(cell.dat, use.cols, cofactor = 1000)
cell.dat.2 <- do.mn(cell.dat, use.cols, asinh.cf = 500, lower.threshold = -2)

make.colour.plot(cell.dat.1, "CD3_asinh", "CD4_asinh", filename = "for test.png")




z.dat <- do.zscore(cell.dat.1, paste0(use.cols, "_asinh"))
z.dat

make.colour.plot(z.dat, "CD11b_asinh", "CD45_asinh")
make.colour.plot(z.dat, "CD3_asinh", "CD4_asinh")

make.colour.plot(cell.dat.2, "CD3_transf", "CD4_transf")








dat.dt <- as.data.table(Spectre::demo.clustered)
dat.df <- as.data.frame(Spectre::demo.clustered)
dat.mt <- as.matrix(Spectre::demo.clustered)

object.size(dat.dt) / 1024 / 1024
object.size(dat.df) / 1024 / 1024
object.size(dat.mt) / 1024 / 1024


x <- dat.df[1:100, ]




dat.transf <- do.mn(cell.dat, use.cols, asinh.cf = 1000)
dat.transf


quantile(cell.dat.1$CD3_asinh, 0.005)

make.colour.plot(cell.dat.1, "CD3_asinh", "CD4_asinh", "CD3_asinh")
make.colour.plot(dat.transf, "CD3_transf", "CD4_transf", "CD3_transf")
make.colour.plot(dat.transf, "CD11b_transf", "CD45_transf")

dat.transf <- run.flowsom(dat.transf, paste0(use.cols, "_transf"))
dat.transf <- run.fitsne(dat.transf, paste0(use.cols, "_transf"))

make.colour.plot(dat.transf, "FItSNE_X", "FItSNE_Y", "FlowSOM_metacluster", "factor", add.label = TRUE)
