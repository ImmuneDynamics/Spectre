
###

dat <- create.spectre(data = list('cyto' = Spectre::demo.start[,-1]))
dat

###

dat2 <- dat

for(i in 1:5){
  message(i)
  dat@data$cyto <- rbind(dat@data$cyto, dat@data$cyto)
}

###

# >100M cells

dat

dat@data$cyto[dat@data$cyto$`CellID` == 'Cell_002007',]

x <- dat@data$cyto$`CellID` == 'Cell_002007'
y <- c(1,3,5,6)

dat@data$cyto[x,]
dat@data$cyto[y,]

which(x)



dat <- Spectre2::do.asinh(dat, use.cols = c('CD3', 'CD11b'), cofactor = 500, assay = 'cyto')
dat

dat@data$cyto
dat@data$cyto_asinh


dat2 <- Spectre2::do.asinh(dat2, use.cols = c('CD3', 'CD11b'), cofactor = 500, assay = 'cyto')
dat2

dat2@data$cyto
dat2@data$cyto_asinh

