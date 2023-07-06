

###

library('Spectre2')

dat <- new('Spectre')
dat

dat <- Spectre::demo.start
cols <- names(dat)[-1]
dat <- do.asinh(dat = dat, use.cols = cols, cofactor = 500)
dat

###

test <- Spectre::demo.start
test$FileName <- gsub('.csv', '', test$FileName)
# z <- paste0("%0", nchar(nrow(test)), "d")
# x <- sprintf(z,c(1:nrow(test)))
# test$CellID <- paste0(test$FileName, '_', x)
# test <- test[,c(which(names(test) == 'CellID'), which(names(test) != 'CellID')), with = FALSE]
# setkey(test, CellID)
# 
# key(test)
test
# dat <- new('Spectre')
# dat@data$meta <- test[,c(1:2)]
# dat@data$cyto <- test[,-c(2)]
# dat@key <- 'CellID'

library('Spectre2')

dat <- create.spectre(data = list('meta' = test[,1], 
                                  'cyto' = test[,-1]), 
                      name = 'Test')
dat

# dat <- create.spectre(data = list('meta' = test[,1], 
#                                   'cyto' = test[,-1]), 
#                       name = 'Test', 
#                       key = test$FileName)
# dat

dat <- do.asinh(dat = dat, assay = 'cyto')
dat



