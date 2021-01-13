library(Spectre)
Spectre::package.check()    # Check that all required packages are installed
Spectre::package.load()     # Load required packages

## Set PrimaryDirectory
dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
PrimaryDirectory <- getwd()

## Create output directory
dir.create("Output_Spectre", showWarnings = FALSE)
setwd("Output_Spectre")
OutputDirectory <- getwd()

## Start test
setwd(PrimaryDirectory)
cell.dat <- data.table::fread("DimRed_HexPlotdemo.csv")
sample.col <- 'Sample'
group.col <- 'Group'
ColumnNames <- as.matrix(unname(colnames(cell.dat))) # assign reporter and marker names (column names) to 'ColumnNames'
CellularColsNos <- c(5,6,8,9,11:13)
CellularCols <- ColumnNames[CellularColsNos]

cell.dat.sub <- Spectre::do.subsample(dat = cell.dat, 
                                      method = 'random',
                                      samp.col = sample.col,
                                      targets = 1000)

setwd(OutputDirectory)
setwd('Output-plots')
dir.create("hexPlots")
setwd("hexPlots")

## The following two must be the same!
make.colour.plot(dat = cell.dat,
                 x.axis = "UMAP_X",
                 y.axis = "UMAP_Y",
                 col.axis = ColumnNames[30],
                 title = 'HexDemo_Bins=100',
                 hex = TRUE,
                 hex.bins = 100)

make.colour.plot(dat = cell.dat,
                 x.axis = "UMAP_X",
                 y.axis = "UMAP_Y",
                 col.axis = ColumnNames[30],
                 title = 'Actual')

## The following two must be the same!
make.colour.plot(dat = cell.dat,
                 x.axis = "UMAP_X",
                 y.axis = "UMAP_Y",
                 col.axis = ColumnNames[30],
                 title = 'Align',
                 align.xy.by = cell.dat[1:10,],
                 align.col.by = cell.dat[1:10,])

make.colour.plot(dat = cell.dat,
                 x.axis = "UMAP_X",
                 y.axis = "UMAP_Y",
                 col.axis = ColumnNames[30],
                 title = 'HexDemo_align',
                 align.xy.by = cell.dat[1:10,],
                 align.col.by = cell.dat[1:10,],
                 hex = TRUE,
                 hex.bins = 100)


## The following must work!
make.multi.marker.plot(dat = cell.dat,
                       x.axis = "UMAP_X",
                       y.axis = "UMAP_Y",
                       density.plot = TRUE,
                       plot.by = c(CellularCols),
                       align.xy.by = cell.dat.sub,
                       align.col.by = cell.dat.sub,
                       figure.title = "Multi_Normal",
                       dot.size = 0.5,
                       save.each.plot = TRUE)

## The following must work! and look similar to the one above
make.multi.marker.plot(dat = cell.dat,
                       x.axis = "UMAP_X",
                       y.axis = "UMAP_Y",
                       density.plot = TRUE,
                       plot.by = c(CellularCols),
                       align.xy.by = cell.dat.sub,
                       align.col.by = cell.dat.sub,
                       figure.title = "Multi_hex",
                       dot.size = 0.5,
                       save.each.plot = TRUE,
                       hex = TRUE,
                       hex.bins = 100)

## The following must work!
make.multi.plot(dat = cell.dat,
                         type = "colour",
                         x.axis = "UMAP_X",
                         y.axis = "UMAP_Y",
                         col.axis = ColumnNames[30],
                         plot.by = group.col,
                         align.xy.by = cell.dat,
                         align.col.by = cell.dat,
                         dot.size = 0.5,
                         save.each.plot = TRUE)

## The following must work! and look similar to the one above
make.multi.plot(dat = cell.dat,
                type = "colour",
                x.axis = "UMAP_X",
                y.axis = "UMAP_Y",
                col.axis = ColumnNames[30],
                plot.by = group.col,
                align.xy.by = cell.dat,
                align.col.by = cell.dat,
                dot.size = 0.5,
                save.each.plot = TRUE,
                hex = TRUE,
                hex.bins = 50)

