
library('devtools')
devtools::install_github(repo = "sydneycytometry/spectre", ref = 'development')
library("Spectre")

setwd("/Users/Tom/Google Drive (t.ashhurst@centenary.org.au)/_Sydney Cytometry/2019_Synced/GitHub/Public github/Spectre - workflow scripts/")

getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

## Can set manually using these lines, if desired
#PrimaryDirectory <- "/Users/thomasashhurst/Documents/Github/Public Github/Spectre/Other/Demo_dataset/"


library(devtools)
library(roxygen2)










setwd("/Volumes/PRJ-sydneycytometry/Consultation projects/Hansbro/TA270.2 2019-08-06/TA270_BM_norm_Leukocytes/2019-08-16/Output_fileprep/Output_CAPX/")

list.files(getwd(), ".csv")

## Read data in
Spectre::read.files(file.loc = getwd(),
                    file.type = ".csv",
                    do.embed.file.names = FALSE)

names(data.list)

cell.dat.sub <- data.list[[1]]
#cell.dat.sub <- data.list[[2]]

names(cell.dat.sub)

colour.plot(d = cell.dat.sub,
            x.axis = "UMAP_42_X",
            y.axis = "UMAP_42_Y",
            col.axis = "UMAP_42_X",
            colour = "magma",
            title = "")

factor.plot(d = cell.dat.sub,
            x.axis = "UMAP_42_X",
            y.axis = "UMAP_42_Y",
            col.axis = "Batch",
            dot.size = 0.25,
            #colour = "magma",
            title = "")


