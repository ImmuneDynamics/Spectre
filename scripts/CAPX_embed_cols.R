##########################################################################################################
#### DRAFT Cytometry Analysis Pipeline for large and compleX data (CAPX) v3.0 - using the Spectre R package
#### Part X - annotate clusters
##########################################################################################################


### Install/load packages



### Read in cell.dat and cell.dat.sub



### Perform other clustering approaches if desired

Spectre::embed.columns(x = cell.dat,
                       type = "data.frame",
                       base.name = "FlowSOM_metacluster",
                       col.name = "PopName",
                       match.to = c(1:40),
                       new.cols = c(rep("PMN", 20), rep("Tcell", 20)))

cell.dat <- cbind(cell.dat, embed.res)
cell.dat
rm(embed.res)

### Write files


### Sumtables


### Additional plots

