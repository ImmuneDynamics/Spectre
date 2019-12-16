##########################################################################################################
#### 8. Print tSNE/UMAP plots to disk
##########################################################################################################

setwd(PrimaryDirectory)
setwd(OutputDirectory)
dir.create("Output-ColourPlots")
setwd("Output-ColourPlots")

### Loop for cellular markers etc
setwd(OutputDirectory)
setwd("Output-ColourPlots")

plots <- CellularCols

## Plot for all data
for(a in plots){
  p <- Spectre::colour.plot(d = cell.dat.sub,
                            x.axis = "UMAP_42_X",
                            y.axis = "UMAP_42_Y",
                            col.axis = a,
                            title = a,
                            colours = "magma",
                            dot.size = 1)
  
  ggsave(p, filename = paste0("All_samples_", a, ".png"), width = 9, height = 7)
}

## Plot divided by groups

## Plot divided by samples
## Plot by sample
#for(i in unique(cell.dat.sub[["FileName"]])){
#  # subset
#}
## Plot by group

### Loop for clusters etc

setwd(OutputDirectory)
setwd("Output-ColourPlots")

head(cell.dat.sub)
factors <- c("FlowSOM_metacluster")

## Plot for all data
for(a in factors){
  p <- Spectre::labelled.factor.plot(d = cell.dat.sub,
                                     x.axis = "UMAP_42_X",
                                     y.axis = "UMAP_42_Y",
                                     col.axis = "FlowSOM_metacluster",
                                     title = "Cluster",
                                     dot.size = 1) # assumes numeric
  
  ggsave(p, filename = paste0("All_samples_", a, ".png"), width = 9, height = 7)
}

setwd(PrimaryDirectory)

## Plot divided by groups


## Plot divided by samples
