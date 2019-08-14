

setwd("/Users/Tom/Desktop/Output_CellNums/")
list.files(getwd(), ".csv")

fold.hmap <- read.csv("SumTable_CellsPerTissue_FoldChangeLog2.csv")


setwd("/Users/Tom/Desktop/Output_MFI_per_sample/")
list.files(getwd(), ".csv")

MFI.clustVmarker <- read.csv("SumTable_MFI_cluster_x_marker-all_data.csv")


hmap.foldcell <- fold.hmap
hmap.mfi <- MFI.clustVmarker

save(hmap.foldcell, "hmap.foldcell.Rdata")
save(hmap.mfi, "hmap.mfi.Rdata")
