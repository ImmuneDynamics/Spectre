
Parameter <- c("ExperimentID", "ExperimentName", "Date")
Entry <- c("TA170", "BM analysis", "20190812")

dat <- list()

dat[["data"]] <- cell.dat
dat[["metadata"]] <- list()

dat$metadata[["exp.details"]] <- data.frame(Parameter,Entry)
dat$metadata[["sample.details"]] <- sample.table
dat$metadata[["analysis.log"]] <- data.frame("Function" = NA, "Parameter" = NA, "Value" = NA)

fwrite(x = dat[[1]], "dat.csv")
#write.csv(x = dat[[2]], "meta.csv")


#dat[[1]]
dat[[2]]

dat[[2]][[1]]
dat[[2]][[2]]
dat[[2]][[3]]










file.col <- "Filename"
sample.col <- "Sample"
group.col <- "Group"
batch.col <- "Batch"


all.sample.names <- sample.table[sample.col]
all.group.names <- sample.table[group.col]
all.batch.names <- sample.table[batch.col]
cells.per.sample <- sample.table["Cells.per.sample"]


all.file.names == unique(names(data.list))

#### USE SOME KIND OF LOOP BASED ON WHATS IN THE TABLE -- ADD KEYWORDS/NUMS to the samples

# make lists from what's in the table

for(i in c(1:length(all.file.names))){
  #i <- 1
  data.list[[i]][[sample.col]] <- NA # fills a new colum
  data.list[[i]][[sample.col]] <- all.sample.names[i,]

  data.list[[i]][[group.col]] <- NA # fills a new colum
  data.list[[i]][[group.col]] <- all.group.names[i,]

  data.list[[i]][[batch.col]] <- NA # fills a new colum
  data.list[[i]][[batch.col]] <- all.batch.names[i,]
}

head(data.list[[1]])
head(data.list[[6]])
head(data.list[[12]])




################
#dim(cell.dat)
#cell.dat.large <- rbind(cell.dat, cell.dat)
#for(i in c(1:8)){
#  cell.dat.large <- rbind(cell.dat.large, cell.dat) # x8 or so
#}
#cell.dat <- cell.dat.large
################
